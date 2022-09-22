import itertools
import os
import random
import time
from collections import defaultdict

import networkx as nx
import numpy as np
import prody
from scipy.spatial import distance

from utils import *


class Sampler:
    """Class to sample conformations, based on a"""

    def __init__(
        self,
        molecules,
        envrionment,
        dihe_parameters,
        vdw_parameters,
        clash_dist=1.8,
        grid_resolution=1.5,
    ):
        """
        Parameters
            molecules: list of Molecules instances
            environment: AtomGroup that will not be samples (e.g. protein, membrane, etc.)
            dihe_parameters: dictionary of parameters for dihedrals from CHARMMParameters. Atomtype as key and [k, n, d] as values
            vdw_parameters: dictionary of parameters for van der Waals from CHARMMParameters. Atomtype as key and [r, e] as values
            clash_dist = threshold for defining a clash (A)
        """
        self.molecules = molecules
        self.environment = envrionment
        self.clash_dist = clash_dist
        self.cutoff_dist = 10.0
        self.dihe_parameters = dihe_parameters
        self.vdw_parameters = vdw_parameters
        self.energy = {}
        self.energy_lookup = []
        self.nbr_clashes = np.zeros(len(self.molecules))
        self.non_bonded_energy = []
        self.charges = []
        self.vdw = []
        self.exclude1_3 = []
        self.genes = []
        self.sample = []
        # size of environment
        if self.environment:
            self.grid_resolution = grid_resolution
            c = self.environment.select("not resname ASN").getCoords()
            self.c_min = np.min(c, axis=0)
            self.c_max = np.max(c, axis=0)
            c_size = np.round((self.c_max - self.c_min) / self.grid_resolution).astype(
                int
            )
            self.bins = [
                np.linspace(self.c_min[0], self.c_max[0], c_size[0]),
                np.linspace(self.c_min[1], self.c_max[1], c_size[1]),
                np.linspace(self.c_min[2], self.c_max[2], c_size[2]),
            ]
            self.environment_grid = np.full(c_size, False)
            x_idx = np.digitize(c[:, 0], self.bins[0]) - 1
            y_idx = np.digitize(c[:, 1], self.bins[1]) - 1
            z_idx = np.digitize(c[:, 2], self.bins[2]) - 1
            self.environment_grid[(x_idx, y_idx, z_idx)] = True

        self.parse_patches(os.path.join(GLYCOSYLATOR_PATH, "support/topology/pres.top"))
        self.interresidue_torsionals = []

        self.coordinate_idx = np.zeros(len(molecules) + 1, dtype=int)
        self.exclude_nbr_clashes = np.zeros(len(molecules))
        idx = 0
        #        self.coordinate_idx = int(idx)
        for mol_id, molecule in enumerate(self.molecules):
            idx += molecule.atom_group.numAtoms()
            self.coordinate_idx[mol_id + 1] = idx
            self.sample.append(True)
            self.genes.append(None)
            types = nx.get_node_attributes(molecule.connectivity, "type")
            keys = types.keys()
            keys.sort()
            charges = nx.get_node_attributes(molecule.connectivity, "charge")
            names = nx.get_node_attributes(molecule.connectivity, "name")
            charge = []
            vdw = []
            for k in keys:
                charge.append(charges[k])
                e, r = vdw_parameters[types[k]][:2]
                # increase the r_min for hydrogen atoms
                if r < 1.0:
                    r = 1.340
                vdw.append([e, r])

            self.charges.append(charge)
            self.vdw.append(vdw)
            lookup = []

            self.build_1_3_exclude_list(mol_id)
            self.count_self_exclude(mol_id)
            # self.nbr_clashes.append(self.count_total_clashes(mol_id))
            # self.non_bonded_energy.append(self.compute_total_energy(mol_id))
            # self.count_total_clashes(mol_id)
            # self.compute_total_energy(mol_id)

            self.interresidue_torsionals.append(
                molecule.get_interresidue_torsionals(self.patches)
            )
            self.energy["skip"] = self.compute_inv_cum_sum_dihedral([[0.35, 1.0, 0.0]])
            for dihe in molecule.torsionals:
                atypes = []
                for d in dihe:
                    atypes.append(types[d])
                k1 = "-".join(atypes)
                atypes.reverse()
                k2 = "-".join(atypes)
                if k1 in self.energy:
                    lookup.append(k1)
                    continue
                if k2 in self.energy:
                    lookup.append(k2)
                    continue

                if k1 in dihe_parameters:
                    k = k1
                elif k2 in dihe_parameters:
                    k = k2
                else:
                    print("Missing parameters for ") + k1
                    print("This dihedral will be skipped")
                    lookup.append("skip")
                    continue
                par_list = dihe_parameters[k]
                self.energy[k] = self.compute_inv_cum_sum_dihedral(par_list)
                lookup.append(k)

            self.energy_lookup.append(lookup)

        self.molecule_coordinates = np.zeros((idx, 3))
        self.count_total_clashes_fast()
        # for i,molecule in enumerate(self.molecules):
        #    i0,i1 = self.coordinate_idx[i:i+2]
        #    self.molecule_coordinates[i0:i1, :] = molecule.atom_group.getCoords()
        # self.non_bonded_energy = np.array(self.non_bonded_energy)

    def compute_inv_cum_sum_dihedral(self, par_list, n_points=100):
        """Computes the an interpolation of the inverse transform sampling of a CHARMM
        dihedral term: sum(k*[1-cos(n*phi -d)]).
        This allows to map a uniform distribution [0:1[ to the dihedral energy distribution
        Parameters:
            par_list: list of parameters [k, n, d]
            n_points: number of points used for computing the energy function
        Returns:
            inv_cdf: interpolate object
        """
        phi = np.linspace(-180.0, 180.0, n_points)
        energy = np.zeros(phi.shape)
        for p in par_list:
            k, n, d = p
            energy += (0.01 + k) * (1 - np.cos((n * phi - d) * np.pi / 180.0))
        energy = np.max(energy) - energy
        energy = energy - np.min(energy)
        energy = energy / np.sum(energy[1:] * np.diff(phi))
        cum_values = np.zeros(energy.shape)
        cum_values[1:] = np.cumsum(energy[1:] * np.diff(phi))
        # plt.plot(cum_values, phi)
        #  x has to be strickly increasing
        idx = np.full(cum_values.size, True)
        idx[1:] = np.diff(cum_values) > 0.0
        x = cum_values[idx]
        y = phi[idx]
        inv_cdf = prody.InterpolatedUnivariateSpline(x, y)
        return inv_cdf

    def parse_patches(self, fname):
        lines = readLinesFromFile(fname)
        self.patches = defaultdict(list)
        for line in lines:
            line = line.split("\n")[0].split("!")[0].split()  # remove comments and endl
            if line:
                if line[0] == "PRES":
                    patch = line[1]
                if line[0] == "DIHE":
                    dihe = [line[1:5], list(itertools.pairwise(map(float, line[5:])))]
                    self.patches[patch].append(dihe)

    def get_uniform(self, interp_fn, angle):
        """Returns a number between [0:1[ which corresponds to an angle,
        based on the distribution of angle energy.
        Parameters:
            inter_fn: interpolate object
            angle: angle in degrees

        Return:
            root: number between [0:1[
        """
        if angle > 180:
            angle = 180 - angle
        x = np.linspace(0, 1, 10)
        y = [interp_fn(i) - angle for i in x]
        spl = prody.InterpolatedUnivariateSpline(x, y)
        r = spl.roots()
        if not r:
            return 1.0
        else:
            return r[0]

    def build_1_3_exclude_list(self, mol_id):
        """list with set of neighboring atoms"""
        molecule = self.molecules[mol_id]
        G = molecule.connectivity
        exclude_mol = []
        for a in sorted(G.nodes()):
            exclude = set()
            for n in G.neighbors(a):
                exclude.add(n)
                exclude.update(G.neighbors(n))
            exclude_mol.append(exclude)
        self.exclude1_3.append(exclude_mol)

    def count_self_exclude(self, mol_id):
        """Counts the number bonds and 1_3 exclusion for each molecule
        KDTree based
        Parameters:
            mol_id: id of molecule to consider
            increment: should the overwrite or update nbr_clashes
        Returns
            nbr_clashes: the number of clashes
        """
        molecule = self.molecules[mol_id]
        kd = prody.KDTree(molecule.atom_group.getCoords())
        kd.search(self.clash_dist)
        atoms = kd.getIndices()
        exclude_mol = self.exclude1_3[mol_id]
        c = 0
        for a1, a2 in atoms:
            if a1 + 1 in exclude_mol[a2]:
                c += 1
        self.exclude_nbr_clashes[mol_id] = c

    def count_total_clashes_fast(self):
        for i, molecule in enumerate(self.molecules):
            i0, i1 = self.coordinate_idx[i : i + 2]
            self.molecule_coordinates[i0:i1, :] = molecule.atom_group.getCoords()
        self.count_clashes_fast()
        self.count_environment_clashes_grid()

    def count_clashes_fast(self):
        """Counts all the clashes for molecules at ones. KDTree == (nbr_clashes + nbr_bonds).
        Since nbr_bonds is constant, it is not required to recount them each time.
        """
        kd = prody.KDTree(self.molecule_coordinates)
        kd.search(self.clash_dist)
        atoms = kd.getIndices().flatten()
        self.nbr_clashes = (
            np.histogram(atoms, self.coordinate_idx)[0] / 2 - self.exclude_nbr_clashes
        )

    def count_environment_clashes_grid(self):
        if self.environment:
            idx = np.all(
                np.logical_and(
                    self.molecule_coordinates >= self.c_min[np.newaxis, :],
                    self.molecule_coordinates <= self.c_max[np.newaxis, :],
                ),
                axis=1,
            )
            counts = np.zeros(idx.shape)
            if np.sum(idx) > 0:
                x_idx = np.digitize(self.molecule_coordinates[idx, 0], self.bins[0]) - 1
                y_idx = np.digitize(self.molecule_coordinates[idx, 1], self.bins[1]) - 1
                z_idx = np.digitize(self.molecule_coordinates[idx, 2], self.bins[2]) - 1
                counts[idx] = self.environment_grid[x_idx, y_idx, z_idx]
                self.nbr_clashes += np.histogram(
                    np.argwhere(counts), self.coordinate_idx
                )[0]

    def count_total_clashes(self, mol_id, increment=False):
        """Counts the total number of clashes in the system.
        Performed in two steps:
                - count the number of clashes within a molecules (KDTree)
                - count the number of clashes between the molecule and it's environment (Grid)
        Parameters:
            mol_id: id of the molecule
        """
        nbr_clashes = self.count_self_clashes(mol_id, increment=increment)
        nbr_clashes += self.count_environment_clashes(mol_id, increment=increment)
        if not increment:
            self.nbr_clashes[mol_id] = nbr_clashes
        return self.nbr_clashes[mol_id]

    def count_self_clashes(self, mol_id, increment=False):
        """Counts the number of clashes for a molecule
        KDTree based
        Parameters:
            mol_id: id of molecule to consider
            increment: should the overwrite or update nbr_clashes
        Returns
            nbr_clashes: the number of clashes
        """
        molecule = self.molecules[mol_id]
        kd = prody.KDTree(molecule.atom_group.getCoords())
        kd.search(self.clash_dist)
        atoms = kd.getIndices()
        G = molecule.connectivity
        nbr_clashes = 0
        exclude_mol = self.exclude1_3[mol_id]
        for a1, a2 in atoms:
            if a1 + 1 not in exclude_mol[a2]:
                nbr_clashes += 1
        if increment:
            self.nbr_clashes[mol_id] += nbr_clashes
        else:
            self.nbr_clashes[mol_id] = nbr_clashes
        return nbr_clashes

    def count_environment_clashes(self, mol_id, increment=False):
        """Counts the number of a molecule and its environment
        Parameters:
            mol_id: id of molecule
            increment: overwite or increment nbr_clashes. If increment assumes that clashes have already been computed for all molecules < mol_id
        """
        XA = self.molecules[mol_id].atom_group.getCoords()
        nbr_clashes = 0
        if self.environment:
            XB = self.environment.select("not resname ASN").getCoords()
            Y = distance.cdist(XA, XB, "euclidean")
            nbr_clashes += np.sum(Y < self.clash_dist)

        start = 0
        if increment:
            start = mol_id + 1

        for mol in np.arange(start, len(self.molecules)):
            if mol == mol_id:
                continue
            XB = self.molecules[mol].atom_group.getCoords()
            Y = distance.cdist(XA, XB, "euclidean")
            nbr = np.sum(Y < self.clash_dist)
            if increment:
                self.nbr_clashes[mol] += nbr
            nbr_clashes += nbr

        if increment:
            self.nbr_clashes[mol_id] += nbr_clashes
        else:
            self.nbr_clashes[mol_id] = nbr_clashes
        return nbr_clashes

    def compute_total_energy(self, mol_id, fast=False, threshold=1e5):
        """Computes the total non-bonded energy for a system.
        Performed in two steps:
                - computes the non-bonded energy within a molecules (KDTree)
                - computes the non-bonded energy  between the molecule and it's environment (Grid)
        Parameters:
            mol_id: id of the molecule
            threshold: energy below which the environment energy is also computed
        """
        energy_self = self.compute_self_non_bonded_energy(mol_id)
        energy = self.compute_environment_energy(mol_id)
        if energy_self > threshold and fast:
            energy = energy_self
        else:
            energy = self.compute_environment_energy(mol_id)
        return energy_self, energy

    def compute_self_non_bonded_energy(self, mol_id, repulsion=0.0):
        """Computes the van der Waals and electrostatic interaction of glycan
        The non bonded energy will be computed for all atoms, except those connected by a bond
        Parameters:
            mol_id: id of molecule
        Returns:
            energy: non bonded energy
        """
        molecule = self.molecules[mol_id]
        kd = prody.KDTree(molecule.atom_group.getCoords())
        kd.search(self.cutoff_dist)
        atoms = kd.getIndices()
        distances = kd.getDistances()
        G = molecule.connectivity
        energy = 0
        exclude_mol = self.exclude1_3[mol_id]
        for a, r in zip(atoms, distances):
            a1, a2 = a
            if a1 + 1 not in exclude_mol[a2]:
                e1, r1 = self.vdw[mol_id][a1]
                e2, r2 = self.vdw[mol_id][a2]
                rvdw = (r1 + r2 + repulsion) / r
                rvdw = rvdw**6
                # Lennard-Jones
                energy += np.sqrt(e1 * e2) * (rvdw**2 - 2 * rvdw)
                # Electrostatic
                # c1 = self.charges[mol_id][a1]
                # c2 = self.charges[mol_id][a2]
                # energy += c1*c2 / r
        return energy

    def compute_environment_energy(self, mol_id):
        """Computes the non bonded energy between a molecule and its environment. The environment is approximated by CA atom type
        Parameters:
            mol_id: id of molecule
        """
        XA = self.molecules[mol_id].atom_group.getCoords()
        atoms = self.molecules[mol_id].atom_group.getSerials() - 1
        nbr_clashes = 0
        r_env = 1.992400
        e_env = -0.070000

        energy = 0
        if self.environment:
            XB = self.environment.select("not resname ASN").getCoords()
            Y = distance.cdist(XA, XB, "euclidean")
            for a in atoms:
                e1, r1 = self.vdw[mol_id][a]
                e = np.sqrt(e1 * e_env)
                r = Y[a, :]
                rvdw = ((r1 + r_env) / r[r < self.cutoff_dist]) ** 6
                energy += e * np.sum(rvdw**2 - 2 * rvdw)

        for mol in np.arange(len(self.molecules)):
            if mol == mol_id:
                continue
            XB = self.molecules[mol].atom_group.getCoords()
            Y = distance.cdist(XA, XB, "euclidean")
            for a in atoms:
                e1, r1 = self.vdw[mol_id][a - 1]
                e = np.sqrt(e1 * e_env)
                r = Y[a, :]
                rvdw = ((r1 + r_env) / r[r < self.cutoff_dist]) ** 6
                energy += e * np.sum(rvdw**2 - 2 * rvdw)
        return energy

    def compute_dihedral_energy(self, mol_id):
        """Computes the CHARMM torsional angles energy for a molecule
        Parameter:
            mol_id: index of molecule
        Return:
            energy: torsional energy of the molecule
        """
        lookup = self.energy_lookup[mol_id]
        molecule = self.molecules[mol_id]
        angles = molecule.get_all_torsional_angles()
        energy = 0
        for phi, k in zip(angles, lookup):
            par_list = self.dihe_parameters[k]
            for p in par_list:
                k, n, d = p
                energy += k * (1 - np.cos((n * phi - d) * np.pi / 180.0))
        return energy

    def compute_TotalEnergy(self, torsionals):
        """Compute the total energy of a molecule; Torsional + Non-bonded
        Parameters:
            torsionals: list of all torsional angles
        """
        counter = 0
        mol_id = 0
        energy = 0
        for molecule in self.molecules:
            ncounter = counter + len(molecule.torsionals)
            molecule.set_torsional_angles(
                molecule.torsionals, torsionals[counter:ncounter]
            )
            energy += np.sum(self.compute_total_energy(mol_id))
            energy += self.compute_dihedral_energy(mol_id)
            counter = ncounter
            mol_id += 1
        return energy

    def _get_all_torsional_angles(self):
        """Measures and returns all the torsional angle in molecules"""
        angles = []
        n_torsionals = []
        for molecule in self.molecules:
            a = molecule.get_all_torsional_angles()
            n_torsionals.append(len(a))
            angles += a
        return angles, n_torsionals

    def _build_individue_from_angles(self, mol_ids=[]):
        """Generate an individue for GA based on the torsional angles present in molecules"""
        if not len(mol_ids):
            mol_ids = np.arange(len(self.molecules))
        i = 0
        # set torsional angles
        individue = []
        for mol_id in mol_ids:
            molecule = self.molecules[mol_id]
            torsionals = molecule.get_all_torsional_angles()
            thetas = []
            t_id = 0
            for t in torsionals:
                e = self.energy_lookup[mol_id][t_id]
                # if t < 0:
                #    t += 360.
                # individue.append(t/360.)
                individue.append(self.get_uniform(self.energy[e], t))
                t_id += 1
            mol_id += 1
        return individue

    def _eugenics(self, percentage_of_population=0.75, mol_ids=[]):
        """Uses a set of predefined angles to bias the sampled torsional angles. torsional angles defined in pres.top will be biased
        Parameters:
            percentage_of_population: percent of the total population that should be biased
        """
        if not len(mol_ids):
            mol_ids = np.arange(len(self.molecules))
        i = 0

        size = int(self.population.shape[0] * percentage_of_population)
        for mol_id in mol_ids:
            molecule = self.molecules[mol_id]
            interresidue = self.interresidue_torsionals[mol_id]
            n = len(molecule.torsionals)
            torsionals = self.population[:, i : i + n]
            for torsional_ids in interresidue.keys():
                p_id, deltas = interresidue[torsional_ids]
                patch = self.patches[p_id]
                n = len(patch[0][1])
                preferred_angles = np.random.randint(n, size=size)
                # preferred_angles = self.gmm[p_id].sample(size)
                for p, t_id in zip(patch, map(int, torsional_ids.split("-"))):
                    e = self.energy_lookup[mol_id][t_id]
                    angles = []
                    torsional = []
                    for p_angle in preferred_angles:
                        t1 = self.get_uniform(self.energy[e], p[1][p_angle][0])
                        t2 = self.get_uniform(self.energy[e], p[1][p_angle][1])
                        torsional.append(np.random.uniform(t1, t2))

                    torsionals[:size, t_id] = torsional

    def _build_individue(self, individue, mol_ids=[]):
        """Builds and sets torsionals angles of molecules, based on a individue
        Parameter:
            individue: list of corresponding to one individue in population
        """
        # set torsional angles
        if not len(mol_ids):
            mol_ids = np.arange(len(self.molecules))
        i = 0
        for mol_id in mol_ids:
            molecule = self.molecules[mol_id]
            n = len(molecule.torsionals)
            if self.sample[mol_id]:
                torsionals = individue[i : i + n]
                thetas = []
                for t_id, t in enumerate(torsionals):
                    e = self.energy_lookup[mol_id][t_id]
                    thetas.append(self.energy[e](t))
                    # thetas.append(t*360)
                molecule.set_torsional_angles(
                    molecule.torsionals, thetas, absolute=True
                )
            mol_id += 1
            i += n

    def _evaluate_population(self, clash=False, fast=False, mol_ids=[]):
        """Evaluates the fittnest of a population:
        Parameters:
            clash: boolean defing if the nonbonded energy(default) or clashes should be computed
            fast: boolean defining if the fast (and inacurate) implementation of clash/energy algorithms should be considered
            mol_ids: only consider subgroup of molecules with index
        """
        energies = []
        t_build = []
        t_energy = []
        for p in self.population:
            t_1 = time.time()
            self._build_individue(p, mol_ids)
            t_build.append(time.time() - t_1)
            # evaluate energy
            mol_id = 0
            energy = 0
            t_1 = time.time()
            self.count_total_clashes_fast()
            energies.append(np.sum(self.nbr_clashes))
            t_energy.append(time.time() - t_1)
        print("Average build time: ", np.mean(t_build))
        print("Average evaluation time: ", np.mean(t_energy))
        ee = np.argsort(energies)
        print(
            "Best energy: ",
            f"{energies[ee[0]]:e}",
            "|| Median energy: ",
            f"{np.median(energies):e}",
            "|| Worst energy: ",
            f"{energies[ee[-1]]:e}",
        )
        return ee

    def _immortalize_fittest(self, threshold=0.01):
        """Immortalizes individue with clashes/energy contributing to less than threshold (fraction) of the total energy/clashes
        Parameters:
            threshold: defines which individue should be immortalized
        """
        total_clashes = np.sum(self.nbr_clashes)
        highlanders = (self.nbr_clashes / total_clashes) < threshold
        fittest = self.population[0, :]
        i = 0
        mol_id = 0
        for molecule in self.molecules:
            n = len(molecule.torsionals)
            # immortalize new fit
            if highlanders[mol_id] and self.sample[mol_id]:
                self.sample[mol_id] = False
                self.genes[mol_id] = fittest[i : i + n]
                thetas = []
                t_id = 0
                for t in self.genes[mol_id]:
                    e = self.energy_lookup[mol_id][t_id]
                    thetas.append(self.energy[e](t))
                    t_id += 1
                molecule.set_torsional_angles(molecule.torsionals, thetas)
            # Make unfit immortal perishable
            elif not self.sample[mol_id] and not highlanders[mol_id]:
                fittest[i : i + n] = self.genes[mol_id]
                self.genes[mol_id] = None
                self.sample[mol_id] = True
            mol_id += 1
            i += n

    def _select_fit(self, sorted_population, nbr_survivors=0.3, lucky_few=0.2):
        """Selection of the fittest individue for mating
        Parameters:
            sorted_population: sorted indices of population fitness
            nbr_survivors: number of surviving individues
            luck_few: fraction of less fit individue that survive
        """
        n = int(len(sorted_population) * nbr_survivors)
        mates = np.empty([n, self.population.shape[1]])
        fittest = int(n * (1 - lucky_few))
        mates[0:fittest, :] = self.population[sorted_population[0:fittest], :]
        mates[fittest:, :] = random.sample(self.population, n - fittest)
        return mates

    def _create_new_population(self, mates, pop_size, mutation_rate, crossover_rate):
        """Generates a new populations based on a list of mates
        Parameters:

        """
        # keep best
        self.population[0, :] = mates[0, :]
        i = 1
        while i < pop_size:
            # randomly select a mate
            mate1 = self._tournament(mates)
            if random.random() < crossover_rate:
                mate2 = self._tournament(mates)
                offsprings = self._crossover(mate1, mate2)
            else:
                offsprings = [mate1.copy()]
            self._mutate(offsprings, mutation_rate)

            for offspring in offsprings:
                self.population[i, :] = offspring
                i += 1
                if i >= pop_size:
                    break

    def _tournament(self, mates, size=5, choose_best=0.9):
        competitors = np.random.randint(mates.shape[0], size=size)
        competitors = np.sort(competitors)
        if random.random() < choose_best:
            return mates[competitors[0], :]
        else:
            return mates[random.choice(competitors[1:]), :]

    def _crossover(self, mate1, mate2):
        n = len(mate1)
        offsprings = np.empty([2, n])
        i1, i2 = np.sort(np.random.randint(n, size=2))
        if random.random() < 0.75:
            offsprings[0, 0:i1] = mate1[0:i1]
            offsprings[0, i1:i2] = mate2[i1:i2]
            offsprings[0, i2:] = mate1[i2:]
            offsprings[1, 0:i1] = mate2[0:i1]
            offsprings[1, i1:i2] = mate1[i1:i2]
            offsprings[1, i2:] = mate2[i2:]
        else:
            i = random.choice([i1, i2])
            offsprings[0, 0:i] = mate1[0:i]
            offsprings[0, i:] = mate2[i:]
            offsprings[1, 0:i] = mate2[0:i]
            offsprings[1, i:] = mate1[i:]

        return offsprings

    def _mutate(self, offsprings, mutation_rate):
        for offspring in offsprings:
            idx = np.random.rand(len(offspring)) < mutation_rate
            if idx.any():
                offspring[idx] = np.squeeze(np.random.rand(np.sum(idx)))

    def remove_clashes_GA_iterative(
        self,
        n_iter=10,
        n_individues=5,
        n_generation=50,
        pop_size=40,
        mutation_rate=0.01,
        crossover_rate=0.9,
    ):
        fast = False
        clash = True
        n_individues = np.min((n_individues, len(self.molecules)))
        for iter_cnt in np.arange(n_iter):
            # Select individue with highest number of clashes
            print("Iteration", iter_cnt)
            cnt = 0
            selected_molecules = []
            for idx in np.argsort(self.nbr_clashes)[::-1]:
                if self.sample[idx] and self.nbr_clashes[idx] > 0.0:
                    cnt += 1
                    selected_molecules.append(idx)
                if cnt > n_individues:
                    break

            #            selected_molecules = np.sort(np.argsort(self.nbr_clashes)[-n_individues:])
            selected_molecules = np.sort(selected_molecules)
            print("Selected Molecules", selected_molecules)
            torsionals = []
            n_torsionals = []
            for mol_id in selected_molecules:
                a = self.molecules[mol_id].get_all_torsional_angles()
                n_torsionals.append(len(a))
                torsionals += a

            length = len(torsionals)
            self.population = np.random.rand(pop_size, length)
            # self._eugenics(mol_ids = selected_molecules)
            # set input structure to first structure
            self.population[0, :] = self._build_individue_from_angles(
                mol_ids=selected_molecules
            )
            gen_cnt = 0
            while gen_cnt < n_generation:
                print("Generation:", gen_cnt)
                t1 = time.time()
                sorted_population = self._evaluate_population(
                    clash=clash, fast=fast, mol_ids=selected_molecules
                )
                t2 = time.time()
                print("Evaluation time: ", t2 - t1)

                t1 = time.time()
                mates = self._select_fit(sorted_population)
                t2 = time.time()
                print("Selection time: ", t2 - t1)
                t1 = time.time()
                self._create_new_population(
                    mates,
                    pop_size,
                    mutation_rate=mutation_rate,
                    crossover_rate=crossover_rate,
                )
                t2 = time.time()
                print("New population time: ", t2 - t1)
                gen_cnt += 1
                print("=" * 70)
            sorted_population = self._evaluate_population(
                clash=clash, fast=fast, mol_ids=selected_molecules
            )
            self._build_individue(
                self.population[sorted_population[0]], mol_ids=selected_molecules
            )

        sorted_population = self._evaluate_population(
            clash=clash, fast=fast, mol_ids=selected_molecules
        )
        self._build_individue(
            self.population[sorted_population[0]], mol_ids=selected_molecules
        )

    def remove_clashes_GA(
        self, n_generation=50, pop_size=40, mutation_rate=0.01, crossover_rate=0.9
    ):
        torsionals, n_torsionals = self._get_all_torsional_angles()
        length = len(torsionals)
        self.population = np.random.rand(pop_size, length)
        # self._eugenics()
        # set input structure to first structure
        self.population[0, :] = self._build_individue_from_angles()
        i = 0
        fast = False
        clash = True
        while i < n_generation:
            print("Generation:", i)
            t1 = time.time()
            sorted_population = self._evaluate_population(clash=clash, fast=fast)
            t2 = time.time()
            print("Evaluation time: ", t2 - t1)

            t1 = time.time()
            mates = self._select_fit(sorted_population)
            t2 = time.time()
            if i and not i % 10 and False:
                clash = not clash
                if clash:
                    pop_size = int(pop_size * 3)
                    self.population = np.zeros([pop_size, length])
                else:
                    pop_size = pop_size / 3
                    self.population = np.zeros([pop_size, length])

            if False and i == n_generation / 2:
                print("Evaluating self energy and reducing the size of the population")
                clash = False
                fast = True
                pop_size = pop_size / 3
                mutation_rate = 0.01
                self.population = np.zeros([pop_size, length])
            elif False and i == 3 * n_generation / 4:
                print("Evaluation full energy")
                fast = False

            print("Selection time: ", t2 - t1)
            t1 = time.time()
            self._create_new_population(
                mates,
                pop_size,
                mutation_rate=mutation_rate,
                crossover_rate=crossover_rate,
            )
            t2 = time.time()
            print("New population time: ", t2 - t1)
            i += 1
            print("=" * 70)
        sorted_population = self._evaluate_population(clash=clash, fast=fast)
        self._build_individue(self.population[sorted_population[0]])


#####################################################################################
#                               PSO Sampler                                         #
#####################################################################################
class SamplerPSO:
    """Class to sample conformations, based on a PSO optimization"""

    def __init__(
        self,
        molecules,
        envrionment,
        dihe_parameters,
        vdw_parameters,
        clash_dist=1.8,
        grid_resolution=1.5,
    ):
        """
        Parameters
            molecules: list of Molecules instances
            environment: AtomGroup that will not be samples (e.g. protein, membrane, etc.)
            dihe_parameters: dictionary of parameters for dihedrals from CHARMMParameters. Atomtype as key and [k, n, d] as values
            vdw_parameters: dictionary of parameters for van der Waals from CHARMMParameters. Atomtype as key and [r, e] as values
            clash_dist = threshold for defining a clash (A)
        """
        self.molecules = molecules
        self.environment = envrionment
        self.clash_dist = clash_dist
        self.cutoff_dist = 10.0
        self.dihe_parameters = dihe_parameters
        self.energy = {}
        self.energy_lookup = []
        self.nbr_clashes = np.zeros(len(self.molecules))
        self.exclude1_3 = []
        self.swarm = []
        # size of environment
        if self.environment:
            self.grid_resolution = grid_resolution
            c = self.environment.select("not resname ASN").getCoords()
            self.c_min = np.min(c, axis=0)
            self.c_max = np.max(c, axis=0)
            c_size = np.round((self.c_max - self.c_min) / self.grid_resolution).astype(
                int
            )
            self.bins = [
                np.linspace(self.c_min[0], self.c_max[0], c_size[0]),
                np.linspace(self.c_min[1], self.c_max[1], c_size[1]),
                np.linspace(self.c_min[2], self.c_max[2], c_size[2]),
            ]
            self.environment_grid = np.full(c_size, False)
            x_idx = np.digitize(c[:, 0], self.bins[0]) - 1
            y_idx = np.digitize(c[:, 1], self.bins[1]) - 1
            z_idx = np.digitize(c[:, 2], self.bins[2]) - 1
            self.environment_grid[(x_idx, y_idx, z_idx)] = True

        self.parse_patches(os.path.join(GLYCOSYLATOR_PATH, "support/topology/pres.top"))
        self.interresidue_torsionals = []

        self.coordinate_idx = np.zeros(len(molecules) + 1, dtype=int)
        self.exclude_nbr_clashes = np.zeros(len(molecules))
        idx = 0
        #        self.coordinate_idx = int(idx)
        for mol_id, molecule in enumerate(self.molecules):
            idx += molecule.atom_group.numAtoms()
            self.coordinate_idx[mol_id + 1] = idx
            types = nx.get_node_attributes(molecule.connectivity, "type")
            keys = types.keys()
            keys.sort()

            self.build_1_3_exclude_list(mol_id)
            self.count_self_exclude(mol_id)

            self.interresidue_torsionals.append(
                molecule.get_interresidue_torsionals(self.patches)
            )
            self.energy["skip"] = self.compute_inv_cum_sum_dihedral([[0.35, 1.0, 0.0]])
            lookup = []
            for dihe in molecule.torsionals:
                atypes = []
                for d in dihe:
                    atypes.append(types[d])
                k1 = "-".join(atypes)
                atypes.reverse()
                k2 = "-".join(atypes)
                if k1 in self.energy:
                    lookup.append(k1)
                    continue
                if k2 in self.energy:
                    lookup.append(k2)
                    continue

                if k1 in dihe_parameters:
                    k = k1
                elif k2 in dihe_parameters:
                    k = k2
                else:
                    print(f"Missing parameters for {k1}")
                    print("This dihedral will be skipped")
                    lookup.append("skip")
                    continue
                par_list = dihe_parameters[k]
                self.energy[k] = self.compute_inv_cum_sum_dihedral(par_list)
                lookup.append(k)

            self.energy_lookup.append(lookup)

        self.molecule_coordinates = np.zeros((idx, 3))
        self.count_total_clashes_fast()

    def compute_inv_cum_sum_dihedral(self, par_list, n_points=100):
        """Computes the an interpolation of the inverse transform sampling of a CHARMM
        dihedral term: sum(k*[1-cos(n*phi -d)]).
        This allows to map a uniform distribution [0:1[ to the dihedral energy distribution
        Parameters:
            par_list: list of parameters [k, n, d]
            n_points: number of points used for computing the energy function
        Returns:
            inv_cdf: interpolate object
        """
        phi = np.linspace(-180.0, 180.0, n_points)
        energy = np.zeros(phi.shape)
        for p in par_list:
            k, n, d = p
            energy += (0.01 + k) * (1 - np.cos((n * phi - d) * np.pi / 180.0))
        energy = np.max(energy) - energy
        energy = energy - np.min(energy)
        energy = energy / np.sum(energy[1:] * np.diff(phi))
        cum_values = np.zeros(energy.shape)
        cum_values[1:] = np.cumsum(energy[1:] * np.diff(phi))

        # plt.plot(cum_values, phi)
        idx = np.full(cum_values.size, True)
        idx[1:] = np.diff(cum_values) > 0.0
        x = cum_values[idx]
        y = phi[idx]

        inv_cdf = prody.InterpolatedUnivariateSpline(x, y)
        return inv_cdf

    def parse_patches(self, fname):
        lines = readLinesFromFile(fname)
        self.patches = defaultdict(list)
        for line in lines:
            line = line.split("\n")[0].split("!")[0].split()  # remove comments and endl
            if line:
                if line[0] == "PRES":
                    patch = line[1]
                if line[0] == "DIHE":
                    dihe = [line[1:5], list(itertools.pairwise(map(float, line[5:])))]
                    self.patches[patch].append(dihe)

    def get_uniform(self, interp_fn, angle):
        """Returns a number between [0:1[ which corresponds to an angle,
        based on the distribution of angle energy.
        Parameters:
            inter_fn: interpolate object
            angle: angle in degrees

        Return:
            root: number between [0:1[
        """
        if angle > 180:
            angle = 180 - angle
        x = np.linspace(0, 1, 10)
        y = [interp_fn(i) - angle for i in x]
        spl = prody.InterpolatedUnivariateSpline(x, y)
        r = spl.roots()
        if not r:
            return 1.0
        else:
            return r[0]

    def build_1_3_exclude_list(self, mol_id):
        """list with set of neighboring atoms"""
        molecule = self.molecules[mol_id]
        G = molecule.connectivity
        exclude_mol = []
        for a in sorted(G.nodes()):
            exclude = set()
            for n in G.neighbors(a):
                exclude.add(n)
                exclude.update(G.neighbors(n))
            exclude_mol.append(exclude)
        self.exclude1_3.append(exclude_mol)

    def count_self_exclude(self, mol_id):
        """Counts the number bonds and 1_3 exclusion for each molecule
        KDTree based
        Parameters:
            mol_id: id of molecule to consider
            increment: should the overwrite or update nbr_clashes
        Returns
            nbr_clashes: the number of clashes
        """
        molecule = self.molecules[mol_id]
        kd = prody.KDTree(molecule.atom_group.getCoords())
        kd.search(self.clash_dist)
        atoms = kd.getIndices()
        exclude_mol = self.exclude1_3[mol_id]
        c = 0
        for a1, a2 in atoms:
            if a1 + 1 in exclude_mol[a2]:
                c += 1
        self.exclude_nbr_clashes[mol_id] = c

    def count_total_clashes_fast(self):
        for i, molecule in enumerate(self.molecules):
            i0, i1 = self.coordinate_idx[i : i + 2]
            self.molecule_coordinates[i0:i1, :] = molecule.atom_group.getCoords()
        self.count_clashes_fast()
        self.count_environment_clashes_grid()

    def count_clashes_fast(self):
        """Counts all the clashes for molecules at ones. KDTree == (nbr_clashes + nbr_bonds).
        Since nbr_bonds is constant, it is not required to recount them each time.
        """
        kd = prody.KDTree(self.molecule_coordinates)
        kd.search(self.clash_dist)
        atoms = kd.getIndices().flatten()
        self.nbr_clashes = (
            np.histogram(atoms, self.coordinate_idx)[0] / 2 - self.exclude_nbr_clashes
        )

    def count_environment_clashes_grid(self):
        if self.environment:
            idx = np.all(
                np.logical_and(
                    self.molecule_coordinates >= self.c_min[np.newaxis, :],
                    self.molecule_coordinates <= self.c_max[np.newaxis, :],
                ),
                axis=1,
            )
            counts = np.zeros(idx.shape)
            if np.sum(idx) > 0:
                x_idx = np.digitize(self.molecule_coordinates[idx, 0], self.bins[0]) - 1
                y_idx = np.digitize(self.molecule_coordinates[idx, 1], self.bins[1]) - 1
                z_idx = np.digitize(self.molecule_coordinates[idx, 2], self.bins[2]) - 1
                counts[idx] = self.environment_grid[x_idx, y_idx, z_idx]
                self.nbr_clashes += np.histogram(
                    np.argwhere(counts), self.coordinate_idx
                )[0]

    def count_self_clashes(self, mol_id):
        """Counts the number of clashes for a molecule
        KDTree based
        Parameters:
            mol_id: id of molecule to consider
            increment: should the overwrite or update nbr_clashes
        Returns
            nbr_clashes: the number of clashes
        """
        molecule = self.molecules[mol_id]
        kd = prody.KDTree(molecule.atom_group.getCoords())
        kd.search(self.clash_dist)
        atoms = kd.getIndices()
        distances = kd.getDistances()
        G = molecule.connectivity
        nbr_clashes = 0
        exclude_mol = self.exclude1_3[mol_id]
        for ((a1, a2), d) in zip(atoms, distances):
            if a1 + 1 not in exclude_mol[a2]:
                nbr_clashes += 1.0
        return nbr_clashes

    def _get_all_torsional_angles(self):
        """Measures and returns all the torsional angle in molecules"""
        angles = []
        n_torsionals = []
        for molecule in self.molecules:
            a = molecule.get_all_torsional_angles()
            n_torsionals.append(len(a))
            angles += a
        return angles, n_torsionals

    def _build_position_from_angles(self, mol_ids=[]):
        """Generate the position of a particle based on the torsional angles present in molecules"""
        if not len(mol_ids):
            mol_ids = np.arange(len(self.molecules))
        i = 0
        # set torsional angles
        position = []
        for mol_id in mol_ids:
            molecule = self.molecules[mol_id]
            torsionals = molecule.get_all_torsional_angles()
            thetas = []
            t_id = 0
            for t in torsionals:
                e = self.energy_lookup[mol_id][t_id]
                position.append(self.get_uniform(self.energy[e], t))
                t_id += 1
            mol_id += 1
        return position

    def _build_molecule(self, position, mol_ids=[]):
        """Builds and sets torsionals angles of molecules, based on the position of a particle
        Parameter:
            position: list of corresponding to the position of a particle
        """
        # set torsional angles
        if not len(mol_ids):
            mol_ids = np.arange(len(self.molecules))
        i = 0
        for mol_id in mol_ids:
            molecule = self.molecules[mol_id]
            n = len(molecule.torsionals)
            torsionals = position[i : i + n]
            thetas = []
            for t_id, t in enumerate(torsionals):
                e = self.energy_lookup[mol_id][t_id]
                thetas.append(self.energy[e](t))
            molecule.set_torsional_angles(molecule.torsionals, thetas, absolute=True)
            i += n

    def _evaluate_swarm(self, mol_ids=[]):
        """Evaluates the fittnest of the swar,:
        Parameters:
            mol_ids: only consider subgroup of molecules with index
        """
        energies = []
        t_build = []
        t_energy = []
        for p in self.swarm:
            t_1 = time.time()
            self._build_molecule(p.position, mol_ids)
            t_build.append(time.time() - t_1)
            t_1 = time.time()
            self.count_total_clashes_fast()
            p.update_energy(np.sum(self.nbr_clashes))
            energies.append(p.energy)
            t_energy.append(time.time() - t_1)
        print("Average build time: ", np.mean(t_build))
        print("Average evaluation time: ", np.mean(t_energy))
        ee = np.argsort(energies)
        print(
            "Best energy: ",
            f"{energies[ee[0]]:e}",
            "|| Median energy: ",
            f"{np.median(energies):e}",
            "|| Worst energy: ",
            f"{energies[ee[-1]]:e}",
        )
        return ee

    class Particle:
        """This class implements the functions needed for taking care of the Particles
        Parameters:
            position: position of the particle [0;1]
            velocity: velocity of the particle [-1;1]
            pos_best: position of the best energy
            lowest_energy: lowest energy of the particle
            energy: current energy of the particle
        """

        def __init__(self, x0, inertia=0.75, cognitive_prm=1.75, social_prm=2.0):
            self.position = x0  # particle position
            self.velocity = 2 * np.random.rand(x0.shape[0]) - 1  # particle velocity
            self.pos_best = []  # best position individual
            self.lowest_energy = np.Inf  # lowest energy individual
            self.energy = np.Inf  # current energy individual
            self.w = inertia  # inertia constant
            self.c1 = cognitive_prm  # cognitive parameter
            self.c2 = social_prm  # social parameter

        def update_energy(self, energy):
            self.energy = energy
            # check to see if the current position is an individual best
            if self.energy < self.lowest_energy:
                self.pos_best = self.position
                self.lowest_energy = self.energy

        # update new particle velocity
        def update_velocity(self, pos_best_global):
            r1, r2 = np.random.rand(2)
            vel_cognitive = self.c1 * r1 * (self.pos_best - self.position)
            vel_social = self.c2 * r2 * (pos_best_global - self.position)
            self.velocity = self.w * self.velocity + vel_cognitive + vel_social

        # update the particle position
        def update_position(self, pos_best_global):
            self.update_velocity(pos_best_global)
            self.position = self.position + self.velocity
            # check bounds and reflect particle
            idx = self.position >= 1.0
            self.position[idx] = 1.0
            self.velocity[idx] = -0.1 * self.velocity[idx]
            idx = self.position <= 0
            self.position[idx] = 0
            self.velocity[idx] = -0.1 * self.velocity[idx]

    def remove_clashes_PSO(
        self,
        n_generation,
        n_molecules,
        n_particles,
        n_iter,
        inertia=0.75,
        cognitive_prm=1.5,
        social_prm=2.0,
        save_trajectory=False,
    ):

        n_molecules = np.min((n_molecules, len(self.molecules)))

        if save_trajectory:
            natoms = 0
            for i, m in enumerate(self.molecules):
                if i:
                    molecule_trajectory += m.atom_group
                else:
                    molecule_trajectory = m.atom_group
                natoms += m.atom_group.numAtoms()

        for gen_cnt in np.arange(n_generation):
            # Select molecules with highest number of clashes
            print("Generation", gen_cnt)
            cnt = 0
            selected_molecules = []
            for idx in np.argsort(self.nbr_clashes)[::-1]:
                if self.nbr_clashes[idx] > 0.0:
                    cnt += 1
                    selected_molecules.append(idx)
                if cnt > n_molecules:
                    break

            selected_molecules = np.sort(selected_molecules)
            print("Selected Molecules", selected_molecules)
            torsionals = []
            n_torsionals = []
            for mol_id in selected_molecules:
                a = self.molecules[mol_id].get_all_torsional_angles()
                n_torsionals.append(len(a))
                torsionals += a

            length = len(torsionals)
            # Build the swarm
            self.swarm = []
            positions = np.random.rand(n_particles, length)
            positions[0, :] = self._build_position_from_angles(
                mol_ids=selected_molecules
            )
            for x0 in positions:
                self.swarm.append(self.Particle(x0, inertia, cognitive_prm, social_prm))

            pos_best_global = None
            lowest_energy_global = np.Inf
            # set input structure to first structure
            iter_cnt = 0

            while iter_cnt < n_iter:
                print("Iteration:", iter_cnt)
                t1 = time.time()
                sorted_population = self._evaluate_swarm(mol_ids=selected_molecules)
                t2 = time.time()
                print("Total time for one iteration: ", t2 - t1)
                best_particle = self.swarm[sorted_population[0]]
                if best_particle.energy < lowest_energy_global:
                    pos_best_global = best_particle.position
                    lowest_energy_global = best_particle.energy
                    if save_trajectory:
                        self._build_molecule(
                            pos_best_global, mol_ids=selected_molecules
                        )
                        coords = np.zeros((natoms, 3))
                        i = 0
                        j = 0
                        for m in self.molecules:
                            j = i + m.atom_group.numAtoms()
                            coords[i:j, :] = m.atom_group.getCoords()
                            i = j
                        molecule_trajectory.addCoordset(coords)

                # update velocities and position
                max_vel = 0
                for particle in self.swarm:
                    particle.update_position(pos_best_global)
                iter_cnt += 1
                print("=" * 70)
            self._build_molecule(pos_best_global, mol_ids=selected_molecules)

        print("Best energy", lowest_energy_global)
        if save_trajectory:
            prody.writePDB("PSO_trajectory.pdb", molecule_trajectory)
