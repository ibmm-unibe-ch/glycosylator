import itertools

import networkx as nx
import prody

from .charmm_parameters import CHARMMParameters
from .charmm_topology import CHARMMTopology
from .utils import *


class MoleculeBuilder:
    """Class for building/modifying molecule"""

    def __init__(self, topofile, paramfile, force_field="charmm"):
        """
        Parameters:
            topofile: path to topology file
            paramfile: path to parameters file
            force_field: force field name. Currently only CHARMM
        """
        if force_field == "charmm":
            self.Topology = CHARMMTopology(topofile)
            self.Parameters = CHARMMParameters(paramfile)
        else:
            raise NotImplementedError(
                "Non-CHARMM force fields have not been implemented yet"
            )

    def init_new_residue(self, resnum, resname, chain, segname, i=1):
        """Initializes a residue from scratch
        Parameters:
            resid: residue id (int)
            resname: name of residue (str)
            chain: chain id (char)
            segn: segname (str)
        Returns:
            residue: AtomGroup with all atoms initialized (from Topology)
            atoms_name: name of of all the atoms in residue
            bonds: list of bonds (segn,chid,resi,atom_name)
        """
        residue = prody.AtomGroup(f"{resname}{resnum}")
        # add all the atoms
        atoms = self.Topology.get_atoms(resname)
        coords = np.zeros((len(atoms), 3))
        residue.setCoords(coords)

        # each entry in atoms is of form ['C1', 'CC3162', 0.34]
        atom_names = [atom[0] for atom in atoms]
        resnames = [resname] * len(atoms)
        resnums = [resnum] * len(atoms)
        chid = [chain] * len(atoms)
        segnames = [segname] * len(atoms)
        occupancies = [1] * len(atoms)
        betas = [0] * len(atoms)
        serials = list(range(i, len(atoms) + i))
        elements = [atom[0][0] for atom in atoms]
        icodes = [""] * len(atoms)
        altlocs = [""] * len(atoms)

        residue.setNames(atom_names)
        residue.setResnums(resnums)
        residue.setResnames(resnames)
        residue.setChids(chid)
        residue.setSegnames(segnames)
        residue.setOccupancies(occupancies)
        residue.setBetas(betas)
        residue.setSerials(serials)
        residue.setIcodes(icodes)
        residue.setElements(elements)
        residue.setAltlocs(altlocs)
        top_bonds = self.Topology.get_bonds(resname)

        #####
        id_r = f"{segname},{chain},{int(resnum)},,"
        bonds = [(id_r + a1, id_r + a2) for a1, a2 in itertools.pairwise(top_bonds)]

        return residue, atom_names, bonds

    def copy_atom(self, src_atom, dst_atom):
        """copies all the attributes of one atom to another
        Parameters:
            src_atom: original atom
            dst_atom: copy atom
        """
        dst_atom.setCoords(src_atom.getCoords())
        dst_atom.setNames(src_atom.getNames())
        # dst_atom.setResnums(src_atom.getResnums())
        dst_atom.setResnames(src_atom.getResnames())
        # dst_atom.setChids(src_atom.getChids())
        # dst_atom.setSegnames(src_atom.getSegnames())
        dst_atom.setOccupancies(src_atom.getOccupancies())
        dst_atom.setBetas(src_atom.getBetas())
        dst_atom.setSerials(src_atom.getSerials())
        dst_atom.setIcodes(src_atom.getIcodes())
        dst_atom.setElements(src_atom.getElements())
        dst_atom.setAltlocs(src_atom.getAltlocs())

    def find_missing_atoms(self, residue, resnum=None, chain=None, segname=None):
        """Add all missing atoms to a ProDy residue from topology
        Parameters:
            residue: ProDy residue (AtomGroup)
        Retruns:
            complete_residue: Completed ProDy residue, with coordinates of missing atoms set to (0, 0, 0)
            missing_atoms: list of missing atom names
            bonds: list of new bonds
        """
        if not resnum:
            resnum = residue.getResnums()[0]
        if not chain:
            chain = residue.getChids()[0]
        if not segname:
            segname = residue.getSegnames()[0]
        complete_residue, atom_names, bonds = self.init_new_residue(
            resnum, residue.getResnames()[0], chain, segname
        )
        missing_atoms = []
        res_atom_names = residue.getNames()
        for name in atom_names:
            if name in res_atom_names:
                atom = residue.select(f"name {name}")
                catom = complete_residue.select(f"name {name}")
                self.copy_atom(atom, catom)
            else:
                missing_atoms.append(name)
        complete_residue.setResnums([resnum] * len(complete_residue))
        return complete_residue, missing_atoms, bonds

    def apply_patch(self, patch, residue1, residue2):
        """
        Parameters:
            patch: name of patch (str)
            residue1: first residue in patch (ProDy AtomGroup)
            residue2: second residue in patch (ProDy AtomGroup)
        Return:
            bonds: list of all new bonds
            dele_atoms: list of atoms which should be deleted
        """
        bonds = []
        dele_atoms = []
        if patch:
            dele_atoms = self.dele_atoms(patch, residue1, residue2)
            bonds = self.patch_bonds(patch, residue1, residue2)
        return dele_atoms, bonds

    def build_IC_graph(self, atoms, ics):
        """Extracts ICs to build missing atoms
        Parameters:
            atoms: list of missing atoms
            ics: list of internal coordinates
        Returns:
            unsorted_graph: dictionay representation of dependency of ics
                            key: atom name
                            value: list of atoms connected to key
            required_atoms: list atoms required to build missing atoms
        """
        unsorted_graph = {}
        required_atoms = []
        atomsIC = [atom.replace("*", "") for ic in ics for atom in ic[0:4]]
        # Build graph
        for a in atoms:
            if a in atomsIC:
                required_atoms.append(a)
                if a not in unsorted_graph:
                    unsorted_graph[a] = []
                for ic in self.Topology.get_IC(ics, a):
                    for aic in ic[0:3]:
                        aic = aic.replace("*", "")
                        if aic in unsorted_graph:
                            unsorted_graph[aic].append(a)
                        else:
                            unsorted_graph[aic] = [a]
        return unsorted_graph, required_atoms

    def build_from_patch(self, link_residue, resid, resname, chain, segname, patch):
        """Build residue from a patch
        Parameters:
            link_residue: residue that the patch will use to build new residue
            resid: residue number
            resname: residue name
            chain: residue chain
            segname: residue segname
            patch: name of patch (str)
        Returns:
            denovo_residue: complete residue (AtomGroup)
            dele_atoms: list of atoms which should be deleted
            bonds: list of all new bonds
        """
        denovo_residue, missing_atoms, bonds = self.init_new_residue(
            resid, resname, chain, segname
        )
        ics = self.Topology.patches[patch]["IC"]
        # patch_atoms = sorted(set([atom.replace('*', '')[1:] for ic in ics for atom in ic[0:4] if atom.replace('*', '')[0]=='2']))
        patch_atoms = sorted(
            {
                atom.replace("*", "")
                for ic in ics
                for atom in ic[0:4]
                if atom.replace("*", "")[0] == "2"
            }
        )
        self.build_patch_missing_atom_coord(
            link_residue, denovo_residue, patch_atoms, ics
        )
        missing_atoms = [a for a in missing_atoms if f"2{a}" not in patch_atoms]
        ics = self.Topology.topology[resname]["IC"]
        self.build_missing_atom_coord(denovo_residue, missing_atoms, ics)

        dele_atoms, b = self.apply_patch(patch, link_residue, denovo_residue)
        bonds.extend(b)
        return denovo_residue, dele_atoms, bonds

    def patch_bonds(self, patch, residue1, residue2=None):
        """
        Parameters:
            patch: name of patch (str)
            residue1: first residue in patch
            residue2: second residue in patch. None if not required
        Returns:
            bonds: list of bonds
        """
        bonds = []
        segn1 = residue1.getSegnames()[0]
        chid1 = residue1.getChids()[0]
        resi1 = str(residue1.getResnums()[0])
        ic1 = str(residue1.getIcodes()[0])
        if residue2:
            segn2 = residue2.getSegnames()[0]
            chid2 = residue2.getChids()[0]
            resi2 = str(residue2.getResnums()[0])
            ic2 = str(residue2.getIcodes()[0])
        for a1, a2 in itertools.pairwise(self.Topology.patches[patch]["BOND"]):
            b = []
            for a in [a1, a2]:
                if a[0] == "1":
                    b.append(f"{segn1},{chid1},{resi1},{ic1},{a[1:]}")
                if a[0] == "2":
                    if residue2:
                        b.append(f"{segn2},{chid2},{resi2},{ic2},{a[1:]}")
                    else:
                        print(
                            f"Warning BOND: missing residue2 required for patch {patch}"
                        )
            bonds.append(b)
        return bonds

    def dele_atoms(self, patch, residue1, residue2=None):
        """
        Parameters:
            patch: name of patch (str)
            residue1: first residue in patch
            residue2: second residue in patch. None if not required
        Returns:
            dele_atoms: list of atoms to delete. (segn, chid, resi, atom_name)
        """
        dele_atoms = []
        if patch:
            atoms = self.Topology.patches[patch]["dele"]
            segn1 = residue1.getSegnames()[0]
            chid1 = residue1.getChids()[0]
            resi1 = str(residue1.getResnums()[0])
            ic1 = str(residue1.getIcodes()[0])
            if residue2:
                segn2 = residue2.getSegnames()[0]
                chid2 = residue2.getChids()[0]
                resi2 = str(residue2.getResnums()[0])
                ic2 = str(residue2.getIcodes()[0])

            for a in atoms:
                if a[0] == "1":
                    dele_atoms.append(f"{segn1},{chid1},{resi1},{ic1},{a[1:]}")
                if a[0] == "2":
                    if residue2:
                        dele_atoms.append(f"{segn2},{chid2},{resi2},{ic2},{a[1:]}")
                    else:
                        print(f"Warning: missing residue2 required for patch {patch}")
        return dele_atoms

    def delete_atoms(self, molecule, dele_atoms):
        """Removes all atoms that should be delete by patch
        Parameters:
            molecule: AtomGroup defining the molecule
            dele_atoms: list of atoms to be deleted. (segn, chid, resi, atom_name)
        """
        sel = []
        prefix = ["segment", "chain", "resid", "icode", "name"]
        for a in dele_atoms:
            s1 = []
            for p, s in zip(prefix, a.split(",")):
                if s:
                    s1.append(f"{p} {s}")

            sel.append(" and ".join(s1))
        sel = f"not ({' or '.join(sel)})"

        return molecule.select(sel).copy()

    def build_from_DUMMY(
        self,
        resid,
        resname,
        chain,
        segname,
        dummy_patch,
        dummy_coords=[[0, 0, 0], [0, 0, 1], [0, 1, 1]],
    ):
        """Builds residue from DUMMY atoms
        Parameters:
            resid: residue id (int)
            chain: residue chain id (chr)
            segname: residue segname (str)
            dummy_patch: patch to build new residue
            dummy_coords: coordinated of dummy atoms
        Returns:
            denovo_residue: complete residue (AtomGroup)
            bonds: list of all bonds
        """
        dummy_residue, dummy_atoms, bonds = self.init_new_residue(
            0, "DUMMY", "D", "DUM"
        )
        counter = 0
        for a in dummy_atoms:
            dummy_residue.select(f"name {a}").setCoords([dummy_coords[counter]])
            counter += 1
        denovo_residue, dele_atoms, bonds = self.build_from_patch(
            dummy_residue, resid, resname, chain, segname, dummy_patch
        )
        del dummy_residue
        # bonds.extend(bonds_p)
        return denovo_residue, dele_atoms, bonds

    def build_patch_missing_atom_coord(self, link_residue, residue, missing_atoms, ICs):
        """Builds all missing atoms in residue from a patch linking it to link_residue
        Parameters:
            link_residue: first residue in patch (AtomGroup)
            residue: second residue in patch (AtomGroup)
            missing_atoms: list of missing atom in second residue
            ICs: list of internal coordinate to build missing atoms
        """
        unsorted_graph, required_atoms = self.build_IC_graph(missing_atoms, ICs)
        sorted_graph = topological_sort(unsorted_graph)
        atoms = [g[0] for g in sorted_graph if g[0] in required_atoms]
        for a in atoms:
            atomname = a[1:]
            atom = residue.select(f"name {atomname}")
            ic = self.Topology.get_IC(ICs, a)
            if ic:
                ic = ic[0]
                xa_list = []
                for atom_ic in ic[0:4]:
                    if atom_ic.replace("*", "")[0] == "2":
                        sel = residue.select("name " + atom_ic.replace("*", "")[1:])
                        xa_list.append(sel.getCoords()[0])
                    else:
                        sel = link_residue.select(
                            "name " + atom_ic.replace("*", "")[1:]
                        )
                        xa_list.append(sel.getCoords()[0])
                atom.setCoords(self.build_cartesian(*xa_list, ic[8], ic[7], ic[6]))

    def build_missing_atom_coord(self, residue, missing_atoms, ICs):
        """Builds all missing atoms based on the provided internal coordinates
        Parameters:
            residue: Prody residue (AtomGroup)
            missing_atoms: list with all missing atom name
            ICs: list of internal coordinates for building missing atoms
        """
        unsorted_graph, required_atoms = self.build_IC_graph(missing_atoms, ICs)
        sorted_graph = topological_sort(unsorted_graph)
        atoms = [g[0] for g in sorted_graph if g[0] in required_atoms]
        for a in atoms:
            atom = residue.select(f"name {a}")
            ic = self.Topology.get_IC(ICs, a)
            if ic:
                ic = ic[0]
                xa1 = residue.select(f"name {ic[0]}").getCoords()[0]
                xa2 = residue.select(f"name {ic[1]}").getCoords()[0]
                xa3 = residue.select(f"name {ic[2].replace('*', '')}").getCoords()[0]
                atom.setCoords(self.build_cartesian(xa1, xa2, xa3, ic[8], ic[7], ic[6]))

    def build_cartesian(self, a1, a2, a3, r, theta, phi):
        """Builds missing atom from internal coordinates
        Parameters:
            a1: coordinates of atom1
            a2: coordinates of atom2
            a3: coordinates of atom3
            r: distance from atom3
            theta: angle between a2 a3 and missing atom
            phi: torsional angle formed by a1, a2, a3 and missing atom
        Returns:
            coordinates of missing atom
        """

        theta = np.radians(theta)
        phi = np.radians(phi)
        cost = np.cos(theta)
        sint = np.sin(theta)
        cosp = np.cos(phi)
        sinp = np.sin(phi)
        rjk = a2 - a3
        rjk /= np.linalg.norm(rjk)
        rij = a1 - a2
        cross = np.cross(rij, rjk)
        cross /= np.linalg.norm(cross)
        cross2 = np.cross(rjk, cross)
        cross2 /= np.linalg.norm(cross2)
        wt = [r * cost, r * sint * cosp, r * sint * sinp]
        newc = rjk * wt[0] + cross2 * wt[1] + cross * wt[2]
        return a3 + newc

    def get_bonds(self, residue):
        rn = residue.getRenames()[0]
        bonds = []
        for a1, a2 in itertools.pairwise(self.Topology[rn]["BOND"]):
            i1 = residue.select(f"name {a1}").getSerials[0]
            i2 = residue.select(f"name {a2}").getSerials[0]
            bonds += (i1, i2)
        return bonds
