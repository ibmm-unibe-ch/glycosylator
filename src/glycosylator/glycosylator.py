"""
----------------------------------------------------------------------------

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>

2016 Thomas Lemmin
----------------------------------------------------------------------------
"""


import copy
import os
import re
import sqlite3

import networkx as nx
import numpy as np
import prody

from .molecule import Molecule
from .molecule_builder import MoleculeBuilder
from .utils import *


class Glycosylator:
    def __init__(self, topofile, paramfile, force_field="charmm"):
        """
        Parameters:
            topofile: path to topology file
            paramfile: path to parameter file
            force_field: name of force field. Default charmm
        Initializes:
            builder: MoleculeBuilder
            connect_topology: dictionary describing the topology of known glycans
            glycan_keys: dictionary for identifying glycans (built from connect_topology)
            glycoprotein: Atomgroup representing the glycoprotein
            protein: Atomgroup without the identified glycans
            sequences: dictionary with protein chain as keys and sequences as value
            sequons: dictionary with sequon residue id as keys and sequon triplet as value
            glycanMolecules: dictionary with protein residue as keys and Molecules instances as values
            glycans: dictionary with protein residue as keys and graph of connected glycan as values
            names: dictionary of residue names for linked glycans
        """
        if force_field != "charmm":
            raise NotImplementedError(
                "Non-CHARMM force fields have not been implemented yet"
            )
        self.topofile = topofile
        self.builder = MoleculeBuilder(topofile, paramfile)
        self.connect_topology = {}
        self.glycan_keys = {}
        self.glycoprotein = None
        self.protein = None
        self.sequences = {}
        self.sequons = {}
        self.glycanMolecules = {}
        self.glycans = {}
        self.names = {}
        self.prefix = ["segment", "chain", "resid", "icode"]

    def _reset_glycoprotein(self):
        """Initializes all the variables"""
        self.glycoprotein = None
        self.protein = None
        self.sequences = {}
        self.sequons = {}
        self.glycanMolecules = {}
        self.glycans = {}
        self.names = {}

    def read_connectivity_topology(self, connect_file_path):
        """Parse file defining the topology of connectivity trees.
        This function will initialize connect_topology and glycan_keys

        Each connectivity is defined by
            RESI: name of the polyglycan
            UNIT: resname, [list of patches to residue]

        Parameter:
            fileName: path to connectivity file
        """
        with open(connect_file_path, "r") as f:
            lines = [line.strip() for line in f.readlines()]
        residue = {}
        self.connect_topology = {}
        nbr_units = 0
        for line in lines:  # Loop through each line
            line = line.split("!")[0].split()  # remove comments and endl
            if line:
                if line[0] == "RESI":
                    if residue:
                        residue["#UNIT"] = nbr_units
                        self.connect_topology[resname] = copy.copy(residue)
                    residue["UNIT"] = []
                    resname = line[1]
                    nbr_units = 0
                elif line[0] == "UNIT":
                    self.read_unit(line, residue)
                    nbr_units += 1
        residue["#UNIT"] = nbr_units
        self.connect_topology[resname] = copy.copy(residue)
        self.build_keys()

    def import_connectivity_topology(self, filename):
        """Import connectivity topology from sql database
        This function will initialize connect_topology
        Parameters:
            filename: path to database
        """
        try:
            conn = sqlite3.connect(filename)
        except:
            print(f"Error while connecting to the database {filename}")
            return -1
        cursor = conn.cursor()
        cursor.execute("SELECT * FROM glycans")
        glycans = cursor.fetchall()

        self.connect_topology = {}
        for glycan in glycans:
            name, tree = glycan
            residue = {}
            residue["UNIT"] = []
            nbr_unit = 0
            for unit in tree.split("|"):
                unit = unit.split(" ")
                nbr_unit += 1
                if len(unit) > 2:
                    residue["UNIT"].append([unit[0], unit[1], unit[2:]])
                else:
                    residue["UNIT"].append([unit[0], " ", []])

            residue["#UNIT"] = nbr_unit
            self.connect_topology[name] = residue
        self.build_keys()

    def add_glycan_to_connectivity_topology(self, name, connect_tree, overwrite=True):
        """Add new glycan to connect_topology dictionary
        Parameters:
            name: name of new glycan
            connect_tree: dictionary with connectivity tree
            overwrite: should an existing glycan be overwritten
        """
        if name in self.connect_topology and not overwrite:
            print(
                "Glycan with same name "
                + name
                + "already exists. Please change name or allow overwritting"
            )
            return -1
        self.connect_topology[name] = connect_tree

    def export_connectivity_topology(self, filename):
        """Export connectivity topology to sql database
        This function will
        """
        try:
            conn = sqlite3.connect(filename)
        except:
            print(f"Error while connecting to the database {filename}")
            return -1
        cursor = conn.cursor()
        tn = "glycans"
        gn = "glycan_name"
        gt = "glycan_tree"
        cursor.execute(f"DROP TABLE IF EXISTS {tn}")
        cursor.execute(f"CREATE TABLE {tn} ({gn} text, {gt} text)")

        for key in self.connect_topology.keys():
            units = self.connect_topology[key]["UNIT"]
            glycan = []
            for unit in units:
                v = []
                v.extend(unit[0:2])
                v.extend(unit[2])
                glycan.append(" ".join(v))
            glycan = "|".join(glycan)

            cursor.execute(f"INSERT INTO {tn} VALUES ('{key}', '{glycan}')")

        conn.commit()
        conn.close()

    def build_keys(self):
        self.glycan_keys = {}
        for res in self.connect_topology:
            key = [
                f"{r[0]} {' '.join(r[2])}" for r in self.connect_topology[res]["UNIT"]
            ]
            self.glycan_keys[res] = key

    def read_unit(self, unit, residue):
        if len(unit) > 2:
            residue["UNIT"].append([unit[1], unit[2], unit[3:]])
        else:
            residue["UNIT"].append([unit[1], " ", []])

    def build_glycan_gaph(self, connect_tree):
        """Builds a graph representation of a connectivity tree
        Parameters:
            connect_tree: dictionnary of connectivity of a poly glycan
        Returns:
            unsorted_graph: dictionary repersenting the connectivity graph
        """
        inv_connect_tree = {v: k for k, v in connect_tree.items()}
        unsorted_graph = {}
        for g in connect_tree:
            if not g in unsorted_graph:
                unsorted_graph[g] = []
            key = connect_tree[g]
            if not key:
                continue
            gr = inv_connect_tree[" ".join(key.split()[:-1])]
            if gr in unsorted_graph:
                unsorted_graph[gr].append(g)
            else:
                unsorted_graph[gr] = [g]
        return unsorted_graph

    def load_glycoprotein(self, protein):
        """Load and detects glycans in a glycoprotein
        Parameters:
            protein: Atomgroup of
        """
        self._reset_glycoprotein()
        if type(protein) == str:
            protein = prody.parsePDB(protein)
        self.glycoprotein = protein
        sel = self.glycoprotein.select("name CA")
        segn = sel.getSegnames()
        chids = sel.getChids()

        frag_ids = {f"{s},{c}" for s, c in zip(segn, chids)}
        for i in frag_ids:
            sel = ["name CA"]
            for p, s in zip(["segment ", "chain "], i.split(",")):
                if s:
                    sel.append(p + s)
            sel = " and ".join(sel)
            sel = self.glycoprotein.select(sel)
            self.sequences[i] = sel.getSequence()

            self.sequons.update(self.find_sequons(sel, self.sequences[i]))

        self.glycans, self.names = self.find_glycans("NGLB", "ASN")
        for g in self.glycans.values():
            nx.set_node_attributes(g[1], self.names, "resname")
        self.extract_glycans()

    def write_glycoprotein(self, filename):
        """Numpy array requires a continous memory block. When merging several atomGroup, the memory explodes. This function will alocate the memory first and then transfer propreties of AtomGroup"""
        natoms = self.protein.numAtoms()
        for g in self.glycanMolecules.values():
            natoms += g.atom_group.numAtoms()

        coords = np.zeros((natoms, 3))
        resn = []
        resi = []
        atoms_name = []
        chid = []
        segn = []
        occupancy = []
        beta = []
        serial = []
        element = []
        icode = []
        altloc = []

        # inilizalize variable
        ag = self.protein
        natoms = ag.numAtoms()
        coords[0:natoms, :] = ag.getCoords()
        atoms_name += list(ag.getNames())
        resn += list(ag.getResnames())
        resi += list(ag.getResnums())
        chid += list(ag.getChids())
        segn += list(ag.getSegnames())
        occupancy += list(ag.getOccupancies())
        beta += list(ag.getBetas())
        serials = range(1, natoms + 1)
        serial += serials
        icode += list(ag.getIcodes())
        element += list(ag.getElements())
        altloc += list(ag.getAltlocs())

        for g in self.glycanMolecules.values():
            ag = g.atom_group
            nnatoms = natoms + ag.numAtoms()
            coords[natoms:nnatoms, :] = ag.getCoords()
            resn += list(ag.getResnames())
            resi += list(ag.getResnums())
            chid += list(ag.getChids())
            segn += list(ag.getSegnames())
            occupancy += list(ag.getOccupancies())
            beta += list(ag.getBetas())
            atoms_name += list(ag.getNames())
            serials = range(natoms + 1, nnatoms + 1)
            serial += serials
            icode += list(ag.getIcodes())
            element += list(ag.getElements())
            altloc += list(ag.getAltlocs())
            natoms = nnatoms

        glycoprotein = prody.AtomGroup("Glycoprotein")
        glycoprotein.setCoords(coords)

        glycoprotein.setNames(atoms_name)
        glycoprotein.setResnums(resi)
        glycoprotein.setResnames(resn)
        glycoprotein.setChids(chid)
        glycoprotein.setSegnames(segn)
        glycoprotein.setOccupancies(occupancy)
        glycoprotein.setBetas(beta)
        glycoprotein.setSerials(serial)
        glycoprotein.setIcodes(icode)
        glycoprotein.setElements(element)
        glycoprotein.setAltlocs(altloc)
        prody.writePDB(filename, glycoprotein)

    def write_psfgen(self, dirName, proteinName=None):
        """This function will create the configure file for psfgen. The glycans will be split from the rest of the structure.
        It will be assumed that the topology file for the protein has already been previously built. Each glycan will be saved to a seperate file
        Parameters:
            dirName: name for output directory
            proteinName: path to topology (.psf) and coordinate (.pdb) files of the protein
        """
        try:
            os.mkdir(dirName)
        except:
            print("Directory ", dirName, " already exists. Please provide a new name")
            return -1
        os.chdir(dirName)
        psfbuffer = []
        psfbuffer.append("package require psfgen")
        psfbuffer.append("psfcontext reset")
        psfbuffer.append(f"topology {self.topofile}")

        for g in self.glycanMolecules.values():
            ag = g.atom_group
            segname = g.get_segname()
            prody.writePDB(f"{segname}.pdb", ag)
            psfbuffer.append(f"segment {segname} {{pdb {segname}.pdb}}")
            patches = g.get_patches()
            for patch in patches:
                id1, id2 = patch
                patch_name = patches[patch]
                s1, c, r1, i = id1.split(",")
                s2, c, r2, i = id2.split(",")
                psfbuffer.append(f"patch {patch_name} {s1}:{r1} {s2}:{r2}")
            psfbuffer.append(f"coordpdb {s1}.pdb")
            psfbuffer.append("guesscoord")
            psfbuffer.append("")
        psfbuffer.append("writepsf glycans.psf")
        psfbuffer.append("writepdb glycans.pdb")
        psfbuffer.append("")
        psfbuffer.append("psfcontext reset")
        psfbuffer.append(f"topology {self.topofile}")
        if proteinName:
            psfbuffer.append(f"readpsf {proteinName}.psf")
            psfbuffer.append(f"coordpdb {proteinName}.pdb")
        else:
            psfbuffer.append("#!!!!ADD PROTEIN NAME HERE!!!!")
            psfbuffer.append("readpsf protein.psf")
            psfbuffer.append("coordpdb protein.pdb")
        psfbuffer.append("readpsf glycans.psf")
        psfbuffer.append("coordpdb glycans.pdb")
        psfbuffer.append("")
        # TODO !!!hard coded!!! Should be patch dependent
        for res in self.glycanMolecules:
            s, c, r, i = res.split(",")
            sg, cg, rg, ig = self.glycanMolecules[res].rootRes.split(",")
            psfbuffer.append(f"patch NGLB {s}:{r} {sg}:{rg}")
        psfbuffer.append("regenerate angles dihedrals")
        if proteinName:
            psfbuffer.append(f"writepsf {proteinName}_glycosylated.psf")
            psfbuffer.append(f"writepdb {proteinName}_glycosylated.pdb")
        else:
            psfbuffer.append("writepsf glycoprotein.psf")
            psfbuffer.append("writepdb glycoprotein.pdb")
        psfbuffer.append("exit")

        with open("psfgen.conf", "w") as f:
            f.write("\n".join(psfbuffer))

    def get_residue(self, res_id):
        """Returns an AtomGroup of given atom id; composed of 'segname,chain,resid,icode,atomName'"""
        sel = []
        for p, s in zip(self.prefix, res_id.split(",")):
            if s:
                sel.append(f"{p} {s}")
        sel = " and ".join(sel)
        return self.glycoprotein.select(sel)

    def getSequence(self, sel):
        res = sel.getResnames()
        return "".join(aaa2a[r] for r in res)

    def get_start_resnum(self, chid):
        sel = ["name CA"]
        for p, s in zip(["segment ", "chain "], chid.split(",")):
            if s:
                sel.append(p + s)
        return self.glycoprotein.select(" and ".join(sel)).getResnums()[0]

    def get_glycans(self, sequons):
        """Returns a dictornary of glycans detected for sequons"""
        glycans = {}
        for s in sequons:
            if s in self.glycans:
                glycans[s] = self.glycans[s]
        return glycans

    def get_sequons(self, chid):
        """Returns all sequons of a chain (chid)"""
        return [k for k in self.sequons.keys() if chid in k[: len(chid)]]

    def find_sequons(self, sel, seq=None):
        """finds all sequon in AtomGroup"""
        if not seq:
            seq = sel.getSequence()
        sequons = {}
        res = sel.getResnums()
        icodes = sel.getIcodes()
        segn = sel.getSegnames()
        chid = sel.getChids()
        for m in re.finditer(r"(N[^P][S|T])", seq):
            idx = m.start()
            sequons[
                ",".join([segn[idx], chid[idx], str(res[idx]), icodes[idx]])
            ] = m.group()
        return sequons

    def find_glycans(self, patch, resname):
        """Looks from all molecules (not protein) which are linked (1.7A) to residue
        Parameters:
            patch: name of patch that connects molecule to residue
            resname: residue name to which the molecule is bound
        Returns:
            glycans: dictionary with residue id as key and a list of root unit and graph of glycan
            names: dictionary with glycan id as key and glycan resname as value
        """
        at = self.builder.Topology.patches[patch]["BOND"][0:2]
        for a in at:
            if a[0] == "1":
                a1 = a[1:]
            else:
                a2 = a[1:]
        # a2 is assumed to be the glycan atoms and a1 from protein
        sel = self.glycoprotein.select(
            f"(not protein and name {a2}) or (resname {resname} and name {a1})"
        )
        kd = prody.KDTree(sel.getCoords())
        kd.search(1.7)
        atoms = kd.getIndices()
        if atoms is None:
            return {}, {}
        G = nx.Graph()
        ids, rn = self.get_ids(sel)
        for a1, a2 in atoms:
            id1 = ids[a1]
            rn1 = rn[a1]
            id2 = ids[a2]
            rn2 = rn[a2]

            if id1 != id2:
                e = (id1, id2)
                G.add_edge(*e)

        names = self.connect_all_glycans(G)
        glycans = {}
        for graph in list(nx.connected_component_subgraphs(G)):
            for node in graph.nodes():
                if node not in names:
                    root_glycan = list(graph.neighbors(node))[0]
                    graph.remove_node(node)
                    glycans[node] = [root_glycan, graph]
                    break
        return glycans, names

    def extract_glycans(self):
        sel_all = []
        for k, v in self.glycans.items():
            r, g = v
            sel = []

            # !!!hard coded!!! Should be patch dependent
            selg = []
            for p, s in zip(self.prefix, r.split(",")):
                if s:
                    selg.append(f"{p} {s}")
            rootAtom = self.glycoprotein.select(
                f"{' and '.join(selg)} and name C1"
            ).getSerials()[0]

            for node in g.nodes():
                selg = []
                for p, s in zip(self.prefix, node.split(",")):
                    if s:
                        selg.append(f"{p} {s}")
                sel.append(" and ".join(selg))
            sel = f"({') or ('.join(sel)})"
            sel_all.append(sel)
            glycan = Molecule(self.glycoprotein.select(sel).copy(), rootAtom=rootAtom)
            glycan.interresidue_connectivity = g.copy()
            self.assign_patches(glycan)
            self.glycanMolecules[k] = glycan
        if sel_all:
            self.protein = self.glycoprotein.select(
                f"not (({') or ('.join(sel_all)}))"
            ).copy()
        else:
            self.protein = self.glycoprotein.copy()

    def connect_all_glycans(self, G):
        """Builds a connectivity graph for molecules (not protein) in AtomGroup. Edges with unknown patches will be removed
        Parameters:
            G: undirected graph of connected elements
        Returns:
            names: dictionary with residue id (get_id) as keys and residue name as value
        """
        sel = self.glycoprotein.select("not protein")
        kd = prody.KDTree(sel.getCoords())
        kd.search(1.7)
        atoms = kd.getIndices()
        atom_names = sel.getNames()
        names = {}
        ids, rn = self.get_ids(sel)

        for a1, a2 in atoms:
            id1 = ids[a1]
            rn1 = rn[a1]
            id2 = ids[a2]
            rn2 = rn[a2]
            if id1 in G and id1 not in names:
                names[id1] = rn1

            if id2 in G and id2 not in names:
                names[id2] = rn2

            if id1 != id2:
                e = (id1, id2)
                an1 = atom_names[a1]
                an2 = atom_names[a2]
                patch = self.find_patch(an1, an2)
                if patch:
                    G.add_edge(id1, id2, patch=patch, atoms=f"{an1}:{an2}")
                # G.add_edge(*e)
                names[id1] = rn1
                names[id2] = rn2

        return names

    def get_ids(self, sel):
        segn = sel.getSegnames()
        chid = sel.getChids()
        res = sel.getResnums()
        ic = sel.getIcodes()
        rn = sel.getResnames()
        ids = []
        for s, c, r, i in zip(segn, chid, res, ic):
            ids.append(",".join([s, c, str(r), i]))
        return ids, rn

    def get_id(self, sel):
        segn = sel.getSegnames()[0]
        chid = sel.getChids()[0]
        res = sel.getResnums()[0]
        ic = sel.getIcodes()[0]
        rn = sel.getResnames()[0]
        return ",".join([segn, chid, str(res), ic]), rn

    def glycosylate(
        self,
        glycan_name,
        glycan_molecule=None,
        link_residue=None,
        link_patch=None,
        template_glycan_tree={},
        template_glycan=None,
        chain="X",
        segname="G1",
    ):
        """Builds a polyglycan from a connectivity tree
        Parameters:
            glycan_name: name of poly glycan
            glycan_molecule: Molecule instance
            kwargs:

                segname: segname for new glyan (str)
                chain: chain for new
                link_residue: residue to link new glycan (AtomGroup)
                link_patch: name of patch (str) which should be used to link glycan
                template_glycan_tree: dictionary with the identified glycan that should be modified
                                key: identification of each unit in glycan
                                value: selection for residue; segname, chain, resid, icode
                template_glycan: AtomGroup containing the glycan that should be modified

        Returns:
            glycan: structure of the glycan (AtomGroup)
        """
        # if kwargs is not None:
        #    for key, value in kwargs.ashesteritems():

        resid = 1
        dummy_patch = "DUMMY_MAN"
        prefix = ["segment", "chain", "resid", "icode"]

        # either: retrieve glycan from data by name
        # or use the glycan provided
        if glycan_name in self.connect_topology:
            glycan_topo = self.get_connectivity_tree(glycan_name)
        elif glycan_molecule:
            chain = glycan_molecule.get_chain()
            segname = glycan_molecule.get_segname()
            self.assign_patches(glycan_molecule)
            glycan_topo, template_glycan_tree = self.build_glycan_topo(glycan_molecule)
            template_glycan = glycan_molecule.atom_group
        else:
            raise ValueError("Unknown glycan")

        # Variable for storing the new glycan
        glycan = None
        glycan_bonds = []
        built_glycan = {}
        inv_template_glycan_tree = {}
        dele_atoms = []
        if template_glycan_tree and template_glycan:
            inv_template_glycan_tree = {v: k for k, v in template_glycan_tree.items()}
            # resid = template_glycan.getResnums()[-1] + 1
            # chain = template_glycan.getChids()[0]
            # segname = template_glycan.getSegnames()[0]

        resids = []
        resid = 1
        sorted_units = sorted(glycan_topo.keys(), key=len)
        for unit in sorted_units:

            new_residue = None
            if unit:
                lunit = unit.split(" ")
                previous = " ".join(lunit[:-1])
            else:
                lunit = []
                previous = ""

            del_atom = []
            # check if residue exists
            if unit in inv_template_glycan_tree:
                sel = []

                for p, s in zip(prefix, inv_template_glycan_tree[unit].split(",")):
                    if s:
                        sel.append(f"{p} {s}")
                sel = " and ".join(sel)
                sel_residue = template_glycan.select(sel)
                # autopsf prefers if resids are incremental
                # resid =  sel_residue.getResnums()[0]

                if resid in resids:
                    resid = sorted(resids)[-1] + 1

                resids.append(resid)
                if sel_residue.getResnames()[0] == glycan_topo[unit]:
                    built_glycan[unit] = ",".join(
                        [
                            segname,
                            chain,
                            str(resid),
                        ]
                    )
                    # print resid,inv_template_glycan_tree[unit],','.join([segname, chain, str(resid),])
                    # built_glycan[unit] = inv_template_glycan_tree[unit]
                    new_residue, missing_atoms, bonds = self.builder.add_missing_atoms(
                        sel_residue, resid, chain, segname
                    )
                    ics = self.builder.Topology.topology[glycan_topo[unit]]["IC"]
                    self.builder.build_missing_atom_coord(
                        new_residue, missing_atoms, ics
                    )

                    # Check if first residue is linked to other residue
                    if link_residue and link_patch and not lunit:
                        del_atom, b = self.builder.apply_patch(
                            link_patch, link_residue, new_residue
                        )
                        bonds.extend(b)
                else:
                    print(
                        "The resnames do not match in connect tree. Residue will be build de novo"
                    )
            else:
                if resid in resids:
                    resid = sorted(resids)[-1] + 1
                resids.append(resid)

            # build first unit from DUMMY or from linked residue
            # if resid not in resids:
            #    resids.append(resid)

            if not lunit and not new_residue:
                if link_residue and link_patch:
                    new_residue, del_atom, bonds = self.builder.build_from_patch(
                        link_residue,
                        resid,
                        glycan_topo[unit],
                        chain,
                        segname,
                        link_patch,
                    )
                else:
                    new_residue, del_atom, bonds = self.builder.build_from_DUMMY(
                        resid, glycan_topo[unit], chain, segname, dummy_patch
                    )
            elif previous in built_glycan and lunit:
                patch = lunit[-1]
                sel = []
                for p, s in zip(prefix, built_glycan[previous].split(",")):
                    if s:
                        sel.append(f"{p} {s}")
                sel = " and ".join(sel)
                previous_residue = glycan.select(sel)
                if new_residue:
                    del_atom, b = self.builder.apply_patch(
                        patch, previous_residue, new_residue
                    )
                    bonds.extend(b)
                else:
                    new_residue, del_atom, bonds = self.builder.build_from_patch(
                        previous_residue,
                        resid,
                        glycan_topo[unit],
                        chain,
                        segname,
                        patch,
                    )
            elif lunit:
                print("Error in connect tree!! Glycans will not be built")
                return [], []

            built_glycan[unit] = ",".join(
                [
                    segname,
                    chain,
                    str(resid),
                ]
            )
            dele_atoms += del_atom
            glycan_bonds.extend(bonds)
            if glycan:
                glycan += new_residue
            else:
                glycan = new_residue
            resid += 1

        if dele_atoms:
            glycan = self.builder.delete_atoms(glycan, dele_atoms)
            # remove all non existing bonds
            tmp_bonds = []
            #            print 'Bonds',glycan_bonds
            #            print 'delele',dele_atoms
            for a1, a2 in glycan_bonds:
                if a1 in dele_atoms or a2 in dele_atoms:
                    continue
                else:
                    tmp_bonds.append((a1, a2))
            glycan_bonds = tmp_bonds

        # set serial number
        for i in range(len(glycan)):
            glycan.select(f"index {i}").setSerials([i + 1])

        return glycan, glycan_bonds

    def connect_tree_to_topology(self, connect_tree):
        """Converts a connect tree to a connect topology
        In connect tree the connectivity is represented as a string whereas it is a list in connect topology
        """
        connect_topology = {}
        units = connect_tree["UNIT"]
        unit_list = []
        n_unit = 0
        for unit in units:
            unit = filter(None, unit.split(" "))
            n_unit += 1
            if len(unit) > 1:
                unit_list.append([unit[0], "C1", unit[1:]])
            else:
                unit_list.append([unit[0], "", []])
        connect_topology["UNIT"] = unit_list
        connect_topology["#UNIT"] = n_unit
        return connect_topology

    def build_glycan_topology(self, patch="NGLB", glycanMolecules=None, build_all=True):
        """Builds the topology of all identified glycans
        Parameters:
            patch: patch used to link glycan to protien
            glycanMolecules: dictionary with id as key and Molecule instance as value. Default: glycans detected in by Glycosylator
            build_all: boolean defining if the topology has to be rebuilt
        """
        glycanMolecules = glycanMolecules if glycanMolecules else self.glycanMolecules

        for i, glycan in glycanMolecules.items():
            residue = self.get_residue(i)
            if not glycan.atom_group:
                continue
            if build_all:
                atom_group, bonds = self.glycosylate(
                    None,
                    glycan_molecule=glycan,
                    link_residue=residue,
                    link_patch=patch,
                )
                glycan.set_AtomGroup(
                    atom_group, rootAtom=atom_group.getSerials()[0], bonds=bonds
                )
            self.assign_patches(glycan)
            atom_type = self.assign_atom_type(glycan)
            glycan.set_atom_type(atom_type)
            glycan.define_torsionals(hydrogens=False)

    def get_connectivity_tree(self, glycan_name):
        """Builds connectivity tree of a glycan described in the connectivity topology
        Parameters:
            glycan_name: name of the glycan
        Returns:
            glycan_topo: dictionary with connectivy and resname of glycan
        """
        glycan_topo = {}
        units = self.connect_topology[glycan_name]["UNIT"]
        for unit in units:
            rn, a, p = unit
            # need to add blank at the begining of each key except the root key, which is blank ('')
            if p:
                # key = ' ' + ' '.join(p)
                key = " ".join(p)
            else:
                key = ""
            glycan_topo[key] = rn
        return glycan_topo

    def assign_patches(self, molecule):
        """Assignes patch names to each edge of interresidue digraph
        Parameters:
            molecule: Molecule object
        """
        G = molecule.interresidue_connectivity
        patches = {}
        atoms = nx.get_edge_attributes(G, "atoms")
        unknown_edges = []
        for e in G.edges():
            a1, a2 = atoms[e].split(":")
            patch = self.find_patch(a1, a2)
            if not patch:
                # print 'Unknown inter residue bond', a1, a2
                unknown_edges.append(e)
                continue
            u, v = e
            G[u][v]["patch"] = patch
        for e in unknown_edges:
            G.remove_edge(*e)

    def build_connectivity_tree(self, root_id, G):
        """Defines the connectivity within a glycan polymer
        Parameters:
         root_id: id of first residue of the glycan polymer ('segn,chid,resid')
         G: directed interresidue graph
        Returns:
         connect_tree: dictionary of glycan connectivity
        """
        paths = nx.shortest_path(G, source=root_id)
        # put root residue in dict
        connect_tree = {}
        connect_tree[root_id] = " "
        for n in paths:
            p = paths[n]
            edges = zip(p[1:], p[:-1])
            value = []
            for e2, e1 in edges:
                value += [G[e1][e2]["patch"]]
            connect_tree[n] = " ".join(value)
        connect_tree.items().sort(key=lambda id: len(id[1]))
        return connect_tree

    def build_connect_topology(self, molecule):
        """Builds connectivity topology
        Parameters:
            molecule: Molecule object
        """
        G = molecule.interresidue_connectivity
        connect_tree = self.build_connectivity_tree(molecule.rootRes, G)

        target = []
        for r in connect_tree.keys():
            if connect_tree[r]:
                target.append([G.node[r]["resname"], " ", connect_tree[r].split(" ")])
            else:
                target.append([G.node[r]["resname"], " ", []])

        len_target = len(target)
        residue = {}
        residue["UNIT"] = target
        residue["#UNIT"] = len_target
        return residue

    def build_glycan_topo(self, molecule):
        """Builds connectivity topology
        Parameters:
            molecule: Molecule object
        """
        G = molecule.interresidue_connectivity
        connect_tree = self.build_connectivity_tree(molecule.rootRes, G)

        glycan_topo = {}
        for r in connect_tree.keys():
            glycan_topo[connect_tree[r]] = G.node[r]["resname"]

        return glycan_topo, connect_tree

    def identify_glycan(self, molecule):
        """Identifies glycan name
        Parameters:
            molecule: Molecule object
        """
        G = molecule.interresidue_connectivity
        connect_tree = self.build_connectivity_tree(molecule.rootRes, G)

        target = []
        for r in connect_tree.keys():
            target.append(f"{G.node[r]['resname']} {connect_tree[r]}")
        len_target = len(target)

        for gk in self.glycan_keys:
            if (
                len(set(target) & set(self.glycan_keys[gk])) == len_target
                and len(set(self.glycan_keys[gk])) == len_target
            ):
                break
            else:
                gk = ""

        if gk:
            molecule.id = gk
        else:
            print("Unknown glycan")
            gk = ""
            molecule.id = ""
        return gk

    def write_connectivity_topology(self, glycan_name, connect_tree, fileName):
        """Write connectivity tree of a glycan
        Parameters:
            glycan_name: name of glycan (str)
            fileName: path to output file (str)
        """
        file = open(fileName, "w")  # open the file
        file.write(f"RESI {glycan_name}\n")
        units = connect_tree.items()
        units.sort(key=lambda id: len(id[1]))
        for unit in units:
            file.write(f"UNIT {unit[0].get_resname()} C1 {unit[1]}\n")
        file.close()

    def find_patch(self, atom1, atom2, anomer=""):
        """Finds patch that connects two atoms
        Currently based only on atom names
        Parameters:
            atom1: name of first atom (str)
            atom2: name of second atom (str)
        Returns:
            patch name
        """
        key = "-".join(sorted([atom1, atom2]))
        if key in self.builder.Topology.atomnames_to_patch:
            # for patch in self.builder.Topology.atomnames_to_patch[key]:
            #    self.builder.Topology[patch]['ICs']
            if anomer == "":
                return self.builder.Topology.atomnames_to_patch[key][0]
            else:
                for patch in self.builder.Topology.atomnames_to_patch[key]:
                    if patch[-2:] == anomer:
                        return patch
        else:
            return ""

    def assign_atom_type(self, molecule, connect_tree=None):
        """Returns a dictionary of atom types for a given molecule.
        Parameters:
           molecule: Molecule instance
           connect_tree: connectivity tree generated with "build_connectivity_tree". If not provided will be created in function
        Returns:
           atom_types: dictionary of atom types
                        key: atom serial number
                        value: atom name, atom type, charge, id, element
        """
        if not connect_tree:
            connect_tree = self.build_connectivity_tree(
                molecule.rootRes, molecule.interresidue_connectivity
            )
        inv_connect_tree = {v: k for k, v in connect_tree.items()}
        sorted_units = sorted(inv_connect_tree.keys(), key=len)
        atom_type = {}
        masses = self.builder.Topology.topology["MASS"]
        for unit in sorted_units:
            current = inv_connect_tree[unit]
            cur_res = molecule.get_residue(current)
            cur_rn = cur_res.getResnames()[0]

            cur_atoms = [a.strip() for a in cur_res.getNames()]
            cur_serial = cur_res.getSerials()
            cur_atom_serial = {}
            for a, s in zip(cur_atoms, cur_serial):
                cur_atom_serial[a] = s
            atoms = self.builder.Topology.get_atoms(cur_rn)
            for a in atoms:
                an = a[0].strip()
                if an in cur_atom_serial:
                    atom_type[cur_atom_serial[an]] = {
                        "name": an,
                        "type": a[1],
                        "charge": a[2],
                        "id": f"{current},{an}",
                        "element": masses[a[1]][-1],
                    }
            # correct atom type
            lunit = unit.split(" ")
            patch = lunit[-1]
            if not patch:
                continue
            previous = inv_connect_tree[" ".join(lunit[:-1])]
            pre_res = molecule.get_residue(previous)

            pre_atoms = [a.strip(" ") for a in pre_res.getNames()]
            pre_serial = pre_res.getSerials()
            pre_atom_serial = {}

            for a, s in zip(pre_atoms, pre_serial):
                pre_atom_serial[a] = s

            atoms = self.builder.Topology.patches[patch]["ATOM"]
            for a in atoms:
                an = a[0].strip(" ")
                if an[0] == "1":
                    if an[1:] in pre_atom_serial:
                        atom_type[pre_atom_serial[an[1:]]] = {
                            "name": an[1:],
                            "type": a[1],
                            "charge": a[2],
                            "id": f"{previous},{an[1:]}",
                            "element": masses[a[1]][-1],
                        }
                elif an[0] == "2":
                    if an[1:] in cur_atom_serial:
                        atom_type[cur_atom_serial[an[1:]]] = {
                            "name": an[1:],
                            "type": a[1],
                            "charge": a[2],
                            "id": f"{current},{an[1:]}",
                            "element": masses[a[1]][-1],
                        }
        return atom_type

    def export_patches(self, glycans):
        """Creates list of patches for parameter file for psfgen.
        Parameters:
            glycans: dictionary of glycans with sequon as key and glycan tree as values
        Returns:
            patches: list of patches
        """
        patches = []
        for k in glycans:
            root, tree = glycans[k]
            seg1, chid1, res1, i = k.split(",")
            seg2, chid2, res2, i = root.split(",")
            patches.append(f"patch NGLB {seg1}:{res1} {seg2}:{res2}")

            for e1, e2 in tree.edges():
                link = tree[e1][e2]["patch"]
                seg1, chid1, res1, i = e1.split(",")
                seg2, chid2, res2, i = e2.split(",")
                patches.append(f"patch {link} {seg1}:{res1} {seg2}:{res2}")
        return patches
