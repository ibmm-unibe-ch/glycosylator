import networkx as nx
import prody

from .utils import *


class Molecule:
    """Class for saving a molecule
    Attributes:
            name: structure name
            chain: chain id
            segn: segname
            id: glycosylator id
            key: string representation of connectivity
            atom_group: Prody AtomGroup
            bonds: list of all bonds
            angles: list of all angles
            dihedrals: list of all dihedrals
            connectivity: graph for connectivity of molecule (bonds)
            directed_connectivity: directed acyclique graph of molecule
            interresidue_connectivity: directed acyclique graph representing the interresidue bonds
    """

    def __init__(self, name, chain="X", segn="X"):
        """initialize AtomGroup used to build pdb from scratch
        Parameters:
            name: structure name (str)
            chain: chain id (str)
            segn: segname (str)
        Initializes:
            atom_group: AtomGroup
            rootAtom: serial number of atom used as root for graphs
            bonds: list of all bonds
            angles: list of all angles
            dihedrals: list of all dihedrals
            connectivity: graph for connectivity of molecule (bonds)
            directed_connectivity: directed acyclique graph of molecule
            cycle_id: dictionary where keys are the serial number of atom in cycles and values the corresponding cycle in directed_connectivity
            torsionals: dihedral that can rotate (i.e. not in cycles)
            bond_length: dictionary of bond distance used to guess bonds. Keys are sorted by alphabetical order
        """
        self.name = name
        self.atom_group = prody.AtomGroup(self.name)
        self.chain = chain
        self.segn = segn
        self.rootAtom = ",".join(
            [
                segn,
                chain,
                "O",
            ]
        )
        self.bonds = []
        self.angles = []
        self.dihedrals = []
        self.connectivity = nx.Graph()
        self.directed_connectivity = nx.DiGraph()
        self.interresidue_connectivity = nx.DiGraph()
        self.cycle_id = {}
        self.torsionals = []
        self.bonded_uptodate = False

        self.prefix = ["segment", "chain", "resid", "icode"]
        # Defines distance for bond length between different element used in guess_bonds()
        self.elements = ["C", "H", "N", "O"]
        self.bond_length = {
            "C-H": 1.20,
            "H-N": 1.20,
            "H-O": 1.20,
            "C-C": 1.7,
            "C-N": 1.7,
            "C-O": 1.7,
        }

    def writePDB(self, filename, selection="all"):
        """Saves molecule to a PDB file
        Parameters:
            filename: path to PDB file
            selection: selection of a subset of the molecule (str)
        """
        prody.writePDB(filename, self.atom_group.select(selection))

    def read_molecule_from_PDB(self, filename, rootAtom=1, update_bonds=True, **kwargs):
        """Initialize molecule from a PDB file
        Parameters:
            filename: path to PDB file
            rootAtom: serial number of root atom
            update_bonds: guess bonds, angles, dihedrals and connectivity based on the distance between atoms in PDB
            **kwargs: any of the following which allows the selection of a subset of the PDB file
                subset: selection of a subset of the PDB file
                model: model number (int)
                chain: chain id (str)
        Initializes:
            atom_group
            chain: chain id
            segn: segment name
            connectivity: bonds, angles, dihedrals and graph
        """
        PDBmolecule = prody.parsePDB(filename, **kwargs)
        chain = set(PDBmolecule.getChids())
        segn = set(PDBmolecule.getSegnames())
        self.rootAtom = rootAtom
        if len(chain) == 1 and len(segn) == 1:
            self.atom_group = PDBmolecule
            self.chain = chain.pop()
            self.segn = segn.pop()
            a1 = self.atom_group.select(f"serial {str(self.rootAtom)}")
            segn = a1.getSegnames()
            chid = a1.getChids()
            res = a1.getResnums()
            ics = a1.getIcodes()
            for (
                s,
                c,
                r,
                ic,
            ) in zip(segn, chid, res, ics):
                self.rootRes = ",".join([s, c, str(r), ic])
            if update_bonds:
                self.update_connectivity()
            else:
                self.build_connectivity_graph()
        else:
            print("several chains are present in PDB. Please select only one molecule")
            return -1
        return 0

    def set_id(self):
        segn = self.atom_group.getSegnames()
        chid = self.atom_group.getChids()
        res = self.atom_group.getResnums()
        ics = self.atom_group.getIcodes()
        at = self.atom_group.getNames()
        ser = self.atom_group.getSerials()
        ids = {}

        for i, s, c, r, ic, a in zip(ser, segn, chid, res, ics, at):
            ids[i] = {"id": ",".join([s, c, str(r), ic, a])}

        self.connectivity.add_nodes_from(ids.items())

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

    def get_chain(self):
        return self.chain

    def get_segname(self):
        return self.segn

    def get_residue(self, res_id):
        """Returns an AtomGroup of given atom id; composed of 'segname,chain,resid,icode,atomName'"""

        sel = []
        for p, s in zip(self.prefix, res_id.split(",")):
            if s:
                sel.append(f"{p} {s}")
        sel = " and ".join(sel)
        return self.atom_group.select(sel)

    def get_atom(self, a_id, atom_name):
        """Returns an AtomGroup of given atom id; composed of 'segname,chain,resid,icode'"""
        for p, s in zip(self.prefix, a_id.split(",")):
            if s:
                sel.append(f"{p} {s}")
        sel = " and ".join(sel)
        sel += f" and name {atom_name}"
        return self.atom_group.select(sel)

    def set_atom_type(self, atom_type):
        """Assignes atom name, type and charge to each atom in the connectivity graph"""
        self.connectivity.add_nodes_from(atom_type.items())

    def build_connectivity_graph(self):
        """Builds a connectivity graph for molecules (not protein) in AtomGroup
        Parameters:
            G: undirected graph of connected elements
        Returns:
            names: dictionary with residue id (get_id) as keys and residue name as value
        """
        if not self.bonds:
            kd = prody.KDTree(self.atom_group.getCoords())
            kd.search(1.7)
            atoms = kd.getIndices()
        else:
            atoms = np.array(self.bonds) - 1

        atom_names = self.atom_group.getNames()
        ids, rn = self.get_ids(self.atom_group)
        if len(set(ids)) == 1:
            self.interresidue_connectivity.add_node(ids[0], resname=rn[0])
            return 0

        G = nx.Graph()
        for a1, a2 in atoms:
            id1 = ids[a1]
            rn1 = rn[a1]
            id2 = ids[a2]
            rn2 = rn[a2]

            if id1 != id2:
                e = (id1, id2)
                an1 = atom_names[a1]
                an2 = atom_names[a2]
                G.add_node(id1, resname=rn1)
                G.add_node(id2, resname=rn2)
                G.add_edge(id1, id2, patch="", atoms=f"{an1}:{an2}")
        # create directed graph and remove all unnecessary edges
        if G:
            self.interresidue_connectivity = G.to_directed()
            for edge in nx.dfs_edges(G, self.rootRes):
                edge = list(edge)
                edge.reverse()
                self.interresidue_connectivity.remove_edge(*edge)

    def get_names(self):
        """returns a dictironary with the residue ids as keys and residue names as values"""
        return nx.get_node_attributes(self.interresidue_connectivity, "resname")

    def get_patches(self):
        """returns a dictironary with the residue ids as keys and patches as values"""
        return nx.get_edge_attributes(self.interresidue_connectivity, "patch")

    def get_connectivity_atoms(self):
        """returns a dictironary with the residue ids as keys and patches as values"""
        return nx.get_edge_attributes(self.interresidue_connectivity, "atoms")

    def update_connectivity(self, update_bonds=True):
        """Updates all the connectivity (bond, angles, dihedrals and graphs)"""
        if update_bonds:
            self.guess_bonds()
        self.guess_angles()
        self.guess_dihedrals()
        self.update_graphs()
        self.bonded_uptodate = True

    def set_bonds(self, bonds, update_con=True):
        """Define list of bonds in a molecule
        Parameters:
            bonds: list of bonds
            update_con: update the connectivity with these new bonds
        """
        inv_atom = {
            v: k for k, v in nx.get_node_attributes(self.connectivity, "id").items()
        }
        newbonds = []
        for b1, b2 in bonds:
            if b1 in inv_atom and b2 in inv_atom:
                newbonds.append((inv_atom[b1], inv_atom[b2]))
            # else:
            # print 'Skipping bond', b1,b2
        self.connectivity = nx.Graph()
        self.connectivity.add_edges_from(newbonds)
        self.bonds = self.connectivity.edges()
        self.bonded_uptodate = False
        if update_con:
            self.update_connectivity(update_bonds=False)

    def set_torsional_angles(self, torsionals, angles, absolute=True):
        """Change the torsional angles for a list torsional list
        Parameters:
            torsionals_idx: list of torsional (index or name) angles which should be changed
            angles: list of angles in degrees
            absolute: define if the angles are
        """
        for torsional, theta in zip(torsionals, angles):
            self.rotate_bond(torsional, theta, absolute=absolute)

    def set_AtomGroup(self, AGmolecule, rootAtom=1, bonds=None, update_bonds=False):
        """Creates a Molecule instance from AtomGroup.
        Parameters:
            AGmolecule: prody AtomGroup object
            rootAtom: serial number of rootAtom. Default fist one
            bonds: list of bonds (e.g. generated with MoleculeBuilder)
            update_bonds: if bonds have to be guessed.
        """
        chain = set(AGmolecule.getChids())
        segn = set(AGmolecule.getSegnames())
        self.rootAtom = rootAtom
        if len(chain) == 1 and len(segn) == 1:
            self.atom_group = AGmolecule
            self.chain = chain.pop()
            self.segn = segn.pop()
            a1 = self.atom_group.select(f"serial {str(self.rootAtom)}")
            if not a1:
                self.rootAtom = self.atom_group.getSerials()[0]
                a1 = self.atom_group.select(f"serial {str(self.rootAtom)}")
            segn = a1.getSegnames()
            chid = a1.getChids()
            res = a1.getResnums()
            ics = a1.getIcodes()
            for (
                s,
                c,
                r,
                ic,
            ) in zip(segn, chid, res, ics):
                self.rootRes = ",".join([s, c, str(r), ic])

            if bonds:
                self.set_id()
                self.set_bonds(bonds)

            if update_bonds:
                self.update_connectivity()
            else:
                self.build_connectivity_graph()

        else:
            print(
                "Several chains are present in AtomGroup. Please select only one molecule"
            )
            return -1
        return 0

    def add_residue(self, residue, newbonds, dele_atoms=[]):
        """Add a new residue to a molecule
        Parameters:
            residue: proDy AtomGroup
            newbonds: list of bonds to be added
        """
        if self.atom_group.select(f"resid {ri}and chain {chid}"):
            print(
                "WARNING! A residue with the same id (resid and chain) already exists. The new residue has not been added"
            )
            return -1

        if dele_atoms:
            self.delete_atoms(dele_atoms)

        natoms = self.atom_group.numAtoms()
        self.atom_group += residue
        self.atom_group.setTitle(self.name)
        self.atom_group.setSerials(np.arange(natoms) + 1)

        self.connectivity.add_edges_from(np.array(newbonds) + natoms)
        # self.connectivity.remove_edges_from(delete_bonds)
        self.bonds = self.connectivity.edges()
        self.update_connectivity(update_bonds=False)

    def delete_atoms(self, dele_atoms):
        """removes atoms and bonds from molecule
        Parameter:
            del_atoms: serial number of atoms to be deleted
        """
        newbonds = []
        for a in sorted(dele_atoms, reverse=True):
            for b in self.bonds:
                if a in b:
                    continue
                elif a < b[0]:
                    b[0] -= 1
                elif a < b[1]:
                    b[1] -= 1
                newbonds.append(b)
        self.atom_group = self.atom_group.select(
            f"not serial {dele_atoms.join(' ')}"
        ).copy()
        # renumber atoms
        self.atom_group.setSerial(np.arange(self.atom_group.numAtoms()))
        self.atom_group.setTitle(self.name)
        self.bonds = newbonds
        self.update_connectivity(update_bonds=False)

    def guess_bonds(self, default_bond_length=1.6):
        """Searches for all bonds in molecule
        Parameters:
            default_bond_length: maximum distance between two connected heavy atoms (Angstrom), if not present in bond_length dictionary
        """
        self.connectivity = nx.Graph()
        for a in self.atom_group:
            bonds = []
            sel = ""
            a_elem = a.getElement()
            if a_elem:
                # use predefined bond length for atom pairs
                for e in self.elements:
                    key = "-".join(sorted([a_elem, e]))
                    if key in self.bond_length:
                        if sel:
                            sel += " or "
                        sel += (
                            "((element "
                            + e
                            + ") and (within "
                            + str(self.bond_length[key])
                            + " of serial "
                            + str(a.getSerial())
                            + "))"
                        )
            if not sel:
                sel = (
                    f"within {str(default_bond_length)} of serial {str(a.getSerial())}"
                )
            sel = f"({sel}) and (not serial {str(a.getSerial())})"

            # search for all neighboring atoms
            neighbors = self.atom_group.select(sel)
            if neighbors:
                for aa in neighbors:
                    bonds.append((a.getSerial(), aa.getSerial()))

            self.connectivity.add_edges_from(bonds)
        self.bonds = self.connectivity.edges()

    def guess_angles(self):
        """Searches for all angles in a molecule based on the connectivity"""
        self.angles = []
        for node in self.connectivity.nodes():
            self.angles.extend(self.find_paths(self.connectivity, node, 2))

    def guess_dihedrals(self):
        """Searches for all dihedrals in a molecule based on the connectivity"""
        self.dihedrals = []
        for node in self.connectivity.nodes():
            self.dihedrals.extend(self.find_paths(self.connectivity, node, 3))

    def find_paths(self, G, node, length, excludeSet=None):
        """Finds all paths of a given length
        Parameters:
            G: graph (netwrokx)
            node: starting node
            length: length of path
            excludedSet: set
        Returns:
            paths: list of all paths of a length starting from node
        """
        if excludeSet == None:
            excludeSet = {node}
        else:
            excludeSet.add(node)

        if length == 0:
            return [[node]]
        paths = [
            [node] + path
            for neighbor in G.neighbors(node)
            if neighbor not in excludeSet
            for path in self.find_paths(G, neighbor, length - 1, excludeSet)
        ]
        excludeSet.remove(node)
        return paths

    def set_rootAtom(self, rootAtom):
        """Sets the rootAtom and updates all the directed graph"""
        self.rootAtom = rootAtom
        a1 = self.atom_group.select(f"serial {str(self.rootAtom)}")
        self.rootRes = (
            a1.getSegnames()[0]
            + ","
            + a1.getChids()[0]
            + ","
            + str(a1.getResnums()[0] + "," + a1.getIcodes()[0])
        )
        self.update_graphs()

    def update_graphs(self):
        """Updates connectivity and directed graphs.
        - seaches for all cycles in connectivity graph
        - rebuilts acyclique directed connectivity graph
            starts from rootAtom
            cycle are collapsed
            nodes have the parameters
                    isclycle: True/False
                    cycle: serial number of all atoms in cycle
                    name of node:
                                    not cycle: serial number of atom
                                    cycle: string with all serial number joined by a '-'
        """
        cycles = nx.cycle_basis(self.connectivity, self.rootAtom)
        # flatten cycles
        self.cycle_id = {}
        for cycle in cycles:
            key = "-".join(map(str, cycle))
            for a in cycle:
                self.cycle_id[a] = key
        self.directed_connectivity = nx.DiGraph()

        for edge in nx.dfs_edges(self.connectivity, self.rootAtom):
            directed_edge = []
            atoms = []
            for node in edge:
                atoms.append(self.atom_group.select(f"serial {str(node)}"))
                if node in self.cycle_id:
                    key = self.cycle_id[node]
                    if key not in self.directed_connectivity:
                        self.directed_connectivity.add_node(
                            key,
                            iscycle=True,
                            cycle_id=map(int, self.cycle_id[node].split("-")),
                        )
                    directed_edge.append(key)
                else:
                    self.directed_connectivity.add_node(
                        node, iscycle=False, cycle_id=[]
                    )
                    directed_edge.append(node)
            if directed_edge[0] == directed_edge[1]:
                continue
            self.directed_connectivity.add_edge(directed_edge[0], directed_edge[1])
            a1, a2 = atoms
            if a1.getResnums()[0] != a2.getResnums()[0]:
                r1 = (
                    a1.getSegnames()[0]
                    + ","
                    + a1.getChids()[0]
                    + ","
                    + str(a1.getResnums()[0])
                    + ","
                    + a1.getIcodes()[0]
                )
                r2 = (
                    a2.getSegnames()[0]
                    + ","
                    + a2.getChids()[0]
                    + ","
                    + str(a2.getResnums()[0])
                    + ","
                    + a2.getIcodes()[0]
                )
                self.interresidue_connectivity.add_node(r1, resname=a1.getResnames()[0])
                self.interresidue_connectivity.add_node(r2, resname=a2.getResnames()[0])
                self.interresidue_connectivity.add_edge(
                    r1, r2, patch="", atoms=f"{a1.getNames()[0]}:{a2.getNames()[0]}"
                )

    def define_torsionals(self, hydrogens=True):
        """Builds a list with all the torsional angles that can rotate
        Parameters:
            hydrogens: include torsional involving terminal hydrogens
        Initializes:
            torsionals: a list of serial number of atom defining a torsional angle (quadruplet)
        """
        self.torsionals = []
        # cycles = nx.cycle_basis(self.connectivity, self.rootAtom)
        if not hydrogens:
            elements = nx.get_node_attributes(self.connectivity, "element")
            if not elements:
                print(
                    "Elements have not been defined (use assign_atom_type). Hydrogens cannot be excluded."
                )
                hydrogens = True

        for dihe in self.dihedrals:
            # check that quadruplet is not in cycle
            if dihe[1] in self.cycle_id and dihe[2] in self.cycle_id:
                continue

            d_dihe = []
            for a in dihe[1:-1]:
                if a in self.cycle_id:
                    a = self.cycle_id[a]
                d_dihe.append(a)

            # check direction of dihedral
            if self.directed_connectivity.has_edge(d_dihe[0], d_dihe[1]):
                pass
            elif self.directed_connectivity.has_edge(d_dihe[1], d_dihe[0]):
                dihe.reverse()
            else:
                continue

            # check if hydrogen
            if not hydrogens:
                if elements[dihe[0]] == "H" or elements[dihe[-1]] == "H":
                    continue
            # check if already in torsionals list
            exists = False
            if self.torsionals:
                for t in self.torsionals:
                    if dihe[1] == t[1] and dihe[2] == t[2]:
                        exists = True
                        break

            if exists:
                continue

            self.torsionals.append(dihe)

    def rotate_bond(self, torsional, theta, absolute=False):
        """Rotate the molecule around a torsional angle. Atom affected are in direct graph.
        Parameters:
            torsional: index of torsional angle (in torsionals) or list of serial number of atoms defining the torsional angle.
            theta: amount (degrees)
            absolute: defines if theta is the increment or the absolute value of the angle
        Returns:
            c_angle: angle before the rotation in degrees
        """
        if type(torsional) == int:
            torsional = self.torsionals[torsional]
        else:
            if torsional not in self.torsionals:
                print("Warning torsional")
                # return -1

        atoms = []
        a1 = torsional[-2]

        # check if last atom of torsional angle is in cycle
        if a1 in self.cycle_id:
            a1 = self.cycle_id[a1]
            atoms += a1.split("-")
        else:
            atoms.append(str(a1))

        for n in nx.descendants(self.directed_connectivity, a1):
            if type(n) == str:
                atoms += n.split("-")
            else:
                atoms.append(str(n))
        sel = self.atom_group.select(f"serial {' '.join(atoms)}")
        t = torsional[1:-1]
        v1, v2 = self.atom_group.select(f"serial {' '.join(map(str, t))}").getCoords()[
            np.argsort(t), :
        ]
        axis = v2 - v1
        c_angle = 0.0
        if absolute:
            c_angle = self.measure_dihedral_angle(torsional)
            theta = theta - c_angle

        coords = sel.getCoords()
        M = rotation_matrix(axis, np.radians(theta))
        coords = M.dot(coords.transpose())
        sel.setCoords(coords.transpose() + v2 - np.dot(M, v2))
        return c_angle

    def get_all_torsional_angles(self):
        """Computes all the torsional angles of the molecule
        Return:
            angles: angles in degrees
        """
        angles = []
        for torsional in self.torsionals:
            angle = self.measure_dihedral_angle(torsional)
            angles.append(angle)

        return angles

    def measure_dihedral_angle(self, torsional):
        """Calculates dihedral angle for 4 atoms.
        Parameters:
            torsional: list of atom serial numbers
        Returns:
            angle: dihedral angle in degrees
        """
        if type(torsional) == int:
            torsional = self.torsionals[torsional]
        #        else:
        #            if torsional not in self.torsionals:
        #                print "Warning Unknown torsional:", torsional
        #                return -1

        idx = np.argsort(torsional)
        vec_sel = self.atom_group.select(f"serial {' '.join(map(str, torsional))}")
        c0, c1, c2, c3 = vec_sel.getCoords()[idx, :]

        q1 = c1 - c0
        q2 = c2 - c1
        q3 = c3 - c2

        q1xq2 = np.cross(q1, q2)
        q2xq3 = np.cross(q2, q3)

        n1 = q1xq2 / np.sqrt(np.dot(q1xq2, q1xq2))
        n2 = q2xq3 / np.sqrt(np.dot(q2xq3, q2xq3))

        u1 = n2
        u3 = q2 / (np.sqrt(np.dot(q2, q2)))
        u2 = np.cross(u3, u1)

        cos_theta = np.dot(n1, u1)
        sin_theta = np.dot(n1, u2)
        angle = np.degrees(-np.arctan2(sin_theta, cos_theta))

        return angle

    def get_interresidue_torsionals(self, patches):
        connectivity_patches = self.get_patches()
        ids = nx.get_node_attributes(self.connectivity, "id")
        inv_ids = {v: k for k, v in ids.iteritems()}
        interresidue_torsionals = {}

        for r in connectivity_patches.keys():
            patch = connectivity_patches[r]
            # skip if unknown patch
            if patch not in patches:
                continue
            dihedrals = patches[patch]
            r1, r2 = r
            torsionals_idx = []
            delta_angles = []
            for dihe in dihedrals:
                atoms = dihe[0]
                serials = []
                # get id of central atoms
                for a in atoms:
                    if a[0] == "1":
                        serials.append(inv_ids[f"{r1},{a[1:]}"])
                    else:
                        serials.append(inv_ids[f"{r2},{a[1:]}"])
                # search for torsional angle index and measure delta angle
                for i, t in enumerate(self.torsionals):
                    if serials[1] == t[1] and serials[2] == t[2]:
                        delta_angles.append(
                            self.measure_dihedral_angle(t)
                            - self.measure_dihedral_angle(serials)
                        )
                        torsionals_idx.append(i)
                        break
            interresidue_torsionals["-".join(map(str, torsionals_idx))] = [
                patch,
                delta_angles,
            ]
        return interresidue_torsionals
