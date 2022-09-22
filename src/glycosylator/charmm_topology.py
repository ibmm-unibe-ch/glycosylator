import copy

from utils import *


class CHARMMTopology:
    """Class for parsing and storing CHARMM topology files.
    Attributes:
        topology: dictionary storing topology (RES)
                        key: resname
                        value: dictionary with ATOM, BOND, CHARGE and IC
                        key: MASS contains all masses of atoms in topology
         patches: dictionary storing patches (PRES)
                        key: patchname
                        value: dictionary with dele, ATOM, BOND, CHARGE and IC
         atomnames_to_patch: dictionary storing patch name to connect two atoms
                        key: atom1-atom2
                        value: patchname
    """

    def __init__(self, fileIn):
        self.topology = {}
        self.patches = {}
        self.atomnames_to_patch = {}
        self.read_topology(fileIn)

    def reset(self):
        """Resets previsously read topology"""
        self.topology = {}
        self.patches = {}
        self.atomnames_to_patch = {}

    def read_topology(self, fileIn):
        """Reads CHARMM topology file.
        Parameters:
            fileIn: path to topology file
        Initialize:
            topology: dictionary storing topology (RES)
                        key: resname
                        value: dictionary with ATOM, BOND, CHARGE and IC
                        key: MASS contains all masses of atoms in topology
            patches: dictionary storing patches (PRES)
                        key: patchname
                        value: dictionary with dele, ATOM, BOND, CHARGE and IC
            atomnames_to_patch: dictionary storing patch name to connect two atoms
                        key: atom1-atom2
                        value: patchname
        """
        lines = readLinesFromFile(fileIn)
        topo_type = ""
        residue = {}
        if "MASS" not in self.topology:
            self.topology["MASS"] = {}
        masses = self.topology["MASS"]

        for line in lines:  # Loop through each line
            line = line.split("\n")[0].split("!")[0].split()  # remove comments and endl
            if line:
                if line[0] == "RESI" or line[0] == "PRES":
                    if residue:
                        if topo_type == "RESI":
                            self.topology[resname] = copy.copy(residue)
                        elif topo_type == "PRES":
                            self.patches[resname] = copy.copy(residue)
                            key = "-".join(
                                sorted([residue["BOND"][0][1:], residue["BOND"][1][1:]])
                            )
                            # Allows multiple patches
                            if key in self.atomnames_to_patch:
                                res = self.atomnames_to_patch[key]
                                res.append(resname)
                                self.atomnames_to_patch[key] = res
                            else:
                                self.atomnames_to_patch[key] = [resname]
                    residue["dele"] = []
                    residue["ATOM"] = []
                    residue["BOND"] = []
                    residue["IC"] = []
                    topo_type = line[0]
                    resname = line[1]
                    residue["CHARGE"] = line[2]
                elif line[0] == "ATOM":
                    self.read_atom(line, residue)
                elif line[0] == "BOND":
                    self.read_bond(line, residue)
                elif line[0] == "IC":
                    self.read_ICs(line, residue)
                elif line[0] == "dele":
                    self.read_dele(line, residue)
                elif line[0] == "MASS":
                    self.read_mass(line, masses)

        if topo_type == "RESI":
            self.topology[resname] = copy.copy(residue)
        elif topo_type == "PRES":
            self.patches[resname] = copy.copy(residue)
            key = "-".join(sorted([residue["BOND"][0][1:], residue["BOND"][1][1:]]))
            # Allows multiple patches
            if key in self.atomnames_to_patch:
                res = self.atomnames_to_patch[key]
                res.append(resname)
                self.atomnames_to_patch[key] = res
            else:
                self.atomnames_to_patch[key] = [resname]
            # self.atomnames_to_patch[key] = resname

    def read_mass(self, mass, masses):
        mass[3] = float(mass[3])
        masses[mass[2]] = mass[3:]

    def read_atom(self, atom, residue):
        atom[3] = float(atom[3])
        residue["ATOM"].append(atom[1:])

    def read_dele(self, delatom, residue):
        residue["dele"] += delatom[2:]

    def read_bond(self, bond, residue):
        residue["BOND"] += bond[1:]

    def read_ICs(self, ic, residue):
        ic[5:] = map(float, ic[5:])
        residue["IC"].append(ic[1:])

    def get_IC(self, ics, atom):
        atom_ic = [ic for ic in ics if ic[3] == atom]
        return atom_ic

    def get_atom_name(self, ATOM):
        names = []
        for a in ATOM:
            names.append(a[0])
        return names

    def get_atoms(self, resname):
        return self.topology[resname]["ATOM"]

    def get_ICs(self, resname):
        return self.topology[resname]["IC"]

    def get_bonds(self, resname):
        return self.topology[resname]["BOND"]
