import copy

from utils import *


class CHARMMParameters:
    """Class for parsing and storing CHARMM parameters files.
    Attributes:
        parameters: dictionary storing parameters
            keys: 'BONDS', 'ANGLES', 'DIHEDRALS', 'NONBONDED', 'IMPROPER', 'NBFIX', 'CMAP' and 'ATOMS'
            values: dictionary of parameters
                    BONDS: atom1-atom2                  ->  k0, d0
                    ANGLES: atom1-atom2-atom3           ->  k0, a0, kub, d0
                    DIHEDRALS: atom1-atom2-atom3-atom4  ->  k0, n, dela
                    NONBONDED: atom1                    ->
                    IMPROPER: atom1-atom2-atom3-atom4   ->
                    NBFIX: atom1-atom2                  ->
                    CMAP:
                    ATOM: atom1                         -> mass
    """

    def __init__(self, fileIn):
        self.parameters = {}
        self.read_parameters(fileIn)

    def read_parameters(self, fileIn):
        """Reads CHARMM parameter file.
        Parameters:
            fileIn: path to parameter file
        Initializes:
            parameters: dictionary storing parameters
                keys: 'BONDS', 'ANGLES', 'DIHEDRALS', 'NONBONDED', 'IMPROPER', 'NBFIX', 'CMAP' and 'ATOMS'
                values: dictionary of parameters
                        BONDS: atom1-atom2                  ->  k0, d0
                        ANGLES: atom1-atom2-atom3           ->  k0, a0, kub, d0
                        DIHEDRALS: atom1-atom2-atom3-atom4  ->  k0, n, dela
                        NONBONDED: atom1                    ->  e, r/2, X, e(1-4) , r/2(1-4)
                        IMPROPER: atom1-atom2-atom3-atom4   ->
                        NBFIX: atom1-atom2                  ->
                        CMAP:
                        ATOM: atom1                         -> mass
        """
        lines = readLinesFromFile(fileIn)
        prm = {}
        prm_type = ""
        tags = [
            "BONDS",
            "ANGLES",
            "DIHEDRALS",
            "NONBONDED",
            "IMPROPER",
            "NBFIX",
            "CMAP",
            "ATOMS",
        ]
        # initialize parameter dictionary
        for t in tags:
            if not t in self.parameters:
                self.parameters[t] = {}
        for line in lines:  # Loop through each line
            line = line.split("\n")[0].split("!")[0].split()  # remove comments and endl
            if line:
                if line[0] in tags:
                    if prm:
                        if prm_type in self.parameters:
                            self.parameters[prm_type] = dict(
                                self.parameters[prm_type].items() + prm.items()
                            )
                        else:
                            self.parameters[prm_type] = copy.copy(prm)
                    prm_type = line[0]
                    prm = {}
                    continue
                if prm_type:
                    eval("self.read_" + prm_type + "(line, prm)")
        self.parameters[prm_type] = copy.copy(prm)

    def read_BONDS(self, bond, prm):
        #    CC311D     NC2D1     320.00    1.430
        if len(bond) == 4:
            prm["-".join(bond[0:2])] = map(float, bond[2:])
        else:
            print("Invalid BOND: " + " ".join(bond))

    def read_ANGLES(self, angle, prm):
        # CT1            CC321        HCA2     33.430        110.10     !22.53     2.17900
        if len(angle) == 5 or len(angle) == 7:
            prm["-".join(angle[0:3])] = map(float, angle[3:])
        else:
            print("Invalid ANGLE: " + " ".join(angle))

    def read_DIHEDRALS(self, dihe, prm):
        # CC321C    OC3C61    CC311D    NC2D1        0.62    1        0.0
        if len(dihe) == 7:
            key = "-".join(dihe[0:4])
            if key in prm:
                prm[key].append(map(float, dihe[4:]))
            else:
                prm[key] = [map(float, dihe[4:])]
        else:
            print("Invalid DIHEDRAL: " + " ".join(dihe))

    def read_IMPROPER(self, impr, prm):
        # NC2D1     CC2O1     CC311D    HCP1        20.00    0     0.00
        if len(impr) == 7:
            prm["-".join(impr[0:4])] = map(float, impr[4:])
        else:
            print("Invalid IMPROPER: " + " ".join(impr))

    def read_NONBONDED(self, vdw, prm):
        # CT3            0.0             -0.0780        2.040 ! 0.0 -0.01 1.9 ! alkane, 4/98, yin, adm jr.
        if len(vdw) == 4 or len(vdw) == 7:
            prm[vdw[0]] = map(float, vdw[2:])
        else:
            print("Invalid NONBONDED: " + " ".join(vdw))

    def read_NBFIX(self, nb, prm):
        # SOD        OCL            -0.07502        3.23
        if len(nb) == 4:
            prm["-".join(nb[0:2])] = map(float, nb[2:])
        else:
            print("Invalid NBFIX: " + " ".join(nb))

    def read_CMAP(self, cmap, prm):
        return -1

    def read_ATOMS(self, atom, prm):
        # MASS        31 H            1.00800
        if len(atom) == 4:
            prm[atom[2]] = float(atom[3])
        else:
            print("Invalid ATOM/MASS: " + " ".join(atom))
