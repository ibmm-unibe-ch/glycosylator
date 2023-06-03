"""
This module defines File parsers and data classes to work with the CHARMM force field

Reading CHARMM files
====================

CHARMM files can be read using the `CHARMMTopology` and `CHARMMParameters` classes, depending on the file type.

.. code-block:: python

    from glycosylator.resources import CHARMMTopology, CHARMMParameters

    charmm_topology_file = "top_all36_prot.rtf"
    charmm_parameters_file = "par_all36_prot.prm"

    # load a CHARMM topology file
    top = CHARMMTopology.from_file(charmm_topology_file)

    # load a CHARMM parameters file
    par = CHARMMParameters.from_file(charmm_parameters_file)


Because parsing can be a time intensive task, it is recommended to `save` the parsed objects to a pickle file for later use.
In this case a pre-parsed topology or parameters object can be loaded using the `load` method.

.. code-block:: python

    # save the parsed objects to a pickle file
    top.save("top_all36_prot.pkl")
    par.save("par_all36_prot.pkl")

    # load a CHARMM topology file
    top = CHARMMTopology.load("top_all36_prot.pkl")

    # load a CHARMM parameters file
    par = CHARMMParameters.load("par_all36_prot.pkl")


Working with CHARMM objects
===========================

The `CHARMMTopology` and `CHARMMParameters` classes are include methods to work with the parsed data, 
such as `get_residue` or `get_mass`. They also support adding new data via the corresponding `add_{...}` methods.


Setting default CHARMM objects
==============================

The `glycosylator.utils.defaults` module pre-loads default CHARMM topology and parameters objects for convenience. Many methods that make use of 
these objects such as the `attach` methods of the `Scaffold` and `Molecule` classes also accept arguments for custom topology and parameters objects.
For convenience, a custom object can be set as the default, however, using the `set_default_topology` and `set_default_parameters` functions.

.. code-block:: python

    from glycosylator.utils import defaults
    
    # set a custom topology object as the default
    defaults.set_default_topology(top)


.. warning::

    The `set_default_{...}` functions include an argument `overwrite` that is set to `False` by default. If set to `True` the default object is permanently overwritten
    and will be used for all future sessions. This is not recommended as it may lead to unexpected behavior.

"""

import os
import pickle

import re
import warnings
import glycosylator.utils.abstract as _abstract
import glycosylator.utils.defaults as _defaults

# ===================================================================
# Base Parser
# ===================================================================


def load_topology(filename: str, set_default: bool = True) -> "CHARMMTopology":
    """
    Load a CHARMM topology from a pickle file.

    Parameters
    ----------
    filename: str
        The name of the topology file
    set_default: bool
        If `True`, the topology is set as the default topology

    Returns
    -------
    CHARMMTopology
        The topology object
    """
    top = CHARMMTopology.load(filename)
    if set_default:
        _defaults.set_default_topology(top)
    return top


def read_topology(filename: str, set_default: bool = True) -> "CHARMMTopology":
    """
    Read a CHARMM topology file.

    Parameters
    ----------
    filename: str
        The name of the topology file
    set_default: bool
        If `True`, the topology is set as the default topology

    Returns
    -------
    CHARMMTopology
        The parsed topology object
    """
    top = CHARMMTopology.from_file(filename)
    if set_default:
        _defaults.set_default_topology(top)
    return top


def save_topology(filename: str, topology: "CHARMMTopology" = None):
    """
    Save a CHARMM topology to a pickle file.

    Parameters
    ----------
    filename: str
        The name of the topology file
    topology: CHARMMTopology
        The topology object. If `None`, the default topology is used.
    """
    if topology is None:
        topology = _defaults.get_default_topology()
    topology.save(filename)


def load_parameters(filename: str, set_default: bool = True) -> "CHARMMParameters":
    """
    Load a CHARMM parameters from a pickle file.

    Parameters
    ----------
    filename: str
        The name of the parameters file
    set_default: bool
        If `True`, the parameters are set as the default parameters
    Returns
    -------
    CHARMMParameters
        The parameters object
    """
    par = CHARMMParameters.load(filename)
    if set_default:
        _defaults.set_default_parameters(par)
    return par


def read_parameters(filename: str, set_default: bool = True) -> "CHARMMParameters":
    """
    Read a CHARMM parameters file.

    Parameters
    ----------
    filename: str
        The name of the parameters file
    set_default: bool
        If `True`, the parameters are set as the default parameters.
    Returns
    -------
    CHARMMParameters
        The parsed parameters object
    """
    par = CHARMMParameters.from_file(filename)
    if set_default:
        _defaults.set_default_parameters(par)
    return par


def save_parameters(filename: str, parameters: "CHARMMParameters" = None):
    """
    Save a CHARMM parameters to a pickle file.

    Parameters
    ----------
    filename: str
        The name of the parameters file
    parameters: CHARMMParameters
        The parameters object. If `None`, the default parameters are used.
    """
    if parameters is None:
        parameters = _defaults.get_default_parameters()
    parameters.save(filename)


def has_patch(name: str) -> bool:
    """
    Check if a patch is defined in the CHARMM topology file.

    Parameters
    ----------
    name: str
        The name of the patch

    Returns
    -------
    bool
        `True` if the patch is defined, `False` otherwise
    """
    return name in _defaults.get_default_topology().patches


def get_patch(name: str):
    """
    Get a patch from the CHARMM topology file.

    Parameters
    ----------
    name: str
        The name of the patch

    Returns
    -------
    Patch
        The patch object
    """
    return _defaults.get_default_topology().get_patch(name)


def add_patch(patch, overwrite: bool = False):
    """
    Add a patch to the CHARMM topology file.

    Parameters
    ----------
    patch: Patch
        The patch object
    overwrite: bool
        If `True`, the topology with the added patch is saved to a pickle file and will be used as the default topology for all future sessions.
    """
    top = _defaults.get_default_topology()
    top.add_patch(patch)
    _defaults.set_default_topology(top, overwrite)


class CHARMMParser:
    """
    The base class to parse CHARMM files
    """

    def __init__(self) -> None:
        self.id = None
        self._dict = {
            "masses": {},
        }
        self._file = None
        self._was_loaded_from_pickle = False

    @classmethod
    def load(cls, filename: str):
        """
        Load a pre-parsed object from a pickle file.

        Parameters
        ----------
        filename: str
            The path to the pickle file
        """
        new = cls()
        _data = pickle.load(open(filename, "rb"))
        new._vet_load(_data)
        new._dict = _data
        new._file = filename
        new.id = os.path.basename(filename)
        new._was_loaded_from_pickle = True
        return new

    @classmethod
    def from_file(cls, filename: str):
        """
        Parse a CHARMM file and return a new object.

        Parameters
        ----------
        filename: str
            The path to the CHARMM file
        """
        new = cls()
        new._parse(filename)
        new._file = filename
        new.id = os.path.basename(filename)
        return new

    @property
    def masses(self):
        return self._dict["masses"]

    def save(self, filename: str = None):
        """
        Save the data dictionary as binary file that can be directly loaded again.

        Parameters
        ----------
        filename: str
            The path to the pickle file to save in. By default, this will be the same filename as the one from which the data was loaded or parsed (adding the file-suffix `.pkl`)
        """
        if not filename:
            if not self._file:
                raise ValueError(
                    "No filename was given and no filename from a source file is available!"
                )
            filename = self._file
            if not self._was_loaded_from_pickle and not filename.endswith(".pkl"):
                filename += ".pkl"
        pickle.dump(self._dict, open(filename, "wb"))

    def get_mass(self, id):
        """
        Get a mass by its atom type ID

        Parameters
        ----------
        id : str
            The ID of the atom type

        Returns
        -------
        float
            The mass
        """
        if isinstance(id, (list, tuple)):
            return [self.get_mass(i) for i in id]
        return self._dict["masses"].get(id)

    def add_mass(self, key, mass):
        """
        Add a mass to the topology

        Parameters
        ----------
        key : str
            The ID of the atom type
        mass : float
            The mass of the atom type
        """
        if isinstance(key, (list, tuple)):
            for k, m in zip(key, mass):
                self.add_mass(k, m)
            return
        self._dict["masses"][key] = float(mass)

    def _parse(self, filename: str):
        """
        Parse a CHARMM file and store the data in a dictionary

        Parameters
        ----------
        filename: str
            The path to the CHARMM file
        """
        raise NotImplementedError(
            "This method needs to be implemented by the subclass!"
        )

    def _vet_load(self, _data):
        """
        An optional method to vet the loaded data for specific properties
        """
        pass

    @staticmethod
    def _read_line(line: str):
        """
        Reads a line of the CHARMM file

        Parameters
        ----------
        line : str
            The line of the file

        Returns
        -------
        list
            The line split into its components
        """
        return line.strip().split("!")[0].split()

    def __repr__(self):
        return f"{self.__class__.__name__}({self.id})"


# ===================================================================
# Topology
# ===================================================================


class CHARMMTopology(CHARMMParser):
    """
    A data class and file parser for CHARMM topology files.

    Attributes
    ----------
    id : str
        The ID of the topology file
    residues : list of AbstractResidue
        The residues in the topology file
    patches : list of AbstractPatch
        The patches in the topology file
    masses : list of floats
        The masses in the topology file
    """

    __tags__ = (
        "MASS",
        "RESI",
        "PRES",
        "GROUP",
        "ATOM",
        "BOND",
        "DOUBLE",
        "IMPROPER",
        "IMPHI",
        "DIHE",
        "PATCH",
        "END",
    )

    def __init__(self, id=None):
        super().__init__()
        self._dict["residues"] = {}
        self._dict["patches"] = {}
        self.id = id

    @property
    def residues(self):
        return list(self._dict["residues"].values())

    @property
    def patches(self):
        return list(self._dict["patches"].values())

    def has_residue(self, id):
        """
        Check if a residue is in the topology

        Parameters
        ----------
        id : str
            The ID of the residue

        Returns
        -------
        bool
            True if the residue is in the topology
        """
        if isinstance(id, (list, tuple)):
            return all(self.has_residue(i) for i in id)
        return id in self._dict["residues"]

    def get_residue(self, id):
        """
        Get a residue by its ID

        Parameters
        ----------
        id : str
            The ID of the residue

        Returns
        -------
        AbstractResidue
            The residue
        """
        if isinstance(id, (list, tuple)):
            return [self.get_residue(i) for i in id]
        return self._dict["residues"].get(id)

    def find_residue(self, *atoms):
        """
        Find a residue by its atoms. This will require the
        residue to have all atoms in the list.

        Parameters
        ----------
        atoms : tuple of AbstractAtom
            The atoms to search for

        Returns
        -------
        AbstractResidue
            The residue
        """
        for residue in self.residues:
            if all(residue.has_atom(atom) for atom in atoms):
                return residue
        return None

    def add_residue(self, residue):
        """
        Add a residue to the topology

        Parameters
        ----------
        residue : AbstractResidue
            The residue to add
        """
        if isinstance(residue, (list, tuple)):
            for r in residue:
                self.add_residue(r)
            return
        elif not isinstance(residue, _abstract.AbstractResidue):
            raise TypeError("The residue must be an instance of AbstractResidue")
        self._dict["residues"][residue.id] = residue

    def get_patch(self, id):
        """
        Get a patch by its ID

        Parameters
        ----------
        id : str
            The ID of the patch

        Returns
        -------
        AbstractPatch
            The patch
        """
        if isinstance(id, (list, tuple)):
            return [self.get_patch(i) for i in id]
        return self._dict["patches"].get(id)

    def add_patch(self, patch):
        """
        Add a patch to the topology

        Parameters
        ----------
        patch : AbstractPatch
            The patch to add
        """
        if isinstance(patch, (list, tuple)):
            for p in patch:
                self.add_patch(p)
            return
        elif not isinstance(patch, _abstract.AbstractPatch):
            raise TypeError("The patch must be an instance of AbstractPatch")
        self._dict["patches"][patch.id] = patch

    def _parse(self, filename: str):
        """
        Reads and parses the data from a CHARMM Topology file

        Parameters
        ----------
        filename : str
            The path to the CHARMM Topology file
        """
        with open(filename, "r") as file:
            lines = file.read().split("\n")  # readlines but remove the endlines
            lines = [
                line.strip().split("!") for line in lines
            ]  # get rid of all comments
            lines = [
                line[0] for line in lines if line[0] != ""
            ]  # get rid of all empty lines

        idx = 0
        while idx < len(lines):
            line = lines[idx]
            if line.startswith("MASS"):
                self._parse_mass(line)
            elif line.startswith("RESI"):
                idx = self._parse_residue(lines, idx)
            elif line.startswith("PRES"):
                idx = self._parse_patch(lines, idx)
            idx += 1

        self._adopt_atom_masses()
        self._make_ICs_improper()
        self._file = filename

    def _parse_mass(self, line: str):
        """
        Parses the mass information from a line of the topology file

        Parameters
        ----------
        line : str
            The line of the topology file
        """
        line = self._read_line(line)
        id = line[2]
        mass = float(line[3])
        self._dict["masses"][id] = mass

    def _parse_residue(self, lines: list, idx: int):
        """
        Parses the residue information from a line of the topology file

        Parameters
        ----------
        lines: list
            The list of lines of the topology file
        idx: int
            The index of the current line
        """
        line = self._read_line(lines[idx])

        residue = _abstract.AbstractResidue(id=line[1])

        idx += 1
        while idx < len(lines):
            line = self._read_line(lines[idx])
            start = line[0]

            if start == "" or re.match("GROU(P| )", start):
                idx += 1
                continue

            if start == "ATOM":
                self._parse_atom(line, residue)

            elif start == "BOND":
                self._parse_bond(line, residue)

            elif start == "IC":
                self._parse_ic(line, residue)

            elif start in self.__tags__:
                idx -= 1
                break
            idx += 1

        self._dict["residues"][residue.id] = residue
        return idx

    def _parse_patch(self, lines: list, idx: int):
        """
        Parses the patch information from a line of the topology file

        Parameters
        ----------
        lines: list
            The list of lines of the topology file
        idx: int
            The index of the current line
        """
        line = self._read_line(lines[idx])

        patch = _abstract.AbstractPatch(id=line[1])

        idx += 1
        while idx < len(lines):
            line = self._read_line(lines[idx])
            start = line[0]

            if start == "" or re.match("GROU(P| )", start):
                idx += 1
                continue

            if start == "ATOM":
                self._parse_atom(line, patch)

            elif re.match("dele(te)?", start.lower()):
                self._parse_delete(line, patch)

            elif start == "BOND":
                self._parse_bond(line, patch)

            elif start == "IC":
                self._parse_ic(line, patch)

            elif start in self.__tags__:
                idx -= 1
                break
            idx += 1

        self._dict["patches"][patch.id] = patch
        return idx

    def _adopt_atom_masses(self):
        """
        Adopt the atom masses from the topology file
        """
        for residue in self.residues:
            for atom in residue.atoms:
                if atom.mass is None:
                    atom.mass = self._dict["masses"][atom.type]

    def _make_ICs_improper(self):
        """
        Ensure that improper internal coordinates are also labelled as such
        based on whether their third atom is connected to the first atom
        """
        for residue in self.residues:
            for ic in residue.internal_coordinates:
                if not ic.improper and residue.get_bond(ic.atom1, ic.atom3):
                    ic.improper = True
                    # Actually a bad idea to do this,
                    # because the file still stores 1-2 lengths instead of 1-3 lengths!
                    # if not ic.bond_length_13:
                    #     ic.bond_length_13 = ic.bond_length_12
                    #     ic.bond_length_12 = None

    def _vet_load(self, _data):
        """
        Checks that the data loaded from a file is valid
        """
        if not isinstance(_data, dict):
            raise TypeError("The file must contain a dictionary object")
        if (
            "residues" not in _data.keys()
            or "patches" not in _data.keys()
            or "masses" not in _data.keys()
        ):
            raise KeyError(
                "The dictionary must contain 'residues', 'masses', and 'patches' keys"
            )

    @staticmethod
    def _parse_ic(line: list, obj):
        """
        Parses the internal coordinate information from a line of the topology file

        Parameters
        ----------
        line : list
            The line of the topology file, split into a list
        obj : AbstractResidue or AbstractPatch
            The object to which the internal coordinate should be added
        """
        is_improper = "*" in line[3]
        if is_improper:
            line[3] = line[3][1:]

        atom1 = line[1]
        atom2 = line[2]
        atom3 = line[3]
        atom4 = line[4]

        if isinstance(obj, _abstract.AbstractResidue):
            atom1 = obj.get_atom(atom1)
            atom2 = obj.get_atom(atom2)
            atom3 = obj.get_atom(atom3)
            atom4 = obj.get_atom(atom4)

            if atom1 is None or atom2 is None or atom3 is None or atom4 is None:
                warnings.warn(
                    f"[ignoring line] Found an invalid internal coordinate in {line}"
                )
                return

        if is_improper:
            _bond_lengths = {
                "bond_length_12": None,
                "bond_length_13": float(line[5]),
                "bond_length_34": float(line[9]),
            }
        else:
            _bond_lengths = {
                "bond_length_12": float(line[5]),
                "bond_length_13": None,
                "bond_length_34": float(line[9]),
            }

        ic = _abstract.AbstractInternalCoordinates(
            atom1=atom1,
            atom2=atom2,
            atom3=atom3,
            atom4=atom4,
            bond_angle_123=float(line[6]),
            dihedral=float(line[7]),
            bond_angle_234=float(line[8]),
            improper=is_improper,
            **_bond_lengths,
        )

        obj.add_internal_coordinates(ic)

    @staticmethod
    def _parse_atom(line: list, obj):
        """
        Parses the atom information from a line of the topology file

        Parameters
        ----------
        line : list
            The line of the topology file, split into a list
        obj : AbstractResidue or AbstractPatch
            The object to which the atom belongs
        """
        atom = _abstract.AbstractAtom(id=line[1], type=line[2], charge=float(line[3]))
        obj.add_atom(atom)

    @staticmethod
    def _parse_bond(line: list, obj):
        """
        Parses the bond information from a line of the topology file

        Parameters
        ----------
        line : list
            The line of the topology file, split into a list
        obj : AbstractResidue or AbstractPatch
            The object to which the bond belongs
        """
        line = line[1:]

        # split line into tuple pairs
        line = [(line[i], line[i + 1]) for i in range(0, len(line), 2)]

        if isinstance(obj, _abstract.AbstractResidue):
            for a1, a2 in line:
                atom1 = obj.get_atom(a1)
                atom2 = obj.get_atom(a2)
                if atom1 is None or atom2 is None:
                    warnings.warn(f"[ignoring bond] Found an invalid bond in {line}")
                    return

                bond = _abstract.AbstractBond(atom1, atom2)
                obj.add_bond(bond)

        elif isinstance(obj, _abstract.AbstractPatch):
            for a1, a2 in line:
                bond = _abstract.AbstractBond(a1, a2)
                obj.add_bond(bond)

    @staticmethod
    def _parse_delete(line: list, patch):
        """
        Parses the delete information from a line of the topology file

        Parameters
        ----------
        line : list
            The line of the topology file, split into a list
        patch : AbstractPatch
            The patch to which the delete belongs
        """
        id = line[2]
        patch.add_delete(id)


# ===================================================================
# Parameters
# ===================================================================


class CHARMMParameters(CHARMMParser):
    """
    A data class and file parser for CHARMM Parameter files.
    """

    __tags__ = (
        "BONDS",
        "ANGLES",
        "DIHEDRALS",
        "NONBONDED",
        "IMPROPER",
        "NBFIX",
        "CMAP",
        "ATOMS",
    )
    """
    The valid tags within the Parameters file
    """

    def __init__(self, id=None):
        super().__init__()
        self._dict["bonds"] = {}
        self._dict["angles"] = {}
        self._dict["dihedrals"] = {}
        self._dict["impropers"] = {}
        self._dict["nonbonded"] = {}

        self.id = id

    @property
    def bonds(self):
        return list(self._dict["bonds"].values())

    @property
    def angles(self):
        return list(self._dict["angles"].values())

    @property
    def dihedrals(self):
        return list(self._dict["dihedrals"].values())

    @property
    def impropers(self):
        return list(self._dict["impropers"].values())

    @property
    def nonbonded(self):
        return list(self._dict["nonbonded"].values())

    def get_angle(self, atoms):
        """
        Gets an angle from the parameters

        Parameters
        ----------
        atoms : tuple
            The atoms in the angle

        Returns
        -------
        AbstractAngle
            The angle
        """
        return self._dict["angles"].get(tuple(atoms))

    def add_angle(self, angle):
        """
        Adds an angle to the parameters

        Parameters
        ----------
        angle : AbstractAngle
            The angle to add
        """
        self._dict["angles"][angle.atoms] = angle

    def get_improper(self, atoms):
        """
        Gets an improper from the parameters

        Parameters
        ----------
        atoms : tuple
            The atoms in the improper

        Returns
        -------
        AbstractImproper
            The improper
        """
        return self._dict["impropers"].get(tuple(atoms))

    def add_improper(self, improper):
        """
        Adds an improper to the parameters

        Parameters
        ----------
        improper : AbstractImproper
            The improper to add
        """
        self._dict["impropers"][improper.atoms] = improper

    def get_dihedral(self, atoms):
        """
        Gets a dihedral from the parameters

        Parameters
        ----------
        atoms : tuple
            The atoms in the dihedral

        Returns
        -------
        AbstractDihedral
            The dihedral
        """
        return self._dict["dihedrals"].get(tuple(atoms))

    def add_dihedral(self, dihedral):
        """
        Adds a dihedral to the parameters

        Parameters
        ----------
        dihedral : AbstractDihedral
            The dihedral to add
        """
        self._dict["dihedrals"][dihedral.atoms] = dihedral

    def get_bond(self, *atoms):
        """
        Get a bond by its atom ids

        Parameters
        ----------
        atoms : tuple
            The atom ids

        Returns
        -------
        AbstractBond
            The bond
        """
        if len(atoms) != 2:
            raise IndexError(
                "Bonds can only be retrieved through a two-long tuple of participating atoms."
            )
        return self._dict["bonds"].get(tuple(atoms))

    def add_bond(self, bond):
        """
        Adds a bond to the parameters

        Parameters
        ----------
        bond : AbstractBond
            The bond to add
        """
        self._dict["bonds"][bond.atoms] = bond

    def get_nonbonded(self, atom):
        """
        Gets a nonbonded from the parameters

        Parameters
        ----------
        atom : AbstractAtom or str
            The atom

        Returns
        -------
        AbstractNonbonded
            The nonbonded
        """
        if isinstance(atom, _abstract.AbstractAtom):
            atom = atom.type
        return self._dict["nonbonded"].get(atom)

    def add_nonbonded(self, nonbonded):
        """
        Adds a nonbonded to the parameters

        Parameters
        ----------
        nonbonded : AbstractNonbonded
            The nonbonded to add
        """
        self._dict["nonbonded"][nonbonded.atom] = nonbonded

    def _parse(self, filename):
        """
        Parses a CHARMM parameter file

        Parameters
        ----------
        filename : str
            The path to the parameter file
        """
        with open(filename, "r") as f:
            lines = f.readlines()

        idx = 0
        while idx < len(lines):
            line = self._read_line(lines[idx])
            if len(line) == 0:
                idx += 1
                continue

            start = line[0]

            if start == "BONDS":
                idx = self._parse_bonds(lines, idx)

            elif start == "ANGLES":
                idx = self._parse_angles(lines, idx)

            elif start == "DIHEDRALS":
                idx = self._parse_dihedrals(lines, idx)

            elif start == "MASS":
                self._parse_mass(line)

            elif start == "IMPROPER":
                idx = self._parse_impropers(lines, idx)

            elif start == "NONBONDED":
                idx = self._parse_nonbonded(lines, idx)

            idx += 1

    def _parse_nonbonded(self, lines: list, idx: int):
        """
        Parses the nonbonded parameters from a line of the parameter file

        Parameters
        ----------
        lines : list
            The lines of the parameter file
        idx : int
            The index of the line to start parsing
        """
        idx += 1
        while idx < len(lines):
            line = self._read_line(lines[idx])
            if len(line) == 0:
                idx += 1
                continue
            elif len(line) < 4:
                warnings.warn(f"[ignoring] line {idx+1} is not a valid NonBonded entry")
                idx += 1
                continue

            start = line[0]
            if start in self.__tags__:
                idx -= 1
                break

            atom = line[0]
            epsilon = float(line[2])
            min_radius = float(line[3])
            epsilon_14 = None
            min_radius_14 = None
            if len(line) == 7:
                epsilon_14 = float(line[5])
                min_radius_14 = float(line[6])

            nonbonded = _abstract.AbstractNonBonded(
                atom, epsilon, min_radius, epsilon_14, min_radius_14
            )
            self.add_nonbonded(nonbonded)

            idx += 1

        return idx

    def _parse_impropers(self, lines: list, idx: int):
        """
        Parses the improper parameters from a line of the parameter file

        Parameters
        ----------
        lines : list
            The lines of the parameter file
        idx : int
            The index of the line to start parsing
        """
        idx += 1
        while idx < len(lines):
            line = self._read_line(lines[idx])
            if len(line) == 0:
                idx += 1
                continue
            start = line[0]
            if start in self.__tags__:
                idx -= 1
                break
            if len(line) < 7:
                warnings.warn(f"[ignoring] line {idx+1} is not a valid Improper entry")
                idx += 1
                continue

            atoms = tuple(line[0:4])
            k = float(line[4])
            phi = float(line[6])
            improper = _abstract.AbstractImproper(*atoms, K=k, angle=phi)
            self.add_improper(improper)

            idx += 1

        return idx

    def _parse_mass(self, line):
        """
        Parses the mass from a line of the parameter file

        Parameters
        ----------
        line : list
            The line of the parameter file, split into a list
        """
        if len(line) < 4:
            warnings.warn(f"[ignoring] line is not a valid Mass entry")
            return
        atom_type = line[2]
        mass = float(line[3])
        self.add_mass(atom_type, mass)

    def _parse_dihedrals(self, lines: list, idx: int):
        """
        Parses the dihedral parameters from a line of the parameter file

        Parameters
        ----------
        lines : list
            The lines of the parameter file
        idx : int
            The index of the line to start parsing
        """
        idx += 1
        while idx < len(lines):
            line = self._read_line(lines[idx])
            if len(line) == 0:
                idx += 1
                continue
            start = line[0]
            if start in self.__tags__:
                idx -= 1
                break
            if len(line) < 5:
                warnings.warn(f"[ignoring] line {idx+1} is not a valid Dihedral entry")
                idx += 1
                continue

            _attrs = {}
            if len(line) >= 6:
                _attrs["multiplicity"] = int(line[5])

            if len(line) == 7:
                _attrs["angle"] = float(line[6])

            dihedral = _abstract.AbstractDihedral(
                line[0],
                line[1],
                line[2],
                line[3],
                K=float(line[4]),
                **_attrs,
            )

            self.add_dihedral(dihedral)
            idx += 1

        return idx

    def _parse_angles(self, lines: list, idx: int):
        """
        Parses the angle parameters from a line of the parameter file

        Parameters
        ----------
        lines : list
            The lines of the parameter file
        idx : int
            The index of the line to start parsing
        """
        idx += 1
        while idx < len(lines):
            line = self._read_line(lines[idx])
            if len(line) == 0:
                idx += 1
                continue
            start = line[0]
            if start in self.__tags__:
                idx -= 1
                break
            if len(line) < 5:
                warnings.warn(f"[ignoring] line {idx+1} is not a valid Angle entry")
                idx += 1
                continue

            _urey_bradley = {}
            if len(line) == 7:
                _urey_bradley["urey_bradley_k"] = float(line[5])
                _urey_bradley["urey_bradley_length"] = float(line[6])

            _angle = _abstract.AbstractAngle(
                line[0],
                line[1],
                line[2],
                K=float(line[3]),
                angle=float(line[4]),
                **_urey_bradley,
            )
            self._dict["angles"][_angle.atoms] = _angle

            idx += 1

        return idx

    def _parse_bonds(self, lines: list, idx: int):
        """
        Parses the bond parameters from a line of the parameter file

        Parameters
        ----------
        lines : list
            The lines of the parameter file
        idx : int
            The index of the line to start parsing
        """
        idx += 1
        while idx < len(lines):
            line = self._read_line(lines[idx])
            if len(line) == 0:
                idx += 1
                continue
            start = line[0]
            if start in self.__tags__:
                idx -= 1
                break
            if len(line) < 4:
                warnings.warn(f"[ignoring] line {idx+1} is not a valid Bond entry")
                idx += 1
                continue

            _bond = _abstract.AbstractBond(
                line[0],
                line[1],
                K=float(line[2]),
                length=float(line[3]),
            )
            self.add_bond(_bond)

            idx += 1

        return idx

    def _vet_load(self, _data):
        if not isinstance(_data, dict):
            raise TypeError("The file must contain a dictionary object")
        if "bonds" not in _data.keys() or "angles" not in _data.keys():
            raise KeyError("The dictionary must contain 'bonds' and 'angles' keys")


__all__ = [CHARMMTopology, CHARMMParameters]


if __name__ == "__main__":
    _carbs = "/Users/noahhk/GIT/glycosylator/support/toppar_charmm/carbohydrates.rtf"
    _top = CHARMMTopology.from_file(_carbs)
    print(_top)

    _carbs = "/Users/noahhk/GIT/glycosylator/support/toppar_charmm/carbohydrates.prm"
    _prm = CHARMMParameters.from_file(_carbs)
    print(_prm)

    from glycosylator.utils.defaults import (
        DEFAULT_CHARMM_TOPOLOGY_FILE,
        DEFAULT_CHARMM_PARAMETERS_FILE,
    )

    _save_to = "/Users/noahhk/GIT/glycosylator/glycosylator/resources/"
    _top.save(_save_to + os.path.basename(DEFAULT_CHARMM_TOPOLOGY_FILE))
    _prm.save(_save_to + os.path.basename(DEFAULT_CHARMM_PARAMETERS_FILE))
