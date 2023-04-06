"""
Functions to convert different data formats
"""

# raise DeprecationWarning("This module is currently dropped...")

import os
import numpy as np

from periodictable import elements

import Bio.PDB as bio
from openbabel import pybel

import glycosylator.utils.auxiliary as aux


class PybelBioPythonConverter:
    """
    Convert Pybel (openbabel) data structures to Biopython
    """

    __element_counts__ = {}
    _current_residue = None

    def convert(self, obj, full_structure: bool = False):
        """
        Convert between Biopython and Pybel (openbabel) data structures

        Parameters
        ----------
        obj : object
            The object to convert. The following conversions are supported:
            - A pybel `Molecule` -> biopython `Structure`
            - A pybel `Residue` -> biopython `Residue`
            - A pybel `Atom` -> biopython `Atom`
            - A pybel `Residue` -> biopython `Structure` (if `full_structure=True`)
            - A pybel `Atom` -> biopython `Structure` (if `full_structure=True`)

        full_structure : bool, optional
            If `True` a new full structure object is made instead of simply
            the hierarchical class equivalent.

        Returns
        -------
        object
            The converted object
        """
        if isinstance(obj, pybel.Molecule):
            return self.pybel_molecule_to_biopython(obj)

        elif isinstance(obj, pybel.Residue):
            residue = self.pybel_residue_to_biopython(obj)
            if full_structure:
                struct = self._new_biopython_structure()
                struct[0]["A"].add(residue)
                return struct
            return residue

        elif isinstance(obj, pybel.Atom):
            atom = self.pybel_atom_to_biopython(obj)
            if full_structure:
                struct = self._new_biopython_structure()
                residue = bio.Residue.Residue((" ", 1, " "), "UNL", " ")
                residue.add(atom)
                struct[0]["A"].add(atom.get_parent())
                return struct
            return atom

        else:
            raise ValueError(f"Cannot convert object of type {type(obj)}")

    @staticmethod
    def _new_biopython_structure():
        """
        Create a new biopython structure object with a Model and Chain
        """
        struct = bio.Structure.Structure("New Molecule")
        struct.add(bio.Model.Model(0))
        struct[0].add(bio.Chain.Chain("A"))
        return struct

    # @staticmethod
    # def write_tmp_pdb(obj):
    #     """
    #     Write a temporary pdb file from a source object

    #     Parameters
    #     ----------
    #     obj : object
    #         The object to convert.

    #     Returns
    #     -------
    #     str
    #         The path to the temporary pdb file
    #     """
    #     tmpfile = mktemp(suffix="glyco.tmp.pdb")
    #     if is_biopython(obj):
    #         PybelBioPythonConverter._write_pdb_from_biopython(obj, tmpfile)
    #     elif is_pybel(obj):
    #         PybelBioPythonConverter._write_pdb_from_pybel(obj, tmpfile)
    #     else:
    #         raise ValueError(f"Cannot convert object of type {type(obj)} to pdb")
    #     return tmpfile

    # @staticmethod
    # def _write_pdb_from_biopython(obj, filename):
    #     """
    #     Write a temporary pdb file from a biopython object
    #     """

    #     io = bio.PDBIO()

    #     if not isinstance(obj, bio.Structure.Structure):
    #         _obj = bio.Structure.Structure()

    #         if not isinstance(obj, bio.Model.Model):
    #             _obj.add(bio.Model.Model(0))

    #             if not isinstance(obj, bio.Chain.Chain):
    #                 _obj[0].add(bio.Chain.Chain("A"))

    #                 if isinstance(obj, bio.Residue.Residue):
    #                     _obj[0]["A"].add(obj)

    #                 elif isinstance(obj, bio.Atom.Atom):
    #                     _obj[0]["A"].add(bio.Residue.Residue())
    #                     _obj[0]["A"][0].add(obj)

    #                 else:
    #                     raise ValueError(f"Cannot convert object of type {type(obj)} to pdb")

    #             else:
    #                 _obj[0].add(obj)

    #         else:
    #             _obj.add(obj)

    #     else:
    #         _obj = obj

    #     io.set_structure(_obj)
    #     io.save(filename)

    # @staticmethod
    # def _write_pdb_from_pybel(obj, filename):
    #     """
    #     Write a temporary pdb file from a pybel object
    #     """

    #     io = bio.PDBIO()

    #     if not isinstance(obj, pybel.Molecule):
    #         _obj = pybel.Molecule()

    #         if not isinstance(obj, pybel.Residue):
    #             _obj.add(pybel.Residue())

    #             if isinstance(obj, pybel.Atom):
    #                 _obj[0].add(obj)

    #             else:
    #                 raise ValueError(f"Cannot convert object of type {type(obj)} to pdb")

    #         else:
    #             _obj.add(obj)

    #     else:
    #         _obj = obj

    #     io.set_structure(_obj)
    #     io.save(filename)

    def pybel_atom_to_biopython(self, obj):
        """
        Convert a pybel `Atom` to a biopython `Atom`

        Parameters
        ----------
        obj : pybel.Atom
            The pybel `Atom` to convert

        Returns
        -------
        bio.Atom.Atom
            The converted biopython `Atom`
        """

        # just for testing purposes...
        if obj.atomicnum == 67:
            element = "H"
        else:
            element = str(elements[obj.atomicnum])

        # # due to the issue of Pybel reading the CHARMM force-field
        # # designation of hydrogens by ids relating to their connectivity
        # # e.g. HO2 = Hydrogen at Oxygen2, the element is falsely inferred
        # # as Holmium (Ho) instead of H, so we need to manually correct for
        # # this - thereby, of course, destroying glycosylator capacity to work
        # # with Holmium, should anybody actually wish to use it...
        # if element == "Ho":
        #     element = "H"

        if not element in self.__element_counts__:
            self.__element_counts__[element] = {None: 0}

        self.__element_counts__[element].setdefault(self._current_residue, 0)
        self.__element_counts__[element][self._current_residue] += 1
        name = f"{element}{self.__element_counts__[element][self._current_residue]}"

        new = bio.Atom.Atom(
            name=name,
            coord=np.asarray(obj.coords),
            serial_number=obj.idx,
            bfactor=0.0,
            occupancy=0.0,
            altloc=" ",
            fullname=name,
            element=element,
        )
        return new

    def pybel_residue_to_biopython(self, obj):
        """
        Convert a pybel `Residue` to a biopython `Residue`

        Parameters
        ----------
        obj : pybel.Residue
            The pybel `Residue` to convert

        Returns
        -------
        bio.Residue.Residue
            The converted biopython `Residue`
        """
        self._current_residue = (obj.idx, obj.name)
        new = bio.Residue.Residue(
            id=(f"H_{obj.name}", obj.idx + 1, " "), resname=obj.name, segid=" "
        )
        for atom in obj.atoms:
            atom = self.pybel_atom_to_biopython(atom)
            new.add(atom)

        self._current_residue = None
        return new

    def pybel_molecule_to_biopython(self, obj):
        """
        Convert a pybel `Molecule` to a biopython `Structure`

        Parameters
        ----------
        obj : pybel.Molecule
            The pybel `Molecule` to convert

        Returns
        -------
        bio.Structure.Structure
            The converted biopython `Structure`
        """
        if os.path.exists(obj.title):
            id = aux.filename_to_id(obj.title)
        else:
            id = obj.title

        new = bio.Structure.Structure(id)
        new.add(bio.Model.Model(0))
        new[0].add(bio.Chain.Chain("A"))

        if len(obj.residues) == 0:
            residue = bio.Residue.Residue(("H_UNL", 1, " "), "UNL", " ")
            new[0]["A"].add(residue)
            for atom in obj.atoms:
                residue.add(self.pybel_atom_to_biopython(atom))
        else:
            for residue in obj.residues:
                residue = self.pybel_residue_to_biopython(residue)
                new[0]["A"].add(residue)

        return new


def is_biopython(obj):
    """
    Check if an object is a biopython object

    Parameters
    ----------
    obj : object
        The object to check

    Returns
    -------
    bool
        `True` if the object is a biopython object, `False` otherwise
    """
    _valids = {
        bio.Structure.Structure,
        bio.Model.Model,
        bio.Chain.Chain,
        bio.Residue.Residue,
        bio.Atom.Atom,
    }
    return any([isinstance(obj, _valid) for _valid in _valids])


def is_pybel(obj):
    """
    Check if an object is a pybel object

    Parameters
    ----------
    obj : object
        The object to check

    Returns
    -------
    bool
        `True` if the object is a pybel object, `False` otherwise
    """
    _valids = {
        pybel.Molecule,
        pybel.Residue,
        pybel.Atom,
    }
    return any([isinstance(obj, _valid) for _valid in _valids])


if __name__ == "__main__":

    # MANNOSE = "support/examples/MAN.pdb"
    # _pybel = pybel.readfile("pdb", MANNOSE)
    # _pybel = next(_pybel)

    # converter = PybelBioPythonConverter()

    # _biopython = converter.pybel_molecule_to_biopython(_pybel)
    # print(_biopython)

    mol = pybel.readstring("smi", "OCC1OC(O)C(C(C1O)O)O")
    mol.addh()
    mol.make3D()

    converter = PybelBioPythonConverter()
    _biopython = converter.pybel_molecule_to_biopython(mol)
    print(_biopython)
