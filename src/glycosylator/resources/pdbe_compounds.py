"""
Parsers to extract information from the PDBE compound database.
"""

import warnings
import numpy as np
import pdbecif.mmcif_tools as mmcif_tools
import pickle

import Bio.PDB as bio
from Bio import SVDSuperimposer

import glycosylator.utils.auxiliary as aux
import glycosylator.core.molecule as molecule
import glycosylator.structural.infer as infer

__to_float__ = {
    "_chem_comp": set(
        (
            "formula_weight",
            "pdbx_formal_charge",
        )
    ),
}

__to_delete__ = {
    "_chem_comp": set(
        (
            "pdbx_type",
            "mon_nstd_parent_comp_id",
            "pdbx_synonyms",
            "pdbx_initial_date",
            "pdbx_modified_date",
            "pdbx_release_status",
            "pdbx_ambiguous_flag",
            "pdbx_replaced_by",
            "pdbx_replaces",
            "pdbx_processing_site",
            "pdbx_model_coordinates_details",
            "pdbx_model_coordinates_db_code",
            "pdbx_ideal_coordinates_details",
            "pdbx_subcomponent_list",
            "pdbx_model_coordinates_missing_flag",
            "pdbx_ideal_coordinates_missing_flag",
        )
    ),
}

# a list of categories to ignore while parsing
__to_ignore__ = [
    "_pdbx_chem_comp_audit",
    "_pdbx_chem_comp_feature",
    "_pdbx_chem_comp_identifier",
]


class PDBECompounds:
    """
    This class wraps a dictionary of PDBE compounds.
    It facilitates easy data access and searchability for compounds.

    Parameters
    ----------
    compounds : dict
        A dictionary of PDBE compounds.
    """

    def __init__(self, compounds: dict) -> None:

        self._compounds = {k: None for k in compounds.keys()}
        self._pdb = dict(self._compounds)
        self._setup_dictionaries(compounds)

        self._filename = None
        self._orig = compounds

    @classmethod
    def from_file(cls, filename: str) -> "PDBECompounds":
        """
        Create a PDBECompounds object from a `cif` file.

        Parameters
        ----------
        filename : str
            The path to the file.

        Returns
        -------
        PDBECompounds
            The PDBECompounds object.
        """
        _dict = mmcif_tools.MMCIF2Dict()
        _dict = _dict.parse(filename, ignoreCategories=__to_ignore__)
        new = cls(_dict)
        new._filename = filename
        return new

    @classmethod
    def load(self, filename: str) -> "PDBECompounds":
        """
        Load a PDBECompounds object from a `pickle` file.

        Parameters
        ----------
        filename : str
            The path to the file.

        Returns
        -------
        PDBECompounds
            The PDBECompounds object.
        """
        with open(filename, "rb") as f:
            return pickle.load(f)

    @property
    def ids(self) -> list:
        """
        Get a list of all compound ids.

        Returns
        -------
        list
            A list of all compound ids.
        """
        return list(self._compounds.keys())

    @property
    def formulas(self) -> list:
        """
        Get a list of all compound formulas.

        Returns
        -------
        list
            A list of all compound formulas.
        """
        return [i["formula"] for i in self._compounds.values()]

    def save(self, filename: str = None) -> None:
        """
        Save the PDBECompounds object to a `pickle` file.

        Parameters
        ----------
        filename : str
            The path to the file. By default the same file used to load the object is used.
            If no file is available, a `ValueError` is raised.
        """
        if filename is None:
            if self._filename is None:
                raise ValueError("No filename specified.")
            filename = self._filename

        if not filename.endswith(".pkl"):
            filename = aux.change_suffix(filename, ".pkl")

        self._filename = filename

        with open(filename, "wb") as f:
            pickle.dump(self, f)

    def get(
        self,
        query: str,
        by: str = "id",
        return_type: str = "molecule",
    ):
        """
        Get a compound that matches the given criteria.

        Parameters
        ----------
        query : str
            The query to search for.

        by : str, optional
            The type of query, by default "id". Possible values are:
            - "id": search by compound id
            - "name": search by compound name (must match any available synonym exactly)
            - "formula": search by compound formula
            - "smiles": search by compound smiles (this also works for inchi)

        return_type : str, optional
            The type of object to return, by default "molecule". Possible values are:
            - "molecule": return a glycosylator `Molecule` object
            - "dict": return a dictionary of the compound data
            - "structure": return a biopython `Structure` object
            - "residue": return a biopython `Residue` object

        Returns
        -------
        object
            The object that matches the given criteria.
        """
        _dict = self._get(query, by)
        if len(_dict.keys()) == 0:
            raise ValueError(f"No compound found for query '{query}'.")
        elif len(_dict.keys()) > 1:
            return [self.get(i, "id", return_type) for i in _dict.keys()]

        _dict = list(_dict.values())[0]
        if return_type == "molecule":
            return self._make_molecule(_dict)
        elif return_type == "dict":
            return _dict
        elif return_type == "structure":
            return self._make_structure(_dict)
        elif return_type == "residue":
            res = self._make_residue(_dict)
            self._fill_residue(res, _dict)
            return res
        else:
            raise ValueError(f"Invalid return_type '{return_type}'.")

    def has_residue(self, query: str, by: str = "id") -> bool:
        """
        Check if a compound has a residue definition.

        Parameters
        ----------
        query : str
            The query to search for.

        by : str, optional
            The type of query, by default "id". Possible values are:
            - "id": search by compound id
            - "name": search by compound name (must match any available synonym exactly)
            - "formula": search by compound formula
            - "smiles": search by compound smiles (this also works for inchi)

        Returns
        -------
        bool
            True if the compound has a residue definition, False otherwise.
        """
        _dict = self._get(query, by)
        return len(_dict.keys()) > 0
    
    
    def translate_ids_3_to_1(self, ids: list) -> list:
        """
        Translate a list of 3-letter compound ids to 1-letter ids.
        Any ids that cannot be translated are returned as "X".

        Parameters
        ----------
        ids : list
            A list of 3-letter compound ids.

        Returns
        -------
        list
            A list of 1-letter compound ids.
        """
        new = []
        for i in ids:
            i = self.get(i, "id", "dict")
            i = i.get("one_letter_code")
            if i is None:
                i = "X"
            new.append(i)
        return new

    def translate_ids_1_to_3(self, ids: list) -> list:
        """
        Translate a list of 1-letter compound ids to 3-letter ids.
        Any unknown ids are replaced with "XXX".

        Parameters
        ----------
        ids : list
            A list of 1-letter compound ids.

        Returns
        -------
        list
            A list of 3-letter compound ids.
        """
        new = []
        for i in ids:
            i = self.get(i, "id", "dict")
            i = i.get("three_letter_code")
            if i is None:
                i = "XXX"
            new.append(i)
            return new

    def relabel_atoms(self, structure):
        """
        Relabel the atoms of a `Molecule` object to match the atom names of the given compound.

        Parameters
        ----------
        structure
            A `Molecule` object or biopython object holding at least one residue.

        Returns
        -------
        object
            The object with relabeled atoms.
        """
        imposer = SVDSuperimposer.SVDSuperimposer()
        for residue in structure.get_residues():

            # get a reference
            name = residue.resname
            ref = self.get(name, "id", "residue")
            if ref is None:
                warnings.warn(f"Could not find residue '{name}' in PDBECompounds.")
                continue

            residue_coords = np.asarray([atom.coord for atom in residue.get_atoms()])
            ref_coords = np.asarray([atom.coord for atom in ref.get_atoms()])
            imposer.set(ref_coords, residue_coords)
            imposer.run()

            # get the residue coordinates
            new_coords = imposer.get_transformed()

            ref_atom_ids = [atom.id for atom in ref.get_atoms()]
            for idx, atom in enumerate(list(residue.get_atoms())):
                _new = new_coords[idx]
                min_idx = np.argmin(np.linalg.norm(ref_coords - _new, axis=1))
                atom.id = ref_atom_ids[min_idx]

        return structure

    def _get(self, q: str, by: str = "id"):
        """
        Get a dictionary of compounds that match the given criteria.

        Parameters
        ----------
        q : str
            The query to search for.
        by : str, optional
            The type of query, by default "id". Possible values are:
            - "id": search by compound id
            - "name": search by compound name (must match any available synonym exactly)
            - "formula": search by compound formula
            - "smiles": search by compound smiles (this also works for inchi)

        Returns
        -------
        dict
            A dictionary of compounds that match the given criteria.
        """
        if by == "id":
            _q = self._compounds.get(q)
            if _q is None:
                return {}
            return {q: _q}
        elif by == "name":
            return {k: v for k, v in self._compounds.items() if q.lower() in v["names"]}
        elif by == "formula":
            return {
                k: v
                for k, v in self._compounds.items()
                if q.upper().replace(" ", "") == v["formula"]
            }
        elif by == "smiles":
            return {k: v for k, v in self._compounds.items() if q in v["descriptors"]}
        else:
            raise ValueError(f"Invalid search type: {by}")

    def _setup_dictionaries(self, data_dict):
        """
        Fill in the dictionaries with the appropriate data

        Parameters
        ----------
        data_dict : dict
            A dictionary of data from the cif file.
        """
        for key, value in data_dict.items():

            comp = value["_chem_comp"]

            for k in __to_delete__["_chem_comp"]:
                if k in comp:
                    comp.pop(k)
            for k in __to_float__["_chem_comp"]:
                comp[k] = float(comp[k])
            for k in comp:
                if comp[k] == "?":
                    comp[k] = None

            comp["descriptors"] = value["_pdbx_chem_comp_descriptor"]["descriptor"]

            comp["names"] = [comp["name"]]
            synonyms = value.get("_pdbx_chem_comp_synonyms", None)
            if synonyms:
                comp["names"].extend(synonyms["name"])
            comp["names"] = set(i.lower() for i in comp["names"])

            comp["formula"] = comp["formula"].replace(" ", "")

            self._compounds[key] = comp

            atoms = value["_chem_comp_atom"]
            bonds = value["_chem_comp_bond"]

            atoms["pdbx_model_Cartn_x_ideal"] = [
                i.replace("?", "NaN") for i in atoms["pdbx_model_Cartn_x_ideal"]
            ]
            atoms["pdbx_model_Cartn_y_ideal"] = [
                i.replace("?", "NaN") for i in atoms["pdbx_model_Cartn_y_ideal"]
            ]
            atoms["pdbx_model_Cartn_z_ideal"] = [
                i.replace("?", "NaN") for i in atoms["pdbx_model_Cartn_z_ideal"]
            ]
            atoms["charge"] = [i.replace("?", "NaN") for i in atoms["charge"]]

            pdb = {
                "atoms": {
                    "ids": atoms["pdbx_component_atom_id"],
                    "serials": np.array(atoms["pdbx_ordinal"], dtype=int),
                    "coords": np.array(
                        [
                            (float(i), float(j), float(k))
                            for i, j, k in zip(
                                atoms["pdbx_model_Cartn_x_ideal"],
                                atoms["pdbx_model_Cartn_y_ideal"],
                                atoms["pdbx_model_Cartn_z_ideal"],
                            )
                        ]
                    ),
                    "elements": atoms["type_symbol"],
                    "charges": np.array(atoms["charge"], dtype=float),
                },
                "bonds": [
                    (a, b) for a, b in zip(bonds["atom_id_1"], bonds["atom_id_2"])
                ],
            }
            self._pdb[key] = pdb

    def _make_molecule(self, compound: dict) -> molecule.Molecule:
        """
        Make a glycosylator Molecule from a compound.

        Parameters
        ----------
        compound : dict
            A dictionary of a compound.

        Returns
        -------
        Molecule
            A glycosylator Molecule.
        """
        struct = self._make_structure(compound)
        mol = molecule.Molecule(struct)
        pdb = self._pdb[compound["id"]]
        for bond in pdb["bonds"]:
            mol.add_bond(*bond)

        return mol

    def _make_structure(self, compound: dict) -> "bio.PDB.Structure":
        """
        Make a structure from a compound.

        Parameters
        ----------
        compound : dict
            A dictionary of a compound.

        Returns
        -------
        bio.PDB.Structure
            A structure.
        """
        if len(compound) == 0:
            return None

        struct = bio.Structure.Structure(compound["id"])
        model = bio.Model.Model(0)
        struct.add(model)
        chain = bio.Chain.Chain("A")
        model.add(chain)
        res = self._make_residue(compound)
        chain.add(res)
        # add atoms to residue (at the end to ensure that the residue's full id is propagated)
        self._fill_residue(res, compound)
        return struct

    def _fill_residue(self, residue, compound: dict):
        """
        Fill a residue with atoms from a compound.

        Parameters
        ----------
        residue : bio.PDB.Residue
            A residue.
        compound : dict
            A dictionary of a compound.
        """
        pdb = self._pdb.get(compound["id"], None)
        if pdb is None:
            raise ValueError("No pdb data for compound.")

        atoms = pdb["atoms"]
        for i in range(len(atoms["ids"])):
            atom = self._make_atom(
                atoms["ids"][i],
                atoms["serials"][i],
                atoms["coords"][i],
                atoms["elements"][i],
                atoms["charges"][i],
            )
            residue.add(atom)

    def _make_residue(self, compound: dict) -> "bio.PDB.Residue":
        """
        Make a residue from a compound.

        Parameters
        ----------
        compound : dict
            A dictionary of a compound.

        Returns
        -------
        bio.PDB.Residue
            A residue.
        """
        if len(compound) == 0:
            return None

        res = bio.Residue.Residue(("H_" + compound["id"], 1, " "), compound["id"], " ")

        return res

    def _make_atom(
        self, id: str, serial: int, coords: np.ndarray, element: str, charge: float
    ) -> "bio.PDB.Atom":
        """
        Make an atom.

        Parameters
        ----------
        id : str
            The atom id.
        serial : int
            The atom serial.
        coords : np.ndarray
            The atom coordinates.
        element : str
            The atom element.

        Returns
        -------
        bio.PDB.Atom
            An atom.
        """
        atom = bio.Atom.Atom(
            id,
            coord=coords,
            serial_number=serial,
            bfactor=0.0,
            occupancy=0.0,
            element=element,
            fullname=id,
            altloc=" ",
            pqr_charge=charge,
        )
        return atom

    def __getitem__(self, key):
        return self._compounds[key], self._pdb[key]


if __name__ == "__main__":

    import glycosylator.utils.defaults as defaults

    f = "support/pdbe_compounds/test.cif"
    compounds = PDBECompounds.from_file(f)

    from glycosylator.core import Molecule

    real = Molecule.from_compound("GLC")

    scrambled = Molecule.from_compound("GLC")

    counts = {"C": 0, "H": 0, "O": 0, "N": 0, "S": 0, "P": 0}
    for atom in scrambled.atoms:
        counts[atom.element] += 1
        atom.id = atom.element + str(counts[atom.element])

    print(scrambled.atoms)
    compounds.relabel_atoms(scrambled)
    print(scrambled.atoms)

    from glycosylator.utils.visual import MoleculeViewer3D

    v = MoleculeViewer3D(scrambled)
    v.show()
