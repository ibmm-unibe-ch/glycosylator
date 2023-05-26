"""
This script generates a dataset of images describing the connectivity and labelling of a molecule.

Encoding
--------
Each image is a NxNx2 matrix, where N is the number of atoms in the molecule. The rows and columns always mimic the atoms in the molecule
in the order they appear in the structure. 

The first dimension is the element of the atom, as a diagonal matrix that assignes the atomic number (1=H, 8=O, 6=C etc.).
The second dimension describes the connectivity between the atoms as a symmetric binary amtrix. The value 1 means that the atoms are
connected, and 0 means that they are not.

The desired output images are Nx37x4 matrices where each row represents one of the atoms in a molecule and the columns represent one of the
available characters for atom labelling: ABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789'. The third dimension is a binary matrix that indicates
whether the atom is labelled with the character in the column or not. Hence, a maximal atom label length of 4 characters is assumed.

Output
------
The generated output datafile is a pickle file that contains two lists of numpy arrays: inputs and outputs. The 'inputs' list contains the input images and the
'outputs' list contains the desired output images. The lists are wrapped in a dictionary with the keys 'inputs' and 'outputs'.
"""

__alphabet__ = "ABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789'"

import argparse
import logging
import numpy as np
import pandas as pd
import glycosylator as gl
import periodictable as pt
from alive_progress import alive_bar
import pickle
from datetime import datetime


def setup():
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawTextHelpFormatter
    )
    parser.add_argument(
        "-c",
        "--compounds",
        type=str,
        default=None,
        help="Path to a compounds mmCIF file. By default the pre-loaded compounds from Glycosylator are used.",
    )
    parser.add_argument(
        "-o",
        "--output",
        type=str,
        default=None,
        help="Path to the output file. By default a 'dataset-{timestamp}'.npz file is created in the current directory.",
    )
    parser.add_argument(
        "-v", "--verbose", action="store_true", help="Print progress to stdout."
    )
    return parser.parse_args()


def main(args):
    compounds = load_compounds(args.compounds)
    molecules = make_molecules(compounds)
    images = []
    labels = []
    with _make_bar(args.verbose, len(compounds)) as bar:
        for molecule in molecules:
            image, label = encode(molecule)
            images.append(image)
            labels.append(label)
            bar()
    final = {"inputs": images, "outputs": labels}
    if args.output is None:
        args.output = f"dataset-{datetime.now().strftime('%Y%m%d%H%M')}.pkl"
    with open(args.output, "wb") as f:
        pickle.dump(final, f)
    print(f"Saved dataset to {args.output}")


def load_compounds(compounds_path):
    if compounds_path is None:
        compounds = gl.get_default_compounds()
    else:
        compounds = gl.resources.PDBECompounds.from_file(compounds_path)
    return compounds


def make_molecules(compounds):
    for compound in compounds.ids:
        try:
            yield compounds.get(compound)
        except Exception as e:
            logging.warning(f"Could not load: {compound} (reason: {e})")
            continue


def encode(molecule):
    image_data, image_labels = _blank_images(molecule)
    image_data = _encode_elements(image_data, molecule)
    image_data = _encode_bonds(image_data, molecule)
    image_labels = _encode_labels(image_labels, molecule)
    return image_data, image_labels


def atom_number(atom):
    if len(atom.element) > 1:
        element = atom.element.title()
        return pt.elements.symbol(element).number
    return pt.elements.symbol(atom.element).number


def _blank_images(molecule):
    n_atoms = len(molecule.atoms)
    image_data = np.zeros((n_atoms, n_atoms, 2))
    image_labels = np.zeros((n_atoms, len(__alphabet__), 4))
    return image_data, image_labels


def _encode_elements(image, molecule):
    for i, atom in enumerate(molecule.get_atoms()):
        image[i, i, 0] = atom_number(atom)
    return image


def _encode_bonds(image, molecule):
    for bond in molecule.bonds:
        i = bond[0].serial_number - 1
        j = bond[1].serial_number - 1
        image[i, j, 1] = 1
        image[j, i, 1] = 1
    return image


def _encode_labels(image, molecule):
    for i, atom in enumerate(molecule.get_atoms()):
        for j, char in enumerate(atom.id):
            k = __alphabet__.index(char)
            image[i, k, j] = 1
        return image


def _make_bar(verbose, *args, **kwargs):
    if verbose:
        return alive_bar(*args, **kwargs)
    else:
        return _DummyBar(*args, **kwargs)


class _DummyBar:
    def __init__(self, *args, **kwargs):
        pass

    def __call__(self, *args, **kwargs):
        pass

    def __enter__(self):
        pass

    def __exit__(self, *args, **kwargs):
        pass


if __name__ == "__main__":
    args = setup()
    main(args)
