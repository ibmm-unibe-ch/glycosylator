"""
Basic structure related functions
"""

import numpy as np


def compute_angle(atom1, atom2, atom3):
    """
    Compute the angle between three atoms.

    Parameters
    ----------
    atom1 : Bio.PDB.Atom
        The first atom
    atom2 : Bio.PDB.Atom
        The second atom
    atom3 : Bio.PDB.Atom
        The third atom

    Returns
    -------
    angle : float
        The angle between the three atoms in degrees
    """
    a = atom1.coord - atom2.coord
    b = atom3.coord - atom2.coord
    return np.degrees(np.arccos(np.dot(a, b) / (np.linalg.norm(a) * np.linalg.norm(b))))


def compute_dihedral(atom1, atom2, atom3, atom4):
    """
    Compute the dihedral angle between four atoms.

    Parameters
    ----------
    atom1 : Bio.PDB.Atom
        The first atom
    atom2 : Bio.PDB.Atom
        The second atom
    atom3 : Bio.PDB.Atom
        The third atom
    atom4 : Bio.PDB.Atom
        The fourth atom

    Returns
    -------
    dihedral : float
        The dihedral angle between the four atoms in degrees
    """
    ab = atom1.coord - atom2.coord
    bc = atom3.coord - atom2.coord
    cd = atom4.coord - atom3.coord

    # normalize bc so that it does not influence magnitude of vector
    # rejections that come next
    bc /= np.linalg.norm(bc)

    # vector rejections
    v = ab - np.dot(ab, bc) * bc
    w = cd - np.dot(cd, bc) * bc

    # angle between v and w in radians
    x = np.dot(v, w)
    y = np.dot(np.cross(bc, v), w)
    return np.degrees(np.arctan2(y, x))


def compute_torsional(atom1, atom2, atom3, atom4):
    """
    Compute the torsional angle between four atoms.

    Parameters
    ----------
    atom1 : Bio.PDB.Atom
        The first atom
    atom2 : Bio.PDB.Atom
        The second atom
    atom3 : Bio.PDB.Atom
        The third atom
    atom4 : Bio.PDB.Atom
        The fourth atom

    Returns
    -------
    torsional : float
        The torsional angle between the four atoms in degrees
    """
    ab = atom1.coord - atom2.coord
    bc = atom3.coord - atom2.coord
    cd = atom4.coord - atom3.coord

    # normalize b so that it does not influence magnitude of vector
    # rejections that come next
    bc /= np.linalg.norm(bc)

    # vector rejections
    v = ab - np.dot(ab, bc) * bc
    w = cd - np.dot(cd, bc) * bc

    # angle between v and w in radians
    x = np.dot(v, w)
    y = np.dot(np.cross(bc, v), w)
    return np.degrees(np.arctan2(y, x))



def center_of_gravity(masses, coords):
    """
    Compute the center of gravity of a molecule.

    Parameters
    ----------
    masses : array-like
        The masses of the atoms as an nx1 vector
    coords : array-like
        The coordinates of the atoms as an nx3 array

    Returns
    -------
    cog : array-like
        The center of gravity
    """
    return np.average(coords, axis=0, weights=masses)


def _rotation_matrix(axis, angle):
    """
    Compute the rotation matrix about an arbitrary axis in 3D

    Source
    ------
    Stackoverflow thread: http://stackoverflow.com/questions/6802577/python-rotation-of-3d-vector

    Parameters
    ----------
    axis : array-like
        The axis to rotate around
    angle : float
        The angle to rotate by (in radians)

    Returns
    -------
    rotation_matrix : array-like
        The rotation matrix
    """
    # axis = np.asarray(axis)
    # angle = np.asarray(angle)
    axis = axis / np.sqrt(np.dot(axis, axis))
    a = np.cos(angle / 2.0)
    b, c, d = -axis * np.sin(angle / 2.0)
    aa, bb, cc, dd = a * a, b * b, c * c, d * d
    bc, ad, ac, ab, bd, cd = b * c, a * d, a * c, a * b, b * d, c * d
    return np.array(
        [
            [aa + bb - cc - dd, 2 * (bc + ad), 2 * (bd - ac)],
            [2 * (bc - ad), aa + cc - bb - dd, 2 * (cd + ab)],
            [2 * (bd + ac), 2 * (cd - ab), aa + dd - bb - cc],
        ]
    )


def _IC_to_xyz(a, b, c, anchor, r, theta, dihedral):
    """
    compute the coordinates of a fourth atom from a proper internal coordinate
    system and the other three atom coordinates.

    Parameters
    ----------
    a, b, c : np.ndarray
        coordinates of the other three atoms
    anchor : np.ndarray
        coordinates of the anchor atom relative to which the new coordinate should be calculated
    r : float
        bond length of the new atom relative to the anchor
    theta : float
        bond angle between the new atom and its plane partners
    dihedral : float
        dihedral angle of the internal coordinates
    """
    ab = b - a
    bc = c - b

    # compute normalized bond vectors for available atoms
    ab /= np.linalg.norm(ab)
    bc /= np.linalg.norm(bc)

    # compute plane vector for atoms 1-2-3
    plane_abc = np.cross(ab, bc)
    plane_abc /= np.linalg.norm(plane_abc)

    # rotate the plane vector around the middle bond (2-3) to get the plane 2-3-4
    _rot = _rotation_matrix(bc, dihedral)
    plane_bcd = np.dot(_rot, plane_abc)
    plane_bcd /= np.linalg.norm(plane_bcd)

    # rotate the middle bond around the new plane
    _rot = _rotation_matrix(plane_bcd, theta)
    cd = np.dot(_rot, bc)
    cd /= np.linalg.norm(cd)

    # compute the coordinates of the fourth atom
    d = anchor + r * cd
    return d
