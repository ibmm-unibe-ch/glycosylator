"""
Icons for drawing glycan structures in 2D
"""

import os
import PIL

DIR = os.path.dirname(__file__)

__loaded_icons__ = {}


def get_icon_path(name: str) -> str:
    """
    Get the path to an icon PNG file

    Parameters
    ----------
    name : str

    Returns
    -------
    str
    """
    return os.path.join(DIR, name + ".png")


def get_icon(name: str) -> PIL.Image:
    """
    Get an icon by name

    Parameters
    ----------
    name : str
        The name of the icon

    Returns
    -------
    PIL.Image
        The icon
    """
    if name not in __loaded_icons__:
        __loaded_icons__[name] = PIL.Image.open(get_icon_path(name))
    return __loaded_icons__[name]
