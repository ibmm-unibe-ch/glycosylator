from itertools import izip
import math
import numpy as np
import re
import os

def readLinesFromFile(fileName):
    """Reads all lines in a file
    Parameters:
        fileName: path to file
    Returns:
        lines: list with all the lines in a file
    """
    file = open(fileName,'r')                  # open the file
    lines = file.readlines()                  # read all the lines in the file to the list "lines"
    file.close()                                  # close the file
    return lines

def topological_sort(unsorted_graph):
    """Topological sorting of a graph 
    Parameters:
        unsorted_graph: dictionary representation of a graph
    Returns:
        sorted_graph: list of nodes and corresponding edges
    """
    sorted_graph = []
    #sort graph
    while unsorted_graph:
        acyclic = False
        for node, edges in unsorted_graph.items():
            for edge in edges:
                if edge in unsorted_graph:
                    break
            else:
                acyclic = True
                del unsorted_graph[node]
                sorted_graph.append((node, edges))

        if not acyclic:
            print "WARNING! Cyclique dependency occurred in ICs. Impossible to build residue"
            print unsorted_graph
            print sorted_graph
            return ''
            break
    return sorted_graph[::-1]    

def pairwise(iterable):
    "s -> (s0, s1), (s2, s3), (s4, s5), ..."
    a = iter(iterable)
    return izip(a, a)

def rotation_matrix(axis, theta):
    '''Computes the rotation matrix about an arbitrary axis in 3D
    Code from: http://stackoverflow.com/questions/6802577/python-rotation-of-3d-vector
    Parameters:
        axis: axis
        theta: rotation angle
    Return: 
        rotation matrix

    '''

    axis = np.asarray(axis)
    theta = np.asarray(theta)
    axis = axis/math.sqrt(np.dot(axis, axis))
    a = math.cos(theta/2.0)
    b, c, d = -axis*math.sin(theta/2.0)
    aa, bb, cc, dd = a*a, b*b, c*c, d*d
    bc, ad, ac, ab, bd, cd = b*c, a*d, a*c, a*b, b*d, c*d
    return np.array([[aa+bb-cc-dd, 2*(bc+ad), 2*(bd-ac)],
                     [2*(bc-ad), aa+cc-bb-dd, 2*(cd+ab)],
                     [2*(bd+ac), 2*(cd-ab), aa+dd-bb-cc]])

def rotation_matrix2(angle, direction, point=None):
    """Return matrix to rotate about axis defined by point and direction.

    """
    sina = math.sin(angle)
    cosa = math.cos(angle)
    direction = np.array(direction[:3])
    direction = direction/math.sqrt(np.dot(direction, direction))
    # rotation matrix around unit vector
    R = np.diag([cosa, cosa, cosa])
    R += np.outer(direction, direction) * (1.0 - cosa)
    direction *= sina
    R += np.array([[ 0.0,         -direction[2],  direction[1]],
                   [ direction[2], 0.0,          -direction[0]],
                   [-direction[1], direction[0],  0.0]])
    M = np.identity(4)
    M[:3, :3] = R
    if point is not None:
        # rotation not around origin
        point = np.array(point[:3], dtype=np.float64, copy=False)
        M[:3, 3] = point - np.dot(R, point)
    return M


    
def alphanum_sort(l):
    """Alphanumerical sort of a list from
    https://arcpy.wordpress.com/2012/05/11/sorting-alphanumeric-strings-in-python/
    Parameter:
        l: list
    Returns:
        alphanumerically sorted list
    """
    convert = lambda text: int(text) if text.isdigit() else text
    alphanum_key = lambda key: [convert(c) for c in re.split('([0-9]+)', key)]
    return sorted(l, key = alphanum_key)

aaa2a = {'CYS': 'C', 'ASP': 'D', 'SER': 'S', 'GLN': 'Q', 'LYS': 'K',
     'ILE': 'I', 'PRO': 'P', 'THR': 'T', 'PHE': 'F', 'ASN': 'N', 
     'GLY': 'G', 'HIS': 'H', 'LEU': 'L', 'ARG': 'R', 'TRP': 'W', 
     'ALA': 'A', 'VAL':'V', 'GLU': 'E', 'TYR': 'Y', 'MET': 'M'}

GLYCOSYLATOR_PATH = os.path.dirname(os.path.realpath(__file__))
#SELF_BIN = os.path.dirname(os.path.realpath(sys.argv[0]))
#sys.path.insert(0, SELF_BIN + '/support')
