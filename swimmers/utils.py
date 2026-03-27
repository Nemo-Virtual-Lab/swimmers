import sys
import feelpp.core as fppc
import feelpp.toolboxes.core as tb
import numpy as np


def init_environment(name: str, swimmer: str):
    """Initialize the Feel++ environment for a given configuration file.
    
    Args:
        name (str): Path where solution files are stored.
        swimmer (str): type of swimmer.

    Returns:
        fppc.Environment: The initialized Feel++ environment.
    """
    sys.argv = [name]
    if swimmer == "three_sphere":
        e = fppc.Environment(sys.argv, opts=tb.toolboxes_options("fluid"),config=fppc.globalRepository(name))
    elif swimmer == "magneto":
        e = fppc.Environment(sys.argv, opts=tb.toolboxes_options("fsi"),config=fppc.globalRepository(name))
    return e


class Box:

    def __init__(self, point: np.array, lengths: np.array) -> None:
        """Class that defines a box geometry in 2D.

        Args:
            point (np.array): Array with the x and y coordinates of the bottom left point of the box : [x,y]. Dimensions: (2,)
            lengths (np.array): Array with the lengths of the box : [length_x, length_y]. Dimensions: (2,)
        """

        self.point = point
        self.lengths = lengths

