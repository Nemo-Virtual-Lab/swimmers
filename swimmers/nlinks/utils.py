import matlab
import numpy as np


def rotation_x(theta):
    """Generate a rotation matrix for a rotation around the x-axis by angle theta.

    Args:
        theta (float): The rotation angle in radians.
    """
    return np.array([[1, 0, 0],
                     [0, np.cos(theta), -np.sin(theta)],
                     [0, np.sin(theta), np.cos(theta)]])

def rotation_y(theta):
    """Generate a rotation matrix for a rotation around the y-axis by angle theta.

    Args:
        theta (float): The rotation angle in radians.
    """
    return np.array([[np.cos(theta), 0, np.sin(theta)],
                     [0, 1, 0],
                     [-np.sin(theta), 0, np.cos(theta)]])

def rotation_z(theta):
    """Generate a rotation matrix for a rotation around the z-axis by angle theta.

    Args:
        theta (float): The rotation angle in radians.
    """
    return np.array([[np.cos(theta), -np.sin(theta), 0],
                     [np.sin(theta), np.cos(theta), 0],
                     [0, 0, 1]])

def vector_to_euler_angles(v):
    """Convert a 3D vector to Euler angles (phi_y, phi_z).

    Args:
        v (np.array): The input 3D vector.

    Returns:
        np.array: The Euler angles (phi_y, phi_z).
    """
    normalized_v = v / np.linalg.norm(v)
    phi_y = np.arctan2(-normalized_v[2], np.sqrt(normalized_v[0]**2 + normalized_v[1]**2))
    phi_z = np.arctan2(normalized_v[1], normalized_v[0])

    return np.array([phi_y, phi_z])


def numpy_to_matlab_double(arr):
    """Convert a NumPy array to a MATLAB double array.
    
    Args:
        arr (np.array): The input NumPy array.
    """
    return matlab.double(arr.tolist())
