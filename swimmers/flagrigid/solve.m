
from matlab import double
from matlab import int64

import numpy as np
import matplotlib.pyplot as plt
import mpl_toolkits.mplot3d as mplot3d
import swimmers.flagrigid.utils as utils


def write_msh_flagellum(eng, 
                        N, 
                        L, 
                        rad, 
                        step_length, 
                        height, 
                        rad_section, 
                        section_type, 
                        ke_type, 
                        filename) :
    """
    Write the .msh file of the flagellum using the matlab function write_msh_flagellum.

    Args:
        eng: The matlab engine instance.
        N: The number of segments of the flagellum.
        L: The length of the flagellum.
        rad: The radius of the flagellum.
        step_length: The step length of the flagellum.
        height: first estimate of the total height of the flagellum
        rad_section: The radius of the section of the flagellum.
        section_type: The type of the section of the flagellum (e.g., "circular").
        ke_type: The type of the kinetic energy (e.g., "schum" or "pt").
        filename: The filename to save the .msh file (without extension).
    """
    
    eng.write_msh_flagellum(N, L, rad, step_length, height, rad_section, section_type, ke_type, filename, nargout=0)

def get_xF1(eng, mesh_head, alpha, beta) :

    """
    Get the xF1 value for the flagellum.

    Args:
        eng: The matlab engine instance.
        mesh_head: The mesh of the head of the swimmer.
        alpha: The angle alpha of the flagellum.
        beta: The angle beta of the flagellum.

    Returns:
        The xF1 value for the flagellum.
    """  
    return eng.get_xF1(double(np.array(mesh_head.vertices)), int64(np.array(mesh_head.faces)+1), alpha, beta, nargout=1)

def get_vertices_oriented_flagellum(vertices_tail,
                                gamma,
                                delta) :

    """Get the oriented vertices of the flagellum.

    Args:
        vertices_tail: The vertices of the tail of the swimmer.
        gamma: The angle gamma of the flagellum.
        delta: The angle delta of the flagellum.

    Returns:
        The oriented vertices of the flagellum.
    """
    R_gamma = utils.rotation_matrix_y(gamma)
    R_delta = utils.rotation_matrix_z(delta)
    R = R_delta @ R_gamma

    return (R @ vertices_tail.T).T

def get_vertices_translated_flagellum(vertices_tail,
                                  rad_section,
                                  xF1, 
                                  alpha, 
                                  beta) :
    """Get the translated vertices of the flagellum.

    Args:
        vertices_tail: The vertices of the tail of the swimmer.
        rad_section: The radius of the section of the flagellum.
        xF1: The xF1 value for the flagellum.
        alpha: The angle alpha of the flagellum.
        beta: The angle beta of the flagellum.
    Returns:
        The translated vertices of the flagellum.
    """
    
    joint_distance = 2*rad_section 
    print("joint_distance = 2*rad_section")

    return vertices_tail + xF1 + joint_distance * np.array([np.cos(alpha)*np.cos(beta),np.sin(beta)*np.cos(alpha),np.sin(alpha)])

def get_vertices_oriented_translated_flagellum(vertices_tail,
                                            gamma,
                                            delta,
                                            rad_section,
                                            xF1,
                                            alpha,
                                            beta) :

    """
    Get the oriented and translated vertices of the flagellum.

    Args:
        vertices_tail: The vertices of the tail of the swimmer.
        gamma: The angle gamma of the flagellum.
        delta: The angle delta of the flagellum.
        rad_section: The radius of the section of the flagellum.
        xF1: The xF1 value for the flagellum.
        alpha: The angle alpha of the flagellum.
        beta: The angle beta of the flagellum.

    Returns:
        The oriented and translated vertices of the flagellum.
    """
    vertices_oriented = get_vertices_oriented_flagellum(vertices_tail, gamma, delta)
    vertices_oriented_translated = get_vertices_translated_flagellum(vertices_oriented, rad_section, xF1, alpha, beta)
    return vertices_oriented_translated

class FlagRigidSwimmer_Data :
    """Class to store the data of the flagellated rigid swimmer."""

    def __init__(self,
                 head_parameters : dict,
                 flagella_parameters : dict, 
                 angles_parameters : dict) :

        """Initialize the data of the flagellated rigid swimmer."""

        self.mesh_head = head_parameters['mesh']

        self.N_tail = flagella_parameters['N_tail']
        self.L_tail = flagella_parameters['L_tail']
        self.rad_tail = flagella_parameters['rad_tail']
        self.step_length_tail = flagella_parameters['step_length_tail']
        self.rad_section_tail = flagella_parameters['rad_section_tail']
        self.height_tail = flagella_parameters['height_tail']
        self.section_type_tail = flagella_parameters['section_type_tail']
        self.ke_type_tail = flagella_parameters['ke_type_tail']
        self.nb_flag = flagella_parameters['nb_flag']

        self.alpha = angles_parameters['alpha']
        self.gamma = angles_parameters['gamma']
        self.beta = angles_parameters['beta']
        self.delta = angles_parameters['delta']

        self.mesh_tail = None
    
    def add_mesh_tail(self, mesh_tail) :
        """
        Add the mesh of the tail to the swimmer.
        """
        self.mesh_tail = mesh_tail

class FlagRigidSwimmer_Solver :
    """Class to solve the flagellated rigid swimmer problem.
    """

    def __init__(self, flagrigidswimmer_data : FlagRigidSwimmer_Data) :

        self.mesh_head = flagrigidswimmer_data.mesh_head

        self.N_tail = flagrigidswimmer_data.N_tail
        self.L_tail = flagrigidswimmer_data.L_tail
        self.rad_tail = flagrigidswimmer_data.rad_tail
        self.step_length_tail = flagrigidswimmer_data.step_length_tail
        self.rad_section_tail = flagrigidswimmer_data.rad_section_tail
        self.height_tail = flagrigidswimmer_data.height_tail
        self.section_type_tail = flagrigidswimmer_data.section_type_tail
        self.ke_type_tail = flagrigidswimmer_data.ke_type_tail
        self.nb_flag = flagrigidswimmer_data.nb_flag

        self.alpha = flagrigidswimmer_data.alpha
        self.gamma = flagrigidswimmer_data.gamma
        self.beta = flagrigidswimmer_data.beta
        self.delta = flagrigidswimmer_data.delta

    
    def solve_pb(self, eng, eps, volume0, num_pb) : 
        """Solve the optimization problem for the flagellated rigid swimmer.

        Args:
            eng: The matlab engine instance.
            eps: The epsilon values for the constraints ((eps[0] for the constraints on the cost, eps[1] for the constraints on the volume)).
            volume0: The initial volume of the swimmer.
            num_pb: The number of the problem to solve (e.g., 1 or 2).
        
        Returns:
            The cost and the constraints of the optimization problem.
        """

        vtx = double(np.array(self.mesh_head.vertices))
        elt = int64(np.array(self.mesh_head.faces)+1)
        

        cost, constraints = eng.solve_pb(vtx, 
                                    elt, 
                                    self.N_tail, 
                                    self.L_tail, 
                                    self.rad_tail, 
                                    self.step_length_tail, 
                                    self.height_tail,
                                    self.rad_section_tail, 
                                    self.section_type_tail, 
                                    self.ke_type_tail, 
                                    eps[0], 
                                    self.nb_flag, 
                                    self.alpha, 
                                    self.gamma, 
                                    self.beta, 
                                    self.delta, 
                                    num_pb, 
                                    nargout=2, 
                                    background=True).result()
    

        volume = self.mesh_head.volume

        f = -cost 
        c1 = abs(volume) - volume0 - eps[1]
        c2 = volume0 - abs(volume) - eps[1]

        return f, [c1, c2, constraints[0][0], constraints[1][0], constraints[2][0]]
    

    def solve(self, eng, NT) : 
        """Solve the dynamic of the swimmer for one stroke decomposed into NT rotation angles of the flagella.

        Args:
            eng: The matlab engine instance.
            NT: The number of rotation angles of the flagella to compute.
        
        Returns:
            Matrix A of size (6, 2*NT) where A = [U1, U2, ..., UNT; Omega1, Omega2, ..., OmegaNT] where Ui and Omegai are the linear and angular velocities of the swimmer for the i-th rotation angle of the flagella.
        """

        vtx = double(np.array(self.mesh_head.vertices))
        elt = int64(np.array(self.mesh_head.faces)+1)
        

        A = eng.solve(vtx, 
                                    elt, 
                                    self.N_tail, 
                                    self.L_tail, 
                                    self.rad_tail, 
                                    self.step_length_tail, 
                                    self.height_tail,
                                    self.rad_section_tail, 
                                    self.section_type_tail, 
                                    self.ke_type_tail, 
                                    self.nb_flag, 
                                    self.alpha, 
                                    self.gamma, 
                                    self.beta, 
                                    self.delta, 
                                    NT, 
                                    nargout=1, 
                                    background=True).result()
    
        return A


class FlagRigidSwimmer_Visualiser :
    """Class to visualise the flagellated rigid swimmer.
    """

    def __init__(self, flagrigidswimmer_data : FlagRigidSwimmer_Data) :
        self.mesh_head = flagrigidswimmer_data.mesh_head
        self.N_tail = flagrigidswimmer_data.N_tail
        self.L_tail = flagrigidswimmer_data.L_tail
        self.rad_tail = flagrigidswimmer_data.rad_tail
        self.step_length_tail = flagrigidswimmer_data.step_length_tail
        self.rad_section_tail = flagrigidswimmer_data.rad_section_tail
        self.height_tail = flagrigidswimmer_data.height_tail
        self.section_type_tail = flagrigidswimmer_data.section_type_tail
        self.ke_type_tail = flagrigidswimmer_data.ke_type_tail
        self.nb_flag = flagrigidswimmer_data.nb_flag
        self.alpha = flagrigidswimmer_data.alpha
        self.gamma = flagrigidswimmer_data.gamma
        self.beta = flagrigidswimmer_data.beta
        self.delta = flagrigidswimmer_data.delta
        self.mesh_tail = flagrigidswimmer_data.mesh_tail
    
    def plot(self, ax = None) :

        if ax is None:
            fig = plt.figure()
            ax = fig.add_subplot(111, projection='3d')

        ax.plot_trisurf(self.mesh_head.vertices[:, 0],
                        self.mesh_head.vertices[:, 1],
                        triangles=self.mesh_head.faces,
                        Z=self.mesh_head.vertices[:, 2],
                        linewidth = 0.2,
                        antialiased = True,
                        alpha=0.5,
                        color = 'grey',
                        edgecolor='black')
        
        ax.plot_trisurf(self.mesh_tail.vertices[:, 0],
                        self.mesh_tail.vertices[:, 1],
                        triangles=self.mesh_tail.faces,
                        Z=self.mesh_tail.vertices[:, 2],
                        linewidth = 0.2,
                        antialiased = True,
                        alpha=0.5,
                        color = 'grey',
                        edgecolor='black')
        
        if self.nb_flag > 1 :
            R_sym = utils.rotation_matrix_x(np.pi)
            vertices_sym = (R_sym @ self.mesh_tail.vertices.T).T
            ax.plot_trisurf(vertices_sym[:, 0],
                            vertices_sym[:, 1],
                            triangles=self.mesh_tail.faces,
                            Z=vertices_sym[:, 2],
                            linewidth = 0.4,
                            antialiased = True,
                            alpha=0.5,
                            color = 'grey',
                            edgecolor='black')
        

        return ax

    def animate() :
        #TODO
        pass


if __name__ == "__main__" :

    import trimesh
    import matlab.engine
    eng = matlab.engine.start_matlab()
    path_matlab_function = "swimmers/flagrigid"
    eng.addpath(path_matlab_function, nargout=0)

    mesh = trimesh.load_mesh('swimmers/flagrigid/unitsphere3D.stl')

    N = 100
    L = 3.0
    rad = 0.2
    step_length = 1.0
    height = 1.5
    rad_section = 0.067
    section_type = "circular"
    ke_type = "schum"
    EPS = [1e-3, 1e-3]
    nb_flag = 2
    alpha = 0.5*np.pi
    gamma = 0.
    beta = 0.
    delta = 0.
    num_pb = 2

    head_parameters = {'mesh' : mesh}
    flagella_parameters = {'N_tail' : N, 
                          'L_tail' : L, 
                          'rad_tail' : rad, 
                          'step_length_tail' : step_length, 
                          'rad_section_tail' : rad_section, 
                           'height_tail' : height,
                          'section_type_tail' : section_type, 
                          'ke_type_tail' : ke_type,
                          'nb_flag' : nb_flag}
    angles_parameters = {'alpha' : alpha, 'gamma' : gamma, 'beta' : beta, 'delta' : delta}

    filename = "/user/lpalazzo/home/Documents/These/Swimmers/"

    write_msh_flagellum(eng, N, L, rad, step_length, height, rad_section, section_type, ke_type, filename)
    xF1 = get_xF1(eng, mesh, alpha, beta)
    mesh_tail = trimesh.load_mesh(filename+"mesh_tail.ply")
    vertices_oriented_translated_flagellum = get_vertices_oriented_translated_flagellum(mesh_tail.vertices, gamma, delta, rad_section, xF1, alpha, beta)
    mesh_tail.vertices = vertices_oriented_translated_flagellum

    swimmerdata = FlagRigidSwimmer_Data(head_parameters, flagella_parameters, angles_parameters)
    swimmerdata.add_mesh_tail(mesh_tail)

    # visualiser = FlagRigidSwimmer_Visualiser(swimmerdata)
    # ax = visualiser.plot()
    # ax.set_aspect('equal')
    # plt.show()

    swimmer = FlagRigidSwimmer_Solver(swimmerdata)
    swimmer.solve(eng, NT=4.0)
    # f, c = swimmer.solve(eng, EPS, mesh.volume, num_pb)
    # print("f = ", f)
    # print("c = ", c)

    eng.quit()

