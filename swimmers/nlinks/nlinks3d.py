
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
from scipy.interpolate import interp1d
import matlab
from swimmers.nlinks.utils import numpy_to_matlab_double, rotation_x, rotation_y, rotation_z, vector_to_euler_angles




class Nlink3Dmagn_Data :

    def __init__(self, N, r, l, kLT, kLP, kHT, kHP, kHR, k_el, M):
        """Data class for N-link 3D magnetic swimmer parameters.

        Args:
            N (int): Number of links.
            r (float): Radius of the head.
            l (float): Length of each link.
            kLT (float): Tengential tail RFT coefficient.
            kLP (float): Perpendicular tail RFT coefficient.
            kHT (float): Tengential head RFT coefficient.
            kHP (float): Perpendicular head RFT coefficient.
            kHR (float): Rotational head RFT coefficient.
            k_el (float): Elastic constant.
            M (float): Magnetization.
        """

        self.N = N
        self.r = r
        self.l = l
        self.kLT = kLT
        self.kLP = kLP
        self.kHT = kHT
        self.kHP = kHP
        self.kHR = kHR
        self.k_el = k_el
        self.M = M
    


class Nlink3Dmagn_Solver: 

    def __init__(self, nlink3Dmagn_data, t0, tf, z0, timeinterp, ux, uy, uz):
        """Solver class for N-link 3D magnetic swimmer dynamics.

        Args:
            nlink3Dmagn_data (Nlink3Dmagn_Data): Data object containing swimmer parameters.
            t0 (float): Initial time.
            tf (float): Final time.
            z0 (np.array): Initial state vector.
            timeinterp (np.array): Time points for interpolation of magnetic field.
            ux (np.array): Control input in x-direction.
            uy (np.array): Control input in y-direction.
            uz (np.array): Control input in z-direction.
        """

        self.N = nlink3Dmagn_data.N
        self.r = nlink3Dmagn_data.r
        self.l = nlink3Dmagn_data.l
        self.kLT = nlink3Dmagn_data.kLT
        self.kLP = nlink3Dmagn_data.kLP
        self.kHT = nlink3Dmagn_data.kHT
        self.kHP = nlink3Dmagn_data.kHP
        self.kHR = nlink3Dmagn_data.kHR
        self.k_el = nlink3Dmagn_data.k_el
        self.M = nlink3Dmagn_data.M

        self.t0 = t0
        self.tf = tf
        self.z0 = z0
        self.timeinterp = timeinterp
        self.ux = ux
        self.uy = uy
        self.uz = uz

        self.tout = None
        self.zout = None

        self.tout_interp = None
        self.zout_interp = None
        
        self.paramms = [
            ("N", self.N),
            ("r", self.r),
            ("l", self.l),
            ("kLT", self.kLT),
            ("kLP", self.kLP),
            ("kHT", self.kHT),
            ("kHP", self.kHP),
            ("kHR", self.kHR),
            ("k_el", self.k_el),
            ("M", self.M),
            ("t0", t0),
            ("tf", tf),
            ("z0", z0),
            ("timeinterp.shape", timeinterp.shape),
            ("ux.shape", ux.shape),
            ("uy.shape", uy.shape),
            ("uz.shape", uz.shape),
        ]
       

    def print_parameters(self, verbose=True):
        """
        Print the parameters of the Nlink3Dmagn solver.
        """
        if verbose:
            
            # Calculate column widths
            key_width = max(len(k) for k, _ in self.paramms)
            val_width = max(len(str(v)) for _, v in self.paramms)

            title = " Parameters for Nlink3Dmagn "
            print("\n" + title.center(key_width + val_width + 7, "="))
            # Header
            print("_" * (key_width + val_width + 7))
            print(f"| {'Parameter'.ljust(key_width)} | {'Value'.ljust(val_width)} |")
            print("-" * (key_width + val_width + 7))

            # Rows
            for key, value in self.paramms:
                print(f"| {key.ljust(key_width)} | {str(value).ljust(val_width)} |")

            # Footer
            print("-" * (key_width + val_width + 7))

    

    def solve(self, eng, verbose=True):
        """
        Solve the N-link 3D magnetic swimmer dynamics using MATLAB engine.

        Args:
            eng: MATLAB engine instance.
            verbose (bool): Boolean flag to print parameters.

        Returns:
            tuple: Time points and state vectors as NumPy arrays.
        """
        
        self.print_parameters(verbose)
        
        z0 = numpy_to_matlab_double(self.z0)
        ux = numpy_to_matlab_double(self.ux)
        uy = numpy_to_matlab_double(self.uy)
        uz = numpy_to_matlab_double(self.uz)
        timeinterp = numpy_to_matlab_double(self.timeinterp)
        t0 = matlab.double(self.t0)
        tf = matlab.double(self.tf)
        N = matlab.double(self.N)
        tout, zout = eng.solve(N, self.r, self.l, self.kLT, self.kLP, self.kHT, self.kHP, self.kHR, self.k_el, self.M, t0, tf, z0, timeinterp, ux, uy, uz,  nargout=2)
        self.zout = np.asarray(zout)
        self.tout = np.asarray(tout)[:,0]

        return self.tout, self.zout
    

    def interpolate_results(self, tout_interp) :
        """Interpolate the simulation results at specified time points.

        Args:
            tout_interp (np.array): Time points for interpolation.

        Returns:
            tuple: Interpolated time points and state vectors as NumPy arrays.
        """

        self.tout_interp = tout_interp
        self.zout_interp = np.zeros((tout_interp.shape[0], self.zout.shape[1]))
        for i in range(2*self.N + 6):
            self.zout_interp[:, i] = interp1d(self.tout, self.zout[:, i], bounds_error=False, fill_value='extrapolate')(tout_interp)

        return self.tout_interp, self.zout_interp

class Nlink3Dmagn_Visualiser :

    def __init__(self, nlink3Dmagn_data) :
        """Visualiser class for N-link 3D magnetic swimmer.

        Args:
            nlink3Dmagn_data (Nlink3Dmagn_Data): Data object containing swimmer parameters.
        """ 

        self.N = nlink3Dmagn_data.N
        self.r = nlink3Dmagn_data.r
        self.l = nlink3Dmagn_data.l
        self.kLT = nlink3Dmagn_data.kLT
        self.kLP = nlink3Dmagn_data.kLP
        self.kHT = nlink3Dmagn_data.kHT
        self.kHP = nlink3Dmagn_data.kHP
        self.kHR = nlink3Dmagn_data.kHR
        self.k_el = nlink3Dmagn_data.k_el
        self.M = nlink3Dmagn_data.M
    
    def plot(self, state, ax=None, color='royalblue'): 
        """
        Plot the N-link 3D magnetic swimmer in 3D space.
        state is dimension (2*N+6) where :
        - state[0:3] is the position of the head
        - state[3:6] is the orientation of the head (theta_x, theta_y, theta_z)
        - state[6:6+2*N] is the orientation of the link (phi^1_x, phi^1_y, phi^1_z, ..., phi^N_x, phi^N_y, phi^N_z)

        Args:
            state (np.array): State vector of the swimmer. 
            ax (matplotlib.axes._subplots.Axes3DSubplot, optional): Matplotlib 3D axis. Defaults to None.
            color (str, optional): Color of the swimmer. Defaults to 'royalblue'.
        
        Returns:
            matplotlib.axes._subplots.Axes3DSubplot: The 3D axis with the swimmer plotted.
        """
        e1 = np.array([1, 0, 0])

        X_head = state[0:3]
        theta_head = state[3:6]
        phi_links = state[6:2*self.N+6]

        R_head = rotation_x(theta_head[0]) @ rotation_y(theta_head[1]) @ rotation_z(theta_head[2])
        R_links = np.zeros((self.N, 3, 3))
        for i in range(self.N) :
            phi_i = phi_links[i*2:(i+1)*2]
            R_links[i] = rotation_y(phi_i[0]) @ rotation_z(phi_i[1])
        
        X_links = np.zeros((self.N+1, 3))
        X_links[0] = X_head - self.r * R_head @ e1
        for i in range(1, self.N+1): 
            X_links[i] = X_links[i-1] - self.l * R_head @ R_links[i-1] @ e1
        
        if ax is None:
            fig = plt.figure()
            ax = fig.add_subplot(111, projection='3d')

        
        # Create sphere mesh
        phi_head_plot = np.linspace(0, np.pi, 50)
        theta_head_plot = np.linspace(0, 2 * np.pi, 50)
        phi_head_plot, theta_head_plot = np.meshgrid(phi_head_plot, theta_head_plot)
        x_head_plot = self.r * np.sin(phi_head_plot) * np.cos(theta_head_plot) + X_head[0]
        y_head_plot = self.r * np.sin(phi_head_plot) * np.sin(theta_head_plot) + X_head[1]
        z_head_plot = self.r * np.cos(phi_head_plot) + X_head[2]

        # Plotting
        ax.plot_surface(x_head_plot, y_head_plot, z_head_plot, color=color)
        ax.plot(X_links[:, 0], X_links[:, 1], X_links[:, 2], linewidth=3, color=color)#, marker='o', markersize=5)

        return ax

class Nlink3Dmagn_MultipleVisualiser(Nlink3Dmagn_Visualiser):

    def __init__(self, list_nlink3Dmagn_visualiser):

        self.list_nlink3Dmagn_visualiser = list_nlink3Dmagn_visualiser
        self.nb_of_swimmer = len(list_nlink3Dmagn_visualiser)

    def multiple_plot(self, states, ax=None, colors=[]):
        """
        Plot multiple N-link 3D magnetic swimmers in 3D space.
        states is a list of state vectors, each of dimension (2*N+6)

        Args:
            states (List[np.array]): List of state vectors of the swimmers. 
            ax (matplotlib.axes._subplots.Axes3DSubplot, optional): Matplotlib 3D axis. Defaults to None.
            colors (list, optional): List of colors for the swimmers. Defaults to [].   
        
        Returns:
            matplotlib.axes._subplots.Axes3DSubplot: The 3D axis with the swimmers plotted.
        """

        if len(colors) == 0:
            colors = ['royalblue'] * self.nb_of_swimmer
        
        if ax is None:
            fig = plt.figure()
            ax = fig.add_subplot(111, projection='3d')

        for i, state in enumerate(states):
            visualiser_i = self.list_nlink3Dmagn_visualiser[i]
            ax = visualiser_i.plot(state, ax=ax, color=colors[i])

        return ax
   
    def animate(self, states, ax_lims, view_init, path=None, interval = 100, every_frame = 1, colors=[]):
        """ Animate multiple N-link 3D magnetic swimmers in 3D space.
        Args:
            states (List): List of state vectors of the swimmers.
            ax_lims (tuple): Axis limits for the 3D plot.
            view_init (tuple): Elevation and azimuthal angles for the 3D view.
            path (str, optional): Path to save the animation. Defaults to None.
            interval (int, optional): Interval between frames in milliseconds. Defaults to 100.
            every_frame (int, optional): Show every nth frame. Defaults to 1.
            colors (list, optional): List of colors for the swimmers. Defaults to [].

        Returns:
            FuncAnimation: The animation object.
        """
      
        if len(colors) == 0:
            colors = ['royalblue'] * self.nb_of_swimmer

        fig = plt.figure()
        ax = fig.add_subplot(111, projection='3d')
        ax.set_xlim(ax_lims[0])
        ax.set_ylim(ax_lims[1])
        ax.set_zlim(ax_lims[2])
        ax.view_init(elev=view_init[0], azim=view_init[1])

        def update(frame):
            ax.clear()
            ax.set_xlim(ax_lims[0])
            ax.set_ylim(ax_lims[1])
            ax.set_zlim(ax_lims[2])
            ax.view_init(elev=view_init[0], azim=view_init[1])
            list_states_frame = [state[frame, :] for state in states]
            self.multiple_plot(list_states_frame, ax=ax, colors=colors)
            return ax,

        ani = FuncAnimation(fig, update, frames=range(0, np.shape(states[0])[0], every_frame), interval=interval)

        if path is not None:
            ani.save(path+".mp4",  writer='ffmpeg', fps=150)

        return ani  
        

def which_link(s, N, L) :
    """Determine which link corresponds to the index s.

    Args:
        s (float): arc length along the links.
        N (int): number of links.
        L (float): total length the flagellum.

    Returns:
        int: index of the link corresponding to the arc length s.
    """    

    if s==0 :
        return int(0)
    else :
        return int(np.clip(np.ceil(N*s/L)-1,None, N-1))

def Nlink_to_position(state, N, L, r, s) :
    """Compute the position along the N-link swimmer at arc length s.

    Args:
        state (np.array): State vector of the swimmer.
        N (int): Number of links.
        L (float): Total length of the flagellum.
        r (float): Radius of the head.
        s (float): Arc length along the links.  
    
    Returns:
        np.array: Position vector at arc length s.
    """

    l = L / N


    X_head = state[0:3]
    theta_head = state[3:6]
    phi_links = state[6:2*N+6]
    e1 = np.array([1, 0, 0])
    R_head = np.eye(3)#rotation_x(theta_head[0]) @ rotation_y(theta_head[1]) @ rotation_z(theta_head[2])
    R_links = np.zeros((N, 3, 3))
    for i in range(N) :
        phi_i = phi_links[i*2:(i+1)*2]
        R_links[i] = rotation_y(phi_i[0]) @ rotation_z(phi_i[1])
    
    index_link = which_link(s, N, L)

    position = - r * R_head @ e1
    for i in range(index_link) :
        position -= l * R_head @ R_links[i] @ e1
    
    s_tilde = s - index_link * l
    position -= s_tilde * R_head @ R_links[index_link] @ e1

    return position



def Nlink_to_Mlink(L, r, state_Nlink, N, M) :
    """Convert N-link swimmer state to M-link swimmer state.

    Args:
        L (float): Total length of the flagellum.
        r (float): Radius of the head.
        state_Nlink (np.array): State vector of the N-link swimmer.
        N (int): Number of links in the N-link swimmer.
        M (int): Number of links in the M-link swimmer.

    Returns:
        np.array: State vector of the M-link swimmer.
    """

    state_Mlink = np.zeros(2*M+6)
    state_Mlink[0:3] = np.copy(state_Nlink[0:3])
    state_Mlink[3:6] = np.copy(state_Nlink[3:6])

    R_head =  np.eye(3)#(rotation_x(state_Nlink[3]) @ rotation_y(state_Nlink[4]) @ rotation_z(state_Nlink[5]))
    positions_Mlink = np.zeros((M,3))
    for i in range(1,M+1) :
        s = L*i/M
        positions_Mlink[i-1] = Nlink_to_position(state_Nlink, N, L, r, s)
    
    head_junction = - r * R_head @ np.array([1, 0, 0])
    vectors_Mlink = np.zeros((M,3))
    vectors_Mlink[0] = head_junction - positions_Mlink[0]
    for i in range(1, M) :
        vectors_Mlink[i] = positions_Mlink[i-1] - positions_Mlink[i]
    
    for i in range(M):
        phi_y, phi_z = vector_to_euler_angles(vectors_Mlink[i])
        state_Mlink[i*2+6:(i+1)*2+6] = np.array([phi_y, phi_z])
    
    return state_Mlink

