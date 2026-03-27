
import numpy as np
from swimmers.three_sphere.data import ThreeSphere2D
from swimmers.utils import init_environment, Box

from typing import Callable
import matplotlib.pyplot as plt
import feelpp.core as fppc
import feelpp.core.quality as q
from feelpp.toolboxes.fluid import *
import json
import pandas as pd



def get_center_mass(measures, init_center) :
    """Extract the center of mass trajectory of the swimmer's center sphere from simulation measures.

    Args:
        measures (pd.DataFrame): The simulation measures DataFrame.
        init_center (np.array): The initial center position of the swimmer's center sphere.

    Returns:
        np.array: The center of mass trajectory of the swimmer's center sphere.
    """
    center_mass_0 = measures["Quantities_body_CircleCenter.mass_center_0"].values
    center_mass_1 = measures["Quantities_body_CircleCenter.mass_center_1"].values
    center_mass_array = np.column_stack((center_mass_0, center_mass_1))
    center_mass_array = np.concatenate([init_center[None,:],center_mass_array], axis=0)
    return center_mass_array

def get_center_mass_left(measures, init_center) :
    """Extract the center of mass trajectory of the swimmer's left sphere from simulation measures.

    Args:
        measures (pd.DataFrame): The simulation measures DataFrame.
        init_center (np.array): The initial center position of the swimmer's left sphere.

    Returns:
        np.array: The center of mass trajectory of the swimmer's left sphere.
    """
    center_mass_0 = measures["Quantities_body_CircleLeft.mass_center_0"].values
    center_mass_1 = measures["Quantities_body_CircleLeft.mass_center_1"].values
    center_mass_array = np.column_stack((center_mass_0, center_mass_1))
    center_mass_array = np.concatenate([init_center[None,:],center_mass_array], axis=0)
    return center_mass_array

def get_center_mass_right(measures, init_center) :
    """Extract the center of mass trajectory of the swimmer's right sphere from simulation measures.

    Args:
        measures (pd.DataFrame): The simulation measures DataFrame.
        init_center (np.array): The initial center position of the swimmer's right sphere.

    Returns:
        np.array: The center of mass trajectory of the swimmer's right sphere.
    """
    center_mass_0 = measures["Quantities_body_CircleRight.mass_center_0"].values
    center_mass_1 = measures["Quantities_body_CircleRight.mass_center_1"].values
    center_mass_array = np.column_stack((center_mass_0, center_mass_1))
    center_mass_array = np.concatenate([init_center[None,:], center_mass_array], axis=0)
    return center_mass_array

def get_angle(measures, init_angle) :
    """Extract the rotation angle trajectory of the swimmer's center sphere from simulation measures.

    Args:
        measures (pd.DataFrame): The simulation measures DataFrame.
        init_angle (float): The initial rotation angle of the swimmer's center sphere.

    Returns:
        np.array: The rotation angle trajectory of the swimmer's center sphere.
    """
    angle_0 = measures["Quantities_nba_articulation_CircleCenter.rigid_rotation_angles"].values
    angle_array = angle_0
    angle_array = np.concatenate([np.array([init_angle]), angle_array], axis=0)
    return angle_array






def solve_three_sphere2D(three_sphere: ThreeSphere2D, u: Callable, hasContact: bool) :
    """Solve the 2D three-sphere swimmer FSI problem using Feel++.

    Args:
        three_sphere (ThreeSphere2D): The three-sphere swimmer object.
        u (Callable): The control function for the swimmer's motion.
        hasContact (bool): Flag indicating if the swimmer is in contact with the fluid.
        
    Returns:
        pd.DataFrame: The simulation results DataFrame.
    """
    #res
    all_dicts = []

    # Three sphere swimmer planar 2D
    fppc.Environment.setConfigFile(three_sphere.path_cfg)


    f = fluid(dim=2, orderVelocity=2, orderPressure=1)
    f.init()
    #f.printAndSaveInfo()


    # Three sphere swimmer planar 2D
    hfar =  1.0
    hclose = three_sphere.hsize_swimmer
    cst = 0.



    def remesh_toolbox(f, hclose, hfar, parent_mesh, cst):
        # 3 spheres
        required_facets=["CircleLeft","CircleCenter","CircleRight"]
        required_elts=["CirLeft","CirCenter","CirRight"]
    

        n_required_elts_before=fppc.nelements(fppc.markedelements(f.mesh(),required_elts))
        n_required_facets_before=fppc.nelements(fppc.markedfaces(f.mesh(),required_facets))
        #print(" . [before remesh]   n required elts: {}".format(n_required_elts_before))
        #print(" . [before remesh] n required facets: {}".format(n_required_facets_before))

        new_mesh, cpt = fppc.remesh(
            mesh=f.mesh(), metric="gradedls({},{})".format(hclose, hfar), required_elts=required_elts, required_facets=required_facets, params='{"remesh":{ "verbose":-1}}')
        
        #print(" . [after remesh]  n remeshes: {}".format(cpt))
        n_required_elts_after=fppc.nelements(fppc.markedelements(new_mesh,required_elts))
        n_required_facets_after=fppc.nelements(fppc.markedfaces(new_mesh,required_facets))
        #print(" . [after remesh]  n required elts: {}".format(n_required_elts_after))
        #print(" . [after remesh] n required facets: {}".format(n_required_facets_after))
        f.applyRemesh(f.mesh(),new_mesh)


    parent_mesh=f.mesh()
    remesh_toolbox(f, hclose, hfar, None , cst)
    #f.exportResults()


    nbr_remesh = 0
    time_remesh = []

    # Reset execution time parameters
    f.reset_executionTime()

    #Add collision force 
    if hasContact:
        f.addContactForceModel()
        f.addContactForceResModel()

    f.startTimeStep()

    while not f.timeStepBase().isFinished():

        
        min_etaq = q.etaQ(f.mesh()).min()
        
        if min_etaq < 0.6:
            remesh_toolbox(f, hclose, hfar, None, cst)
            
            if hasContact:
                f.addContactForceModel()
                f.addContactForceResModel()


            nbr_remesh += 1
            time_remesh.append(f.time())
        
    
        if fppc.Environment.isMasterRank():
            #print("============================================================\n")
            #print("time simulation: {}s iteration : {}\n".format(f.time(), f.timeStepBase().iteration()))
            #print("  -- mesh quality: {}s\n".format(min_etaq))
            #print("============================================================\n")
            pass

        
        u_value = u(f.time()) #u(t)=[u1(t), u2(t)]
        u_left = u_value[0] #u_left = u1(t)
        u_right = u_value[1] #u_right = u2(t)
        
        f.addParameterInModelProperties("u_left", u_left) #modifier paramètres dans le json
        f.addParameterInModelProperties("u_right",u_right) #modifier paramètres dans le json
        f.updateParameterValues() #mise à jour des paramètres
        f.solve()
        
        f.exportResults()
        if not f.postProcessMeasures().empty():
            measure = f.postProcessMeasures().values()
            all_dicts.append(measure)
            #print(f"measure = {measure}")
        f.updateTimeStep() #t = t+dt
    
    return pd.DataFrame(all_dicts)




if __name__ == "__main__":

    import mpi4py
    mpi4py.rc.thread_level="single"

    # ============ Geometrical parameters
    center = np.array([4,3])
    radii = np.array([0.125,0.125,0.125])
    arms_length = np.array([10*0.125,10*0.125])
    angles = np.array([1., 0.])

    hsize_box = 0.5
    hsize_swimmer = 0.06
    box = Box(np.array([0,0]), np.array([10,6]))

    path_geo = "three_sphere_2D.geo"
    path_json = "three_sphere_2D.json"
    path_cfg = "three_sphere_2D.cfg"
    dt = 0.05
    Tend = 4
    directory = "."

    geo_param = {"box": box, "hsize_box": hsize_box, "hsize_swimmer": hsize_swimmer}
    json_param = {"forceRange": hsize_swimmer, "epsBody": 0.00000005, "epsWall": 0.00000005}
    cfg_param = {"dt": dt, "Tend": Tend, "directory": directory}

    def u(t) :

        def pulse(t, t0, t1):
            # pulse(t, t0, t1) = 1 if t0+ <= t <= t1 else 0 
            return np.heaviside(t-t0, 1) - np.heaviside(t-t1, 1)
        
        def periodic_pulse(t, t0, t1, T):
            t_mod = (t - t0) % T + t0
            return pulse(t_mod, t0, t1)
        
        u_left = 0.4 * ( periodic_pulse(t, 0, 1, 4) - periodic_pulse(t, 2, 3, 4) )
        u_right = 0.4 * ( periodic_pulse(t, 1, 2, 4) - periodic_pulse(t, 3, 4, 4) )
        return np.array([u_left, u_right])
  

    three_sphere = ThreeSphere2D(center, radii, arms_length, angles)
    three_sphere.generate_files(path_geo, path_json, path_cfg, geo_param, json_param, cfg_param)
    
    e = init_environment(".")
    m = solve_three_sphere2D(three_sphere, u, False)
    print(m)
    print(get_center_mass(m, three_sphere.center))
  

    

