import sys
import numpy as np
from swimmers.utils import init_environment, Box
from swimmers.magneto.data import Magneto2D


from typing import Callable
#import matplotlib.pyplot as plt
import feelpp.core as fppc
from feelpp.toolboxes.fsi import *
import json
import pandas as pd




def get_center_head(measures: pd.DataFrame, init_center: np.array) -> np.array:
    """Extract the center of mass trajectory of the swimmer's head from simulation measures.

    Args:
        measures (pd.DataFrame): The simulation measures DataFrame.
        init_center (np.array): The initial center position of the swimmer's head.

    Returns:
        np.array: The center of mass trajectory of the swimmer's head.
    """

    center_mass_0 = measures["Quantities_body_fsi-wall.mass_center_0"].values
    center_mass_1 = measures["Quantities_body_fsi-wall.mass_center_1"].values
    center_mass_array = np.column_stack((center_mass_0, center_mass_1))
    center_mass_array = np.concatenate([init_center[None,:],center_mass_array], axis=0)
    return center_mass_array

    

def get_measures(csv_path: str) -> pd.DataFrame:
    """Read simulation measures from a CSV file.

    Args:
        csv_path (str): The path to the CSV file.

    Returns:
        pd.DataFrame: The simulation measures DataFrame.
    """
    fd = pd.read_csv(csv_path)
    return fd


def solve_magneto2D(magneto: Magneto2D, u: Callable, hasContact: bool, env: str) -> pd.DataFrame:
    """Solve the 2D magneto-swimmer FSI problem using Feel++.

    Args:
        magneto (Magneto2D): The magneto-swimmer object.
        u (Callable): The control function for the swimmer's motion.
        hasContact (bool): Flag indicating if the swimmer is in contact with the fluid.
        env (str): The environment configuration.

    Returns:
        pd.DataFrame: The simulation results DataFrame.
    """

    #res
    all_dicts = []

    # magneto swimmer planar 2D
    fppc.Environment.setConfigFile(magneto.path_cfg)


    fsi_tb = fsi(dim=2, orderU=2, orderP=1, orderGeo=1)
    fsi_tb.init()
    #fsi_tb.printAndSaveInfo()


    #Add Torque FSI
    fsi_tb.addMagnetoTorqueModelFSI()
    fsi_tb.addMagnetoTroqueResModelFSI()

    fsi_tb.startTimeStep()

    while not fsi_tb.timeStepBase().isFinished():
 
        if fppc.Environment.isMasterRank():
            # print("============================================================\n")
            # print("time simulation: {}s iteration : {}\n".format(fsi_tb.time(), fsi_tb.timeStepBase().iteration()))
            # print("============================================================\n")
            pass
        

        #Update control at time t
        uxt = u(fsi_tb.time())[0]
        uyt = u(fsi_tb.time())[1]
        uzt = u(fsi_tb.time())[2]
        fsi_tb.addParameterInModelProperties("uxt", uxt)
        fsi_tb.addParameterInModelProperties("uyt", uyt)
        fsi_tb.addParameterInModelProperties("uzt", uzt)
        fsi_tb.updateParameterValues()


        #Solve FSI
        fsi_tb.solve()

        #Export results
        fsi_tb.exportResults()
        if not fsi_tb.postProcessMeasures().empty():
            measure = fsi_tb.postProcessMeasures().values()
            all_dicts.append(measure)
         
        #Update time : t <- t+dt
        fsi_tb.updateTimeStep()
   
    
    #return pd.DataFrame(all_dicts)
    print("/!\ Modify to take into account parallelism /!\\")
    return get_measures(env + "/np_1/fluid.measures/values.csv")





if __name__ == "__main__":

    import mpi4py
    mpi4py.rc.thread_level="single"

    hsize_box = 0.001
    hsize_swimmer = 0.0005

    hfar = hsize_swimmer * 5
    hclose = hfar /10.
 

    center = np.array([0,0])
    angle = 0
    head_length_height = np.array([0.0005, 0.0015])
    tail_length = 0.0075

    rho_fluid = 1260
    mu_fluid = 1.52
    E_tail = 8e4
    nu_tail = 0.4
    rho_tail = 1200
    E_head  = 4.1e10
    nu_head = 0.281
    rho_head = 7000
    M_head = np.array([200000, 200000, 1])

  
    box = Box(np.array([0,0]), np.array([0.04,0.02]))

    remesh_param = np.array([hfar, hclose])

    dic_fluid = {"rho" : rho_fluid, "mu" : mu_fluid}
    dic_tail = {"E" : E_tail, "nu" : nu_tail, "rho" : rho_tail}
    dic_head = {"E" : E_head, "nu" : nu_head, "rho" : rho_head, "M" : M_head}

    path_geo = "magneto2D.geo"
    path_json = "magneto2D.json"
    path_cfg = "magneto2D.cfg"

    dt = 0.01
    Tend = 0.03
    directory = "."


    geo_param = {"box": box, "hsize_box": hsize_box, "hsize_swimmer": hsize_swimmer}
    json_param = {"fluid": dic_fluid, "tail" : dic_tail, "head" : dic_head, "remesh" : remesh_param}
    cfg_param = {"dt": dt, "Tend": Tend, "directory": directory}


    #============== Control parameters =======================#
    freq = 0.9
    ux = lambda t : 0.005
    uy = lambda t : 0.005 * np.sin(2*np.pi*freq*t)
    uz = lambda t : 0
    u = lambda t : np.array([ux(t), uy(t), uz(t)])
    #=========================================================#
    

    magneto = Magneto2D(center, angle, head_length_height, tail_length)
    magneto.generate_files(path_geo, path_json, path_cfg, geo_param, json_param, cfg_param)
    
    e = init_environment(".")
    m = solve_magneto2D(magneto, u, False, ".")
    print(m)
    print(get_center_head(m, center))

