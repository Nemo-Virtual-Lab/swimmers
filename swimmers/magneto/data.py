import numpy as np
from swimmers.utils import Box


class Magneto2D:

    def __init__(self, center: np.array, angle: float, head_length_height: np.array, tail_length: float) -> None:
        """Class that defines a 2D magneto-swimmer geometry and generates .geo, .json and .cfg files for Feel++ simulations.

        Args:
            center (np.array): Array with the x and y coordinates of the center of the swimmer : [x,y]. Dimensions: (2,)
            angle (float): Angle of rotation of the swimmer in radians.
            head_length_height (np.array): Array with the length and height of the swimmer's head : [length,height]. Dimensions: (2,)
            tail_length (float): Length of the swimmer's tail.
        """

        self.center = center
        self.angle = angle
        self.head_length_head = head_length_height
        self.tail_length = tail_length
        

    def generate_geomesh(self, box : Box, hsize_box: float, hsize_swimmer: float, path: str) -> None:
        """Generates the .geo file with the 3-sphere geometry in 2D.

        Args:
            box (Box): Box object that defines the geometry of the domain.
            hsize_box (float): Size of the mesh elements in the box.
            hsize_swimmer (float): Size of the mesh elements in the 3-sphere geometry.
            path (str): Path to save the .geo file.

        Returns:
            None
        """

        center_x = self.center[0]
        center_y = self.center[1]

        height_head = self.head_length_head[1]
        length_head = self.head_length_head[0]

        length_tail = self.tail_length

        height_box = box.lengths[1]
        length_box = box.lengths[0]

        theta = self.angle


        self.path_geo = path


        geo_file ="""
        SetFactory("OpenCASCADE");

        h = """ + str(hsize_swimmer) + """;
        h1 = """ + str(hsize_box) + """;

        height = """ + str(height_head) + """;
        length_head = """ + str(length_head) + """;
        length_tail = """ + str(length_tail) + """;

        height_box = """ + str(height_box) + """;
        lenght_box = """ + str(length_box) + """;

        theta = """ + str(theta) + """;

        centerx = """ + str(center_x) + """;
        centery = """ + str(center_y) + """;


        Cs = Cos(theta);
        Sn = Sin(theta);

        P_bottom_left_head_x = -length_head/2 + centerx;
        P_bottom_left_head_y = -height/2 + centery;

        P_bottom_right_head_x = length_head/2 + centerx;
        P_bottom_right_head_y = -height/2 + centery;

        P_top_left_head_x = -length_head/2 + centerx;
        P_top_left_head_y = height/2 + centery;

        P_top_right_head_x = length_head/2 + centerx;
        P_top_right_head_y =  height/2 + centery;



        Point(1) = {Cs * P_bottom_left_head_x + Sn * P_bottom_left_head_y, -Sn * P_bottom_left_head_x + Cs * P_bottom_left_head_y, 0, h}; //bottom left
        Point(2) = {Cs * P_top_left_head_x + Sn * P_top_left_head_y, -Sn * P_top_left_head_x + Cs * P_top_left_head_y, 0, h}; // top left
        Point(3) = {Cs * P_top_right_head_x + Sn * P_top_right_head_y, -Sn * P_top_right_head_x + Cs * P_top_right_head_y, 0, h}; // top right
        Point(4) = {Cs * P_bottom_right_head_x + Sn * P_bottom_right_head_y, -Sn * P_bottom_right_head_x + Cs * P_bottom_right_head_y, 0, h}; // bottom right


        P_bottom_tail_x = -length_head/2 - length_tail + centerx;
        P_bottom_tail_y = -length_head/5 + centery;

        P_top_tail_x = -length_head/2 - length_tail + centerx;
        P_top_tail_y = length_head/5 + centery;

        P_middle_tail_x = -length_head/2 - length_tail + length_head/10 + centerx;
        P_middle_tail_y = centery;

        Point(5) = {Cs * P_bottom_tail_x + Sn * P_bottom_tail_y , -Sn * P_bottom_tail_x + Cs * P_bottom_tail_y , 0, h};  // Point de départ de la queue
        Point(6) = {Cs * P_top_tail_x + Sn * P_top_tail_y , -Sn * P_top_tail_x + Cs * P_top_tail_y , 0, h};  // Point de fin de la queue
        Point(7) = {Cs * P_middle_tail_x + Sn * P_middle_tail_y , -Sn * P_middle_tail_x + Cs * P_middle_tail_y , 0, h};  // Point pour début de l'arrondi


        Line(1) = {2, 1};
        Line(2) = {1, 4};
        Line(3) = {4, 3};
        Line(4) = {3, 2};

        Line(28) = {1, 5};
        Line(29) = {2, 6};

        Circle(30) = {6, 7, 5};  



        Point(11) = {-lenght_box/2, -height_box/2, 0, h1};
        Point(12) = {lenght_box/2., -height_box/2, 0, h1};
        Point(13) = {-lenght_box/2, height_box/2, 0, h1};
        Point(14) = {lenght_box/2., height_box/2, 0, h1};

        Line(5) = {13, 11};
        Line(6) = {11, 12};
        Line(7) = {12, 14};
        Line(8) = {14, 13};

        Curve Loop(1) = {29, 30, -28, -1};
        Plane Surface(1) = {1};
        Curve Loop(2) = {4, 1, 2, 3};
        Plane Surface(2) = {2};
        Curve Loop(3) = {8, 5, 6, 7};
        Curve Loop(4) = {29, 30, -28, 2, 3, 4};
        Plane Surface(3) = {3, 4};

        Physical Surface("Solid", 31) = {1};
        Physical Surface("Head", 32) = {2};
        Physical Surface("Fluid", 33) = {3};
        Physical Curve("fluid-inlet", 34) = {5};
        Physical Curve("fluid-wall", 35) = {8, 6};
        Physical Curve("fluid-outlet", 36) = {7};
        Physical Curve("magneto", 37) = {1};
        Physical Curve("fsi-wall", 38) = {29, 30, 28, 2, 3, 4};
        """


        with open(path, 'w') as f:
            f.write(geo_file)
            

        return None
    

    def generate_json(self, fluid_param : dict, tail_param : dict, head_param : dict, hfar: float, hclose: float, path: str) -> None:
        """Generates the .json file with the 3-sphere geometry in 2D.

        Args:
            fluid_param (dict): Dictionary for the fluid : {"rho" : ..., "mu" : ... }
            tail_param (dict): Dictionary for the tail : {"E" : ..., "nu" : ..., "rho" : ...}
            head_param (dict): Dictionary for the head : {"E" : ..., "nu" : ..., "rho" : ..., "M" : ...}
            path (str): Path to save the .json file.

        Returns:
            None
        """

        self.path_json = path
        self.fluid_param = fluid_param
        self.head_param = head_param

        rho_fluid = fluid_param["rho"]
        mu_fluid = fluid_param["mu"]
        E_tail = tail_param["E"]
        nu_tail = tail_param["nu"]
        rho_tail = tail_param["rho"]
        E_head = head_param["E"]
        nu_head = head_param["nu"]
        rho_head = head_param["rho"]
        M_head = head_param["M"]

        json_file = """
                {
            "Name": "Magneto-swimmer",
            "ShortName":"magneto2",
            "Models":
            {
                "fluid":{
                    "materials":["Fluid"],
                    "equations":"Stokes",
                    "equationsAA":"Navier-Stokes"
                },
                "solid":
                {
                    "materials":["Solid","Head"],
                    "setup":{
                        "equations":"Hyper-Elasticity",
                        "material-model":"StVenantKirchhoff"
                    }
                },
                "multibody":{
                    "name":"toto",
                    "setup":{
                        "bodies":{
                            "name":"fsi-wall",
                            "materials":["Solid","Head"],
                            "mass-center":{
                                "evaluate-on-materials":["Head"]
                            }
                        },
                        "connectors": {""" + """},
                        "constraints":{""" + """}
                    }
                },
                "fsi":
                {
                    "materials":["Fluid","Solid","Head"],
                    "setup":{
                        "interfaces":[{
                            "markers":"fsi-wall",
                            "type":"body"
                        }]
                    }
                }
            },
            "Parameters":
            {
                "hfar": """ + str(hfar) + """,
                "hclose_wall1":"hfar:hfar",
                "hclose_wall2": """ + str(hclose) + """,

                "d2r_bounded_fsi_wall1":"min(meshes_fsi_distanceToRange_wall1_normalized_min_max,0.2)/0.2:meshes_fsi_distanceToRange_wall1_normalized_min_max",
                "mymetric_fsi_wall1":"hclose_wall1+(hfar-hclose_wall1)*d2r_bounded_fsi_wall1:hclose_wall1:hfar:d2r_bounded_fsi_wall1",
                "d2r_bounded_fsi_wall2":"min(meshes_fsi_distanceToRange_wall2_normalized_min_max,0.2)/0.2:meshes_fsi_distanceToRange_wall2_normalized_min_max",
                "mymetric_fsi_wall2":"hclose_wall2+(hfar-hclose_wall2)*d2r_bounded_fsi_wall2:hclose_wall2:hfar:d2r_bounded_fsi_wall2",

                "d2r_bounded_wall1":"min(meshes_fluid_distanceToRange_wall1_normalized_min_max,0.2)/0.2:meshes_fluid_distanceToRange_wall1_normalized_min_max",
                "mymetric_wall1":"hclose_wall1+(hfar-hclose_wall1)*d2r_bounded_wall1:hclose_wall1:hfar:d2r_bounded_wall1",
                "d2r_bounded_wall2":"min(meshes_fluid_distanceToRange_wall2_normalized_min_max,0.2)/0.2:meshes_fluid_distanceToRange_wall2_normalized_min_max",
                "mymetric_wall2":"hclose_wall2+(hfar-hclose_wall2)*d2r_bounded_wall2:hclose_wall2:hfar:d2r_bounded_wall2",

                "ux_t":"uxt:uxt",
                "uy_t":"uyt:uyt",
                "uz_t":"uzt:uzt"

            },
            "Meshes":
            {
                "fsi":{
                    "Import":
                    {
                        "filename":\""""+self.path_geo+"""\"
                    },
                    "Partitioning":{
                        "splitting":[ "Fluid", ["Solid","Head"] ],
                        "constraints":{
                            "no_interprocess_faces":["fsi-wall"]
                        }
                    },
                    "DistanceToRange":
                    {
                        "wall1":
                        {
                            "markers":["fluid-inlet","fluid-wall","fluid-outlet"],
                            "normalization":"min_max"
                        },
                        "wall2":
                        {
                            "markers":["fsi-wall"],
                            "normalization":"min_max"
                        }
                    },
                    "MeshAdaptation":[
                        {
                            "metric":"min(mymetric_fsi_wall1,mymetric_fsi_wall2):mymetric_fsi_wall1:mymetric_fsi_wall2",
                            "events":{
                                "after_import":{""" + """}
                            }
                        }
                    ]
                },
                "fluid":
                {
                    "MeshMotion":
                    {
                        "ComputationalDomain":
                        {
                            "markers":"Fluid"
                        },
                        "Displacement":
                        {
                            "Zero":["fluid-inlet","fluid-wall","fluid-outlet"]
                        }
                    },
                    "DistanceToRange":
                    {
                        "wall1":
                        {
                            "markers":["fluid-inlet","fluid-wall","fluid-outlet"],
                            "normalization":"min_max"
                        },
                        "wall2":
                        {
                            "markers":["fsi-wall"],
                            "normalization":"min_max"
                        }
                    },
                    "MeshAdaptation":[

                        {
                            "metric":"min(mymetric_wall1,mymetric_wall2):mymetric_wall1:mymetric_wall2",
                            "required_markers":["fsi-wall","Solid","Head"],
                            "events":{
                                "each_time_step":{
                                    "frequency":20
                                }
                            }
                        }],
                    "Partitioning":{
                        "splitting":[ "Fluid", ["Solid","Head"] ],
                        "constraints":{
                            "no_interprocess_faces":["fsi-wall"]
                        }
                    
                    }
                }
            },
            "Materials":
            {
                "Fluid": {
                    "rho": """ + str(rho_fluid) + """,
                    "mu": """ + str(mu_fluid) + """
                },
                "Solid": {
                    "E": """ + str(E_tail) + """,
                    "nu": """ + str(nu_tail) + """,
                    "rho": """ + str(rho_tail) + """
                },
                "Head": {
                    "E": """ + str(E_head) + """,
                    "nu":""" + str(nu_head) + """,
                    "rho":""" + str(rho_head) + """
                }

            },
            "BoundaryConditions":
            {
                "fluid":
                {
                    "velocity":
                    {
                        "fluid-wall":
                        {
                            "markers":["fluid-wall"],
                            "expr":"{0,0}"
                        }
                    },
                    "outlet_free":
                    {
                        "fluid-outlet":
                        {
                            "expr":"0"
                        },
                        "fluid-inlet":
                        {
                        "expr":"0"
                        }
                    },
                    "body":
                    {
                        "fsi-wall":
                        {
                            "markers":"fsi-wall",
                            "materials":
                            {
                                "names":["Solid","Head"]
                            }
                        }
                    }
                }
            },
            "MagnetoTorque":
            {
                "body":
                {
                    "setup":
                    {
                        "torqueParam":
                        {
                            "mx": """ + str(M_head[0]) + """,
                            "my": """ + str(M_head[1]) + """,
                            "mz": """ + str(M_head[2]) + """
                        }
                    }
                }
            },
            "PostProcess":
            {
                "fluid":
                {
                    "Exports":
                    {
                        "fields":["velocity","pressure","mesh-displacement","pid"]
                    },
                    "Measures":
                    {
                        "Quantities":
                        {
                            "names":"all"
                        }
                    }
                },
                "solid":
                {
                    "Exports":
                    {
                        "fields":["displacement","velocity","pid"]
                    },
                    "Measures":
                    {
                        "Points":
                        {
                            "pointA":
                            {
                                "coord":"{""" + str(self.center[0]) +""",""" + str(self.center[1]) + """}",
                                "fields":["displacement","velocity"]
                            }
                        }
                    }
                }
            }
        }
        """

        with open(path, 'w') as f:
            f.write(json_file)
    

    def generate_cfg(self, directory: str, dt: float, Tend: float, path: str) -> None:

        """Generates the .cfg file with the 3-sphere geometry in 2D.

        Args:
            directory (str): Directory where the results will be saved.
            dt (float): Time step.
            Tend (float): Final time.
            path (str): Path to save the .cfg file.
        
        Returns:
            None
        """
        
        self.path_cfg = path    

        cfg_file = """
        directory=""" + directory + """
        [fsi]
        # files
        filename= """ + self.path_json + """
       
        mesh-save.tag=M2
        # fluid and solid markers
        fluid-mesh.markers=Fluid
        fluid-mesh.markers=Solid
        fluid-mesh.markers=Head
        solid-mesh.markers=Solid
        solid-mesh.markers=Head
        conforming-interface=true
        fixpoint.tol=1e-6#1e-4#1e-6#1e-8
        #fixpoint.initialtheta=0.7#0.1#0.98#0.1#99#0.05
        fixpoint.min_theta=1e-12#0.005#1e-8#1e-4
        fixpoint.maxit=50##10#1#10#2
        coupling-bc=dirichlet-neumann
        #coupling-type=Semi-Implicit #Implicit #Semi-Implicit
        evaluate-fluid-normal-stress-on-reference-mesh=0
        #solve-rigid=true # Resolution of the fluid rigid problem
        #solve-elastic=true # Resolution of the solid problem
        solve-dirichlet=1#true #false#true # Add dirichlet boundary condition onto the head in solid problem
        [fluid]
        exporter.use-static-mesh=0
        #solver=Oseen #Oseen,Picard,Newton
        #solver=Newton
        #snes-monitor=1
        #ksp-converged-reason=true
        #snes-converged-reason=true
        pc-type=lu#gasm #gasm#lu #asm#fieldsplit #ilu
        ksp-type=preonly
        #time-stepping=Theta
        snes-monitor=1
        [fluid.alemesh]
        type=harmonic
        pc-type=gamg#lu
        #ksp-type=preonly
        #pc-type=gamg
        #ksp-maxit=30
        #reuse-prec=true
        ksp-converged-reason=true
        error-if-solver-not-converged=0
        [fluid.alemesh.ho]
        pc-type=gamg
        ksp-maxit=30
        reuse-prec=true
        [fluid.bdf]
        order=2
        [fluid.alemesh.bdf]
        order=2
        [solid]
        #use-null-space=1
        #error-if-solver-not-converged=0
        #solver.nonlinear.apply-dof-elimination-on-initial-guess=0
        #on.type=elimination_symmetric_keep_diagonal #elimination_symmetric
        material_law=StVenantKirchhoff # StVenantKirchhoff, NeoHookean
        snes-monitor=1
        pc-type=lu#gamg
        ksp-type=preonly
        snes-maxit=100
        time-stepping=BDF#Theta #BDF#Theta
        bdf.order=2
        [ts]
        #restart=true
        time-step= """ + str(dt) + """
        time-final= """ + str(Tend) + """
        #restart.at-last-save=true
        save.freq=100
        [exporter]
        #freq=5
        #export=0
        [fluid]
        #verbose=1
        verbose_solvertimer=1
        [solid]
        #verbose=1
        verbose_solvertimer=1
        [fsi]
        #verbose=1
        verbose_solvertimer=1
        """

        with open(path, 'w') as f:
            f.write(cfg_file)
    

    def generate_files(self, path_geo: str, path_json: str, path_cfg: str, geo_param: dict, json_param: dict, cfg_param: dict) -> None:
        """Generates the .geo, .json and .cfg files with the 3-sphere geometry in 2D.

        Args:
            path_geo (str): Path to save the .geo file.
            path_json (str): Path to save the .json file.
            path_cfg (str): Path to save the .cfg file.
            geo_param (dict): Dictionary with the parameters to generate the .geo file.
            json_param (dict): Dictionary with the parameters to generate the .json file.
            cfg_param (dict): Dictionary with the parameters to generate the .cfg file.
        
        Returns:
            None
        """

        self.box = geo_param["box"]
        self.hsize_box = geo_param["hsize_box"]
        self.hsize_swimmer = geo_param["hsize_swimmer"]

        self.fluid_param = json_param["fluid"]
        self.tail_param = json_param["tail"]
        self.head_param = json_param["head"]
        self.remesh_param = json_param["remesh"]

        self.dt = cfg_param["dt"]
        self.Tend = cfg_param["Tend"]
        self.directory = cfg_param["directory"]

        self.generate_geomesh(self.box, self.hsize_box, self.hsize_swimmer, path_geo)
        self.generate_json(self.fluid_param, self.tail_param, self.head_param, self.remesh_param[0], self.remesh_param[1],path_json)
        self.generate_cfg(self.directory, self.dt, self.Tend, path_cfg)

        print("Files generated successfully!")


        return None




