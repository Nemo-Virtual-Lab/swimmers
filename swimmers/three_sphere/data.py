import numpy as np
from swimmers.utils import Box


class ThreeSphere2D:



    def __init__(self, center: np.array, radii: np.array, arms_length: np.array, angles: np.array) -> None:
        """Class that defines the 3-sphere geometry in 2D.

        Args:
            center (np.array): Array with the x and y coordinates of the center of the 3-sphere geometry : [x,y]. Dimensions: (2,)
            radii (np.array): Array with the radii of the 3 spheres : [radleft,radcenter, radright]. Dimensions: (3,)
            arms_length (np.array): Array with the lengths of the arms that connect the spheres : [arm_left, arm_right]. Dimensions: (2,)
            angles (np.array): Array with the angles (number of Pi) of the arms that connect the spheres : [angle_left, angle_right] such as theta_left = angle_left*Pi and theta_right = angle_right*Pi. Dimensions: (2,)
        
        Attributes:
            center (np.array): Array with the x and y coordinates of the center of the 3-sphere geometry : [x,y]. Dimensions: (2,)
            radii (np.array): Array with the radii of the 3 spheres : [radleft,radcenter, radright]. Dimensions: (3,)
            arms_length (np.array): Array with the lengths of the arms that connect the spheres : [arm_left, arm_right]. Dimensions: (2,)
            angles (np.array): Array with the angles (number of Pi) of the arms that connect the spheres : [angle_left, angle_right] such as theta_left = angle_left*Pi and theta_right = angle_right*Pi. Dimensions: (2,)
            path_geo (str): Path to save the .geo file.
            path_json (str): Path to save the .json file.
            path_cfg (str): Path to save the .cfg file.
        """

        self.center = center
        self.radii = radii
        self.arms_length = arms_length
        self.angles = angles

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

        rad_circle_left = self.radii[0]
        rad_circle_center = self.radii[1]
        rad_circle_right = self.radii[2]

        arm_length_left = self.arms_length[0]
        arm_length_right = self.arms_length[1]

        angle_left = self.angles[0]
        angle_right = self.angles[1]

        box_point_x = box.point[0]
        box_point_y = box.point[1]

        box_length_x = box.lengths[0]
        box_length_y = box.lengths[1]

        self.path_geo = path


        geo_file ="""
        h = """ + str(hsize_box) + """;
        lcCircle = """ + str(hsize_swimmer) + """;
        lcDom = h;

        // The construction of the center circle
        RCircle_center = """ + str(rad_circle_center) + """;
        Centerx = """ + str(center_x) + """;
        Centery = """ + str(center_y) + """;

        Point(9) = {Centerx,Centery,0,lcCircle};
        Point(10) = {Centerx+RCircle_center,Centery,0,lcCircle};
        Point(11) = {Centerx-RCircle_center,Centery,0,lcCircle};
        Circle(7) = {10,9,11};
        Circle(8) = {11,9,10};

        // The construction of the left circle
        theta_left = Pi * """ + str(angle_left) + """;
        RCircle_left = """ + str(rad_circle_left) + """;
        ArmLength_left = """ + str(arm_length_left) + """;
        Centerx_left = Centerx + ArmLength_left*Cos(theta_left);
        Centery_left = Centery + ArmLength_left*Sin(theta_left);

        Point(1) = {Centerx_left,Centery_left,0,lcCircle};
        Point(2) = {Centerx_left+RCircle_left,Centery_left,0,lcCircle};
        Point(3) = {Centerx_left,RCircle_left+Centery_left,0,lcCircle};
        Point(4) = {Centerx_left,-RCircle_left+Centery_left,0,lcCircle};
        Circle(1) = {2,1,3};
        Circle(2) = {4,1,2};
        Circle(3) = {3,1,4};

        // The construction of the right circle
        theta_right = Pi * """ + str(angle_right) + """;
        RCircle_right = """ + str(rad_circle_right) + """;
        ArmLength_right = """ + str(arm_length_right) + """;
        Centerx_right = Centerx + ArmLength_right*Cos(theta_right);
        Centery_right = Centery + ArmLength_right*Sin(theta_right);

        Point(5) = {Centerx_right,Centery_right,0,lcCircle};
        Point(6) = {Centerx_right-RCircle_right,Centery_right,0,lcCircle};
        Point(7) = {Centerx_right,RCircle_right+Centery_right,0,lcCircle};
        Point(8) = {Centerx_right,-RCircle_right+Centery_right,0,lcCircle};
        Circle(4) = {7,5,6};
        Circle(5) = {6,5,8};
        Circle(6) = {8,5,7};

        // Defining the arcs and surfaces of the three circles
        Line Loop(7) = {1,2,3};
        Plane Surface(8) = {7};
        Line Loop(9) = {4,5,6};
        Plane Surface(10) = {9};
        Line Loop(11) = {7,8};
        Plane Surface(11) = {11};

        // Constructing the two arms linking the three spheres
        Line(12) = {11, 2};
        Line(13) = {10, 6};

        // Rectangle vertices
        Point(12) = {""" + str(box_point_x) + """,""" + str(box_point_y) + """,0,lcDom};
        Point(13) = {""" + str(box_point_x + box_length_x) + """,""" + str(box_point_y) + """,0,lcDom};
        Point(14) = {""" + str(box_point_x + box_length_x) + """,""" + str(box_point_y + box_length_y) + """,0,lcDom};
        Point(15) = {""" + str(box_point_x) + """,""" + str(box_point_y + box_length_y) + """,0,lcDom};

        // Rectangle lines
        Line(16) = {12, 13};
        Line(17) = {13, 14};
        Line(18) = {14, 15};
        Line(19) = {15, 12};
        Line Loop(20) = {16,17,18,19};

        // Defining the rectangle surface
        Plane Surface(21) = {20, 7, 9, 11};


        Physical Curve("CircleLeft") = {1, 3, 2};
        Physical Curve("CircleCenter") = {7, 8};
        Physical Curve("CircleRight") = {4, 5, 6};
        Physical Curve("BoxWalls") = {16, 17, 18, 19};

        Physical Surface("CirLeft") = {8};
        Physical Surface("CirCenter") = {11};
        Physical Surface("CirRight") = {10};
        Physical Surface("Fluid") = {21};
        """


        with open(path, 'w') as f:
            f.write(geo_file)
            

        return None
    

    def generate_json(self, forceRange: float, epsBody: float, epsWall: float, path: str) -> None:
        """Generates the .json file with the 3-sphere geometry in 2D.

        Args:
            forceRange (float): Range of the force.
            epsBody (float): Epsilon for the body.
            epsWall (float): Epsilon for the wall.
            path (str): Path to save the .json file.

        Returns:
            None
        """

        self.path_json = path

        json_file = """
        {
            "Name": "three_sphere_2D swimmer",
            "ShortName":"three_sphere_2D swimmer",
            "Models":
            {
                "fluid":
                {
                    "materials":"Fluid",
                    "setup":{
                        "equations":"Stokes"
                    }
                },
                "body":{
                    "materials":["CirLeft","CirCenter","CirRight"]
                }
            },
            "Meshes":
            {
                "fluid":
                {
                    "Import":
                    {
                        "filename":\""""+self.path_geo+"""\",
                        "hsize":3
                    },
                    "MeshMotion":
                    {
                        "ComputationalDomain":
                        {
                            "markers":"Fluid",
                            "method":"harmonic_extension"
                        },
                        "Displacement":
                        {
                            "Zero":["BoxWalls"]
                        }
                    }
                }
            },
            "Materials":
            {
                "Fluid":{
                    "rho":"1",
                    "mu":"1"
                },
                "CirLeft":{
                    "rho":1e-1
                },
                "CirCenter":{
                    "rho":1e-1
                },
                "CirRight":{
                    "rho":1e-1
                }
            },
            "Parameters":
            {
                "eps":1e-10
            },
            "BoundaryConditions":
            {
                "fluid":
                {
                    "velocity":
                    {
                        "BoxWalls":
                        {
                            "expr":"{0,0}"
                        }
                    },
                    "body":
                    {
                        "CircleCenter":
                        {
                            "markers":["CircleCenter"],
                            "materials":
                            {
                                "names":["CirCenter"]
                            }
                        },
                        "CircleRight":
                        {
                            "markers":["CircleRight"],
                            "materials":
                            {
                                "names":["CirRight"]
                            },
                            "articulation":
                            {
                                "body":"CircleCenter",
                                "translational-velocity":"u_right:u_right"
                            }
                        },
                        "CircleLeft":
                        {
                            "markers":["CircleLeft"],
                            "materials":
                            {
                                "names":["CirLeft"]
                            },
                            "articulation":
                            {
                                "body":"CircleCenter",
                                "translational-velocity":"u_left:u_left"
                            }
                        }
                    }
                }
            },
            "CollisionForce":  
            {
                "body":
                {
                    "setup": 
                    {
                        "model":"contactAvoidance",
                        "type":"articulatedBody", 
                        "forceParam":
                        {
                            "forceRange":""" + str(forceRange) + """, 
                            "epsBody": """ + str(epsBody) + """,
                            "epsWall": """ + str(epsWall) + """
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
                        "fields":["velocity","pressure","pid","displacement"]
                    },
                    "Measures":
                    {
                        "Quantities":
                        {
                            "names":"all"
                        },
                        "Forces":["CircleCenter","CircleLeft","CircleRight"]
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

        case.dimension=2

        [fluid]
        exporter.use-static-mesh=0
        filename= """ + self.path_json + """
        #solver=Oseen #Oseen,Picard,Newton
        solver=Oseen#Newton
        ksp-monitor=false
        pc-type=lu
        ksp-type=preonly
        #reuse-prec=1
        ksp-maxit-reuse=20
        snes-monitor=true
        snes-maxit=100
        define-pressure-cst=true
        #define-pressure-cst.method=lagrange-multiplier#algebraic
        verbose_solvertimer=0

        #body.articulation.method=p-matrix

        [fluid.alemesh]
        pc-type=lu
        [fluid.bdf]
        order=2

        [ts]
        time-step= """ + str(dt) + """
        time-final= """ + str(Tend) + """
        restart.at-last-save=true
        time-initial=0

        [exporter]
        freq=1
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

        self.forceRange = json_param["forceRange"]
        self.epsBody = json_param["epsBody"]
        self.epsWall = json_param["epsWall"]

        self.dt = cfg_param["dt"]
        self.Tend = cfg_param["Tend"]
        self.directory = cfg_param["directory"]

        self.generate_geomesh(self.box, self.hsize_box, self.hsize_swimmer, path_geo)
        self.generate_json(self.forceRange, self.epsBody, self.epsWall, path_json)
        self.generate_cfg(self.directory, self.dt, self.Tend, path_cfg)

        print("Files generated successfully!")


        return None






