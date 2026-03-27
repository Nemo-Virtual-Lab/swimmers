# Swimmers

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)

This repository contains implementations of different swimmer models:

* **N-Link** model in MATLAB and Python
* **Three-sphere** model in C++ and Python using [Feel++](https://docs.feelpp.org/home/index.html)
* **Magneto** model in C++ and Python using [Feel++](https://docs.feelpp.org/home/index.html)

**Documentation** is available [here](https://luplz.github.io/Swimmers/).

---

## 🚀 Installation

1. **Clone the repository**

```bash
git clone https://github.com/luplz/Swimmers.git
cd Swimmers
```

2. **Install Python dependencies**

```bash
pip install -r requirements.txt
```

* For the N-Link model, install the MATLAB engine:

```bash
pip install matlabengine
```

* For the Three-sphere and Magneto models, install Feel++:

[Feel++ Installation Guide](https://docs.feelpp.org/user/latest/install/index.html)

3. **Add the repository root to your PYTHONPATH**

```bash
export PYTHONPATH="$PWD:$PYTHONPATH"
```

---

## 🛠️ Usage Examples

### N-Link

```python
import matlab.engine
eng = matlab.engine.start_matlab()
path_matlab_function = "swimmers/nlinks"
eng.addpath(path_matlab_function, nargout=0)

import numpy as np
import matplotlib.pyplot as plt
from swimmers.nlinks.nlinks3d import (
    Nlink3Dmagn_Data,
    Nlink3Dmagn_Solver,
    Nlink3Dmagn_Visualiser,
    Nlink3Dmagn_MultipleVisualiser
)
```

Set up swimmer parameters:

```python
#Number of links
N1, N2 = 5, 5

# RFT coefficients:
# L = links, H = head
# T = tangential, P = perpendicular
# R = rotational, el = elastic
kLT, kLP = 0.344, 0.813
kHT, kHP, kHR = 1.143, 4.364, 0.590
k_el = 0.868

#length of each links
l1, l2 = 7e-3/N1, 7e-3/N2

#head radius
r = 0.3e-3

#Magnetization
M = 1.6851e12

#Initial state
z01, z02 = np.zeros(2*N1+6), np.zeros(2*N2+6)

#Final time and discretization
tf = 10
timeinterp = np.linspace(0, tf, 1000, dtype=np.float64)

#Control
ux_1, uy_1, uz_1 = 0*timeinterp, 1e-2*np.cos(2*np.pi*1.5*timeinterp), 0*timeinterp
ux_2, uy_2, uz_2 = 0*timeinterp+1e-2, uy_1, uz_1
```

Create data, solver, and visualiser objects:

```python
#Data
nlink_data_1 = Nlink3Dmagn_Data(N1, r, l1, kLT, kLP, kHT, kHP, kHR, k_el, M)
nlink_data_2 = Nlink3Dmagn_Data(N2, r, l2, kLT, kLP, kHT, kHP, kHR, k_el, M)

#Solver
solver_1 = Nlink3Dmagn_Solver(nlink_data_1, 0, tf, z01, timeinterp, ux_1, uy_1, uz_1)
solver_2 = Nlink3Dmagn_Solver(nlink_data_2, 0, tf, z02, timeinterp, ux_2, uy_2, uz_2)

#Solve
tout_1, zout_1_old = solver_1.solve(eng)
tout_2, zout_2_old = solver_2.solve(eng)

#Interpolation with discretization
tout_1_interp, zout_1 = solver_1.interpolate_results(timeinterp)
tout_2_interp, zout_2 = solver_2.interpolate_results(timeinterp)

#Visualiser for each swimmer
vis_1 = Nlink3Dmagn_Visualiser(nlink_data_1)
vis_2 = Nlink3Dmagn_Visualiser(nlink_data_2)

#Total visualiser
multiple_vis = Nlink3Dmagn_MultipleVisualiser([vis_1, vis_2])

#ax limits
ax_lims = [[-N1*l1, N1*l1], [-N1*l1, N1*l1], [-N1*l1, N1*l1]]
view_init = [30, 30]

#Animation
ani = multiple_vis.animate([zout_1, zout_2], ax_lims, view_init, interval=10, every_frame=1, colors=['b','r'])
plt.show()

#end matlab engine
eng.quit()
```

---

### Three-Sphere

```python
import os
import mpi4py
mpi4py.rc.thread_level="single"
from swimmers.utils import init_environment, Box
from swimmers.three_sphere.data import ThreeSphere2D
from swimmers.three_sphere.solve import solve_three_sphere2D
import numpy as np

#Define path
path = os.path.expanduser("~/")
print("path:", path)
```

Define swimmer and environment:

```python
#Center of the center sphere
center = np.array([4,3])

#Radius of each sphere [left,center,right]
radii = np.array([0.125,0.125,0.125])

#Arms length [left, right]
arms_length = np.array([10*0.125, 10*0.125])

#Angles of arms [left, right]*pi
angles = np.array([1., 0.])

#Box domain [bottom left point, [x length, y length]]
box = Box(np.array([0,0]), np.array([10,6]))

#Space and time discretizations
hsize_box, hsize_swimmer = 0.5, 0.06
dt, Tend = 0.05, 4

#Path of files
path_geo, path_json, path_cfg = path +"three_sphere_2D.geo", path+"three_sphere_2D.json", path +"three_sphere_2D.cfg"
directory = path

#Params of .geo, .json, .cfg
geo_param = {"box": box, "hsize_box": hsize_box, "hsize_swimmer": hsize_swimmer}
json_param = {"forceRange": hsize_swimmer, "epsBody": 5e-8, "epsWall": 5e-8}
cfg_param = {"dt": dt, "Tend": Tend, "directory": directory}

#Magnetic control
def u(t):
    def pulse(t, t0, t1):
        return np.heaviside(t-t0, 1) - np.heaviside(t-t1, 1)
    def periodic_pulse(t, t0, t1, T):
        t_mod = (t - t0) % T + t0
        return pulse(t_mod, t0, t1)
    u_left = 0.4 * (periodic_pulse(t,0,1,4)-periodic_pulse(t,2,3,4))
    u_right = 0.4 * (periodic_pulse(t,1,2,4)-periodic_pulse(t,3,4,4))
    return np.array([u_left, u_right])
```

Generate files and solve:

```python
#Define three sphere and generate files
three_sphere = ThreeSphere2D(center, radii, arms_length, angles)
three_sphere.generate_files(path_geo, path_json, path_cfg, geo_param, json_param, cfg_param)

#Feel++ environment and solve
env = init_environment(path, "three_sphere")
result = solve_three_sphere2D(three_sphere, u, False)
print(result)
```

---

### Magneto

```python
import os
import mpi4py
mpi4py.rc.thread_level="single"
from swimmers.utils import init_environment, Box
from swimmers.magneto.data import Magneto2D
from swimmers.magneto.solve import solve_magneto2D
import numpy as np

#Define path
path = os.path.expanduser("~/")
print("path:", path)
```

Define swimmer, fluid, and control:

```python
#Center and angle orientation of the heand
center, angle = np.array([0,0]), 0

#Heand size [length, height], tail length
head_length_height, tail_length = np.array([0.0005,0.0015]), 0.0075

#Fluid parameters : rho : density, mu : viscosity
rho_fluid, mu_fluid = 1260, 1.52

#Tail parameters : E :, nu :, rho : density
E_tail, nu_tail, rho_tail = 8e4, 0.4, 1200

#Head parameters : E :, nu :, rho : density
E_head, nu_head, rho_head = 4.1e10, 0.281, 7000

#Head magnetization
M_head = np.array([2e5, 2e5, 1])

#Magnetic control
freq = 0.9
ux, uy, uz = lambda t:0.005, lambda t:0.005*np.sin(2*np.pi*freq*t), lambda t:0
u = lambda t: np.array([ux(t), uy(t), uz(t)])

#Box domain [bottom left point, [x length, y length]]
box = Box(np.array([0,0]), np.array([0.04,0.02]))

#Space and time discratization
hsize_box = 0.001
hsize_swimmer = 0.0005
hfar, hclose = 0.0025, 0.00025 #for remeshing according distance to swimmer
remesh_param = np.array([hfar,hclose])
dt, Tend = 0.01, 0.03

#Path of generated files
path_geo, path_json, path_cfg = path+"magneto2D.geo", path+"magneto2D.json", path+"magneto2D.cfg"
directory = path

#Parameters for .geo, .json and .cfg
dic_fluid = {"rho": rho_fluid, "mu": mu_fluid}
dic_tail = {"E": E_tail, "nu": nu_tail, "rho": rho_tail}
dic_head = {"E": E_head, "nu": nu_head, "rho": rho_head, "M": M_head}
geo_param = {"box": box, "hsize_box": hsize_box, "hsize_swimmer":hsize_swimmer}
json_param = {"fluid": dic_fluid, "tail": dic_tail, "head": dic_head, "remesh": remesh_param}
cfg_param = {"dt": dt, "Tend": Tend, "directory": directory}
```

Generate files and solve:

```python
#Define magneto and generate files
magneto = Magneto2D(center, angle, head_length_height, tail_length)
magneto.generate_files(path_geo, path_json, path_cfg, geo_param, json_param, cfg_param)

#Feel++ environment and solve
env = init_environment(path, "magneto")
result = solve_magneto2D(magneto, u, False, ".")
print(result)
```
