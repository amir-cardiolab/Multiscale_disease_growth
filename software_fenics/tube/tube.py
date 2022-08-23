from dolfin import *
from mshr import *

'''
Example of Tube.  Defining the geometry, material properties and number of cycles for growth and system of ODEs coupling. if you face numerical challenges, you should increase the number
of cycles and reduce the growth rate in the final loop. 
'''
x0=0
y0 = 0
z0 = -5
x1 = .25
y1 = .25
z1 = 5
n_x = 20
n_y = 20
n_z = 50
geometry1 = Cylinder(Point(x0, y0, z0),Point(x0,y0,z1), .65,.65)
geo2=Cylinder(Point(x0, y0, z0),Point(x0,y0,z1),0.5,0.5)
geometry=geometry1-geo2
mesh = generate_mesh(geometry, 120)

left_x = CompiledSubDomain("abs(pow(x[0],2)+pow(x[1],2)-0.25)<1e-1  && on_boundary")
right_x = CompiledSubDomain("abs(pow(x[0],2)+pow(x[1],2)-.4225)<1e-1 && on_boundary")
top_z = CompiledSubDomain("near(x[2],-5) && on_boundary")
bottom_z = CompiledSubDomain("near(x[2],5) && on_boundary")


#### Defining the boundary of the tube surfaces. 
boundary_indices = {'left_x': 1,
                    'right_x': 2,
                    'top_z': 3,
                    'bottom_z': 4,
}

### Mark the surface boundaries.
boundary_markers = MeshFunction("size_t", mesh, dim=2, value=0)
left_x.mark(boundary_markers, boundary_indices["left_x"])
right_x.mark(boundary_markers, boundary_indices["right_x"])
top_z.mark(boundary_markers, boundary_indices["top_z"])
bottom_z.mark(boundary_markers, boundary_indices["bottom_z"])



eps = 5e-10      # diffusion coefficient for system of PDEs. 

######### Material parameters   ################

rho = 1100.0  #material density
Kappa_incomp = 1e6 # Penalty term for the nearly incompressible material. 
delta = 0.558
C10Y = 1e7 / 6.  # Stiffness of the tissue. 
num_cycle=100


xdmffile_u = XDMFFile('tube'+'.xdmf')