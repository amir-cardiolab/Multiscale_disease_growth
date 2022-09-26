from dolfin import *
from mshr import *

'''
Example of Box.  Defining the geometry, material properties and number of cycles for growth and system of ODEs coupling. if you face numerical challenges, you should increase the number
of cycles and reduce the growth rate in the final loop. 
'''



x0=-0.5
y0 = -.5
z0 = -.5

#End mesh values
x1 =0.5
y1 = 0.5
z1 = 0.5


n_x = 7
n_y = 7
n_z = 7
mesh = BoxMesh(Point(x0, y0, z0),Point(x1,y1,z1), n_x,n_y,n_z)

left_x = CompiledSubDomain("near(x[0],-0.5) && on_boundary")
right_x = CompiledSubDomain("near(x[0],0.5) && on_boundary")
left_y = CompiledSubDomain("near(x[1],-0.5) && on_boundary")
right_y = CompiledSubDomain("near(x[1],0.5) && on_boundary")
left_z = CompiledSubDomain("near(x[2],-0.5) && on_boundary")
right_z = CompiledSubDomain("near(x[2],0.5) && on_boundary")


boundary_indices = {'left_x': 1,
                    'right_x': 2,
                    'left_y': 3,
                    'right_y': 4,
                    'left_z': 5,
                    'right_z': 6

}

boundary_markers = MeshFunction("size_t", mesh, dim=2, value=0)
left_x.mark(boundary_markers, boundary_indices["left_x"])
right_x.mark(boundary_markers, boundary_indices["right_x"])
left_y.mark(boundary_markers, boundary_indices["left_y"])
right_y.mark(boundary_markers, boundary_indices["right_y"])
left_z.mark(boundary_markers, boundary_indices["left_z"])
right_z.mark(boundary_markers, boundary_indices["right_z"])



eps = 5e-10      # diffusion coefficient for system of PDEs. 

######### Material parameters   ################

rho = 1100.0  #material density
Kappa_incomp = 1e6 # Penalty term for the nearly incompressible material. 
delta = 0.558
C10Y = 1e7 / 6.  # Stiffness of the tissue. 
num_cycle=100   # number of cycles.


xdmffile_u = XDMFFile('box'+'.xdmf')


