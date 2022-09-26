import numpy as np 
from scipy.integrate import ode 
from Geo import *
#Defining the system of ODEs. 




def ODEfunc(t,y,z): 
	A = 0 
	B = 1 
	C = 2 
	dydt = np.zeros(3) 
	dydt[A] = z-y[A]*y[B]
	dydt[B] = -y[A]*y[B]
	dydt[C] =y[A]*y[B]-y[C]
	return dydt 

# utility functions


def initial(z):

    y=np.array([0.1,0.1,0])
    return y




### Time span of system of ODEs to be integrated in one cycle. 
year=365*24*3600/2
tspan = [0, year] 


####  initial Condition for all species in system of ODEs. 
init=np.zeros((loc_node(mesh),4))




#### Final results of System of ODEs in the time timespan. 
F_1=np.zeros((loc_node(mesh)))
F_2=np.zeros((loc_node(mesh)))
F_3=np.zeros((loc_node(mesh)))



### Growth rate function. 
growth_rate_numpy=np.zeros(loc_node(mesh))
Z_POINT=np.zeros(loc_node(mesh))
coordinate = coord(mesh)
for i in range (local_dim):
	x_point=coordinate[i][0]
	y_point=coordinate[i][1]
	z_point=coordinate[i][2]
	growth_rate_numpy[i]=1*exp(-abs(((abs(y_point)+abs(x_point)))**2/(.3)))
	Z_POINT[i]=z_point








