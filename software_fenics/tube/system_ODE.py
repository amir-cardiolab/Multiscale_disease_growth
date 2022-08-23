import numpy as np 
from scipy.integrate import ode 
from Geo import *
#Defining the system of ODEs. 





nn=1.1
tau1=1
ymax1=1
w1=1
EC=0.5
def loadParams(): 
	#species parameters 
	speciesNames = ['MAPK','Smad','MMP','TIMP',] 
	tau = np.array([tau1,tau1,tau1,tau1]) 
	ymax = np.array([ymax1,ymax1,ymax1,ymax1]) 
	y0 = np.array([0, 0., 0., 0.,])
 
	# reaction parameters 
	w = np.array([w1,w1,w1,w1, ]) 
	n = np.array([nn,nn,nn,nn,]) 
	EC50 = np.array([EC,EC,EC,EC ]) 
	return speciesNames, tau, ymax, y0, w, n, EC50 

 




def ODEfunc(t,y,tau,ymax,w,n,EC50,z): 
	A = 0
	B = 1
	C = 2
	D = 3
	dydt = np.zeros(4) 
	dydt[A] = (act(z,w[A],n[A],EC50[A])*ymax[A] - y[A])/tau[A]
	dydt[B] = (act(z,w[A],n[A],EC50[A])*ymax[B] - y[B])/tau[B] 
	dydt[C] = (ANDNOT(act(y[D],w[D],n[D],EC50[D]),act(y[A],w[A],n[A],EC50[A]))*ymax[C]*w[B])/tau[C]
	dydt[D] = (act(y[B],w[C],n[B],EC50[B])*ymax[D]-y[D])/tau[D]
	return dydt 

# utility functions
 
def act(x, w, n, EC50): 
	# hill activation function with parameters w (weight), n (Hill coeff), EC50 
	beta = ((EC50**n)-1)/(2*EC50**n-1) 
	K = (beta-1)**(1/n) 
	fact = w*(beta*x**n)/(K**n+x**n) 
	# if fact > w: 
	# 	fact = w 
	return fact
 
def inhib(x, w, n, EC50): 
	# inverse hill function with parameters w (weight), n (Hill coeff), EC50 
	finhib = w - act(x, w, n, EC50) 
	return finhib
 
def OR(x, y): 
	# OR logic gate 
	z = x + y - x*y 
	return z
 
def AND(w, reactList): 
	# AND logic gate, multiplying all of the reactants together 
	if w == 0: 
		z = 0 
	else: 
		p = np.array(reactList).prod() 
		z = p/w**(len(reactList)-2) 
	return z

def ANDNOT(x,y):
	z=y*(1-x)
	return z

[speciesNames, tau, ymax, y0, w, n, EC50] = loadParams()



### Time span of system of ODEs to be integrated in one cycle. 
year=365*24*3600/2
tspan = [0, year] 


####  initial Condition for all species in system of ODEs. 
init=np.zeros((loc_node(mesh),4))




#### Final results of System of ODEs in the time timespan. 
F_1=np.zeros((loc_node(mesh)))
F_2=np.zeros((loc_node(mesh)))
F_3=np.zeros((loc_node(mesh)))
F_4=np.zeros((loc_node(mesh)))
F_5=np.zeros((loc_node(mesh)))



### Growth rate function. 
growth_rate_numpy=np.zeros(loc_node(mesh))
Z_POINT=np.zeros(loc_node(mesh))
coordinate = coord(mesh)
for i in range (loc_node(mesh)):
    x_point=coordinate[i][0]
    y_point=coordinate[i][1]
    z_point=coordinate[i][2]
    growth_rate_numpy[i]=exp(-abs(z_point**3)/0.3)/2
    Z_POINT[i]=z_point








