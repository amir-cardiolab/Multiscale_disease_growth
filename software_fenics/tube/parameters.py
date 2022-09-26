from Flags import *
from dolfin import *


#time stepping parameters
beta2 = 0.25 #Newmark
gamma = 0.5 #Newmark
dt = Constant(5e-4) #2e-4
rho_inf = 0.0 # 0.5 #0 annihilates high freq in one step. 1: high freq preserved


#generalaized alpha method time integration
alpha_m =  (2.0* rho_inf - 1.0 ) / (1.0 + rho_inf)
alpha_f = rho_inf / ( 1.0 + rho_inf)
if (Flag_generalized_alpha):
  gamma = 0.5 - alpha_m + alpha_f
  beta2 = 0.25*( ( 1.0 - alpha_m + alpha_f )**2 )



Kappa_incomp=1e5
######### Material parameters   ################


def density(rho):
	return Constant(rho)


def C_1(ph):
	return Constant(ph)


def Dev(mytensor, C_e):
    return mytensor - 1.0/3.0*(inner(mytensor, C_e))*inv(C_e)
########################################
def update(ub,ub0):
  ub_vec, ub0_vec = ub.vector() , ub0.vector()
  ub0.vector()[:] = ub.vector()











