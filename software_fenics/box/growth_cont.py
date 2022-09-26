from Flags import *
from parameters import *
from Geo import *

'''
Continuum equation implementation. Weak form of all the equations of Continuum mechanics. 
'''

# Trial and test functions


du_b = TrialFunction(FS_vec(mesh))   # displacement
vb = TestFunction(FS_vec(mesh))
ub = Function(FS_vec(mesh)) # the most recently computed solution
vel0, a0 , ub0 = Function(FS_vec(mesh)), Function(FS_vec(mesh)), Function(FS_vec(mesh))   ### initial condition of velocity displacement and acceleration. 
n_function = FacetNormal(mesh)
h = CellDiameter(mesh)


vel = (gamma/(beta2*dt))*(ub - ub0) - (gamma/beta2 - 1.0)*vel0 - dt*(gamma/(2.0*beta2) - 1.0)*a0 #vel and acceleration at t_(n+1)
a  = (1.0/(beta2*dt**2))*(ub - ub0 - dt*vel0) - (1.0/(2.0*beta2) - 1.0) * a0


if (Flag_generalized_alpha ): #generalized alpha integrator
 a = (1.0 - alpha_m) * a + alpha_m * a0
 vel = (1.0 - alpha_f) * vel + alpha_f * vel0



I = Identity(3)   # define second order identity
if (Flag_generalized_alpha ): #generalized alpha integrator
 F_t = I + grad( (1.0 - alpha_f)* ub + alpha_f * ub0)   # Total deformation gradient tensor. 
else: #Newmark or typical integrators
 F_t = I + grad(ub)


growth_rate=as_tensor([[0, 0, 0 ],[0, 0, 0],[0,0,0]])  ###### Growth tensor definition. for any problem you should change the growth tensor rate and direction yourself. you can update this tensor in the final loop.


F_growth=I+growth_rate
F_growth_inv=inv(F_growth)
F_e = F_t*F_growth_inv
C_e = dot(F_e.T,F_e)   # right Cauchy's strain tensor
J_e = det(F_e)    # Jacobian of elastic tensor
J_t=det(F_t)
I1 = tr(C_e)  # first invariant of C
I2 =  0.5*( tr(C_e)**2 - tr(C_e*C_e) )  # second invariant of C
E = 0.5*(C_e - I)   ##  strain energy based on the right Cauchy's strain tensor. 


#### The nearly incompressible material defination with the penalty term Kappa####
Kappa = Function(FS_elem(mesh))
temp_array_Kappa = Kappa.vector().get_local() 
temp_array_Kappa[:] = Kappa_incomp
Kappa.vector().set_local(temp_array_Kappa)
dU_dJ = Kappa*(J_e-1)



S_vol = J_e*dU_dJ*inv(C_e)  ### Volumetric part of Second Piola stress tensor


d_psi_bar_d_C_bar = ( C_1(C10Y)) * I     ### Constitutive equation of material properties. simple neohookean hyperelastic is defined here. 


S_isc = 2.0*J_e**(-2.0/3.0)*Dev(d_psi_bar_d_C_bar, C_e)  # Isochoric part of second Piola stress tensor

second_PK_stress =  S_vol + S_isc

#########Making weak form##############

Functional_b_isotropic = inner(dot(F_e,second_PK_stress),grad(vb))*dx()
J_b_isotropic = derivative(Functional_b_isotropic, ub, du_b)

bcs_b=[]
if MPI.rank(MPI.comm_world) == 0:
    print ('Assembling solver...')
problem_isotropic = NonlinearVariationalProblem(Functional_b_isotropic, ub, bcs_b, J_b_isotropic)
solver_isotropic = NonlinearVariationalSolver(problem_isotropic)


q_degree = 5
dx = dx(metadata={'quadrature_degree': q_degree})


#### Defining Cauchy's stress and trace of Cauchy's stress
stress_tensor = (1./J_e)* F_e*second_PK_stress*F_e.T
Tensor_F=TensorFunctionSpace(mesh, 'CG', 1)
eqSigma = Function(Tensor_F)
eqSigma.assign(project(stress_tensor,Tensor_F))
eqSigma.rename("stress","Label")
trace=tr(eqSigma)
tra=project(trace,FS_node(mesh))







