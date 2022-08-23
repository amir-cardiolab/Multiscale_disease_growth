
from growth_cont import solver_isotropic





#solver.parameters.update(snes_solver_parameters)
prm = solver_isotropic.parameters
#print prm["newton_solver"]["maximum_iterations"]
prm["newton_solver"]["maximum_iterations"] = 200
#prm["newton_solver"]["relaxation_parameter"] = 0.6 #.7
prm["newton_solver"]["krylov_solver"]["nonzero_initial_guess"] = False
# prm["newton_solver"]["krylov_solver"]["absolute_tolerance"] = 1e-8 #default is 1e-10
prm['newton_solver']['absolute_tolerance'] = 1E-9 #1E-5
#prm['newton_solver']['relative_tolerance'] = 1E-1
prm['newton_solver']['linear_solver'] = 'tfqmr' #'gmres' bicgstab
prm['newton_solver']['preconditioner'] = 'petsc_amg'
prm["newton_solver"]["krylov_solver"]["monitor_convergence"] =False
prm['newton_solver']['krylov_solver']['maximum_iterations'] = 20000

