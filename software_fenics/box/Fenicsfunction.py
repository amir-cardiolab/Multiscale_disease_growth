from growth_cont import *
from system_ODE import *
from scipy.integrate import ode 


count=0

dt = tspan[1]/1000. 


while count<num_cycle:
	## saving files per cycle
	if count%2==0:
		tra.rename('trace','Label')
		ub.rename('disp','Label')
		eqSigma.rename("stress","Label")
		xdmffile_u.write(ub, count)
		xdmffile_u.write(tra,count)
		xdmffile_u.write(eqSigma,count)

	# Solving system of ODEs, for all nodes of the domain seprately. 
	for i in range (loc_node(mesh)):
		t=[]
		r = ode(ODEfunc).set_integrator('vode', method='adams',\
		 order=10, rtol=0, atol=1e-6, with_jacobian=True) 
		r.set_initial_value(init[i,:],tspan[0]).set_f_params(tau,ymax,w,n,EC50,growth_rate_numpy[i]) 
		results = np.empty([0,4]) 
		while r.successful() and r.t <= tspan[1]: 
			r.integrate(r.t + dt) 
			results = np.append(results,[r.y],axis=0) 
			t.append(r.t) 
		init[i,0]=results[-1,0]
		init[i,1]=results[-1,1]
		init[i,2]=results[-1,2]
		init[i,3]=results[-1,3]
		F_1[i]=(results[-1,0])
		F_2[i]=(results[-1,1])
		F_3[i]=(results[-1,2])
		F_4[i]=(results[-1,3])


	growth_rate.vector().set_local(F_3)
	growth_rate.vector().apply('insert')
	solver_isotropic.solve()
	stress_tensor = (1./J_e)* F_e*second_PK_stress*F_e.T
	Tensor_space=TensorFunctionSpace(mesh, 'CG', 1)
	eqSigma = Function(Tensor_space)
	eqSigma.assign(project(stress_tensor,Tensor_space))
	trace=tr(eqSigma)


	count+=1




B_sp=Function(Vbb)
C_sp=Function(Vbb)
D_sp=Function(Vbb)
E_sp=Function(Vbb)

while count<20:
	if count%1==0:
		xdmffile_u.write(ub, tt )
		xdmffile_u.write(B_sp,tt)
		xdmffile_u.write(C_sp,tt)
		xdmffile_u.write(D_sp,tt)
		xdmffile_u.write(E_sp,tt)
		xdmffile_u.write(eqSigma, tt)
		xdmffile_u.write(tra, tt)
	for i in range (local_dim):
		t=[]
		r = ode(ODEfunc).set_integrator('vode', method='adams',\
		 order=10, rtol=0, atol=1e-6, with_jacobian=True) 
		r.set_initial_value(init[i,:],tspan[0]).set_f_params(tau,ymax,w,n,EC50,growth_rate_numpy[i]) 
		results = np.empty([0,4]) 
		while r.successful() and r.t <= tspan[1]: 
			r.integrate(r.t + dt) 
			results = np.append(results,[r.y],axis=0) 
			t.append(r.t) 
		init[i,0]=results[-1,0]
		init[i,1]=results[-1,1]
		init[i,2]=results[-1,2]
		init[i,3]=results[-1,3]
	tt+=t[-1]
	
	B_sp.vector().set_local(init[:,0])
	B_sp.vector().apply('insert')
	C_sp.vector().set_local(init[:,1])
	C_sp.vector().apply('insert')
	D_sp.vector().set_local(init[:,2])
	D_sp.vector().apply('insert')
	E_sp.vector().set_local(init[:,3])
	E_sp.vector().apply('insert')
	if MPI.rank(MPI.comm_world) == 0:
		print("time=",count)
	growth_rate.assign(D_sp/1e8)
	growth_rate.vector().apply('insert')
	solver_isotropic.solve()
	stress_tensor = (1./J_e)* F_e*second_PK_stress*F_e.T
	TTTT=TensorFunctionSpace(bmesh, 'CG', 1)
	eqSigma = Function(TTTT)
	eqSigma.assign(project(stress_tensor,TTTT))
	# eqSigma.rename("stress","Label")
	trace=tr(eqSigma)
	tra=project(trace,Vbb)
	eqSigma.rename("stress","Label")
	tra.rename('trace','Label')
	# tra.rename('trace','Label')
