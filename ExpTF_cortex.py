from brian2 import *
prefs.codegen.target = "numpy"

start_scope()

DT=0.1 # time step
defaultclock.dt = DT*ms
N_inh = 2000 # number of inhibitory neurons
N_exc = 8000 # number of excitatory neurons

TotTime=1000 #Simulation duration (ms)
duration = TotTime*ms


FRout_inh=[]
FRout_exc=[]
Npts = 20  #Resolution of the measured TF
i = 0
j = 1
for rate_exc in linspace(0, 150, Npts):
	FRout_inh.append([])
	FRout_exc.append([])
	for rate_inh in linspace(0, 150, Npts):
		
		eqs='''
		dv/dt = (-GsynE*(v-Ee)-GsynI*(v-Ei)-gl*(v-El)+ gl*Dt*exp((v-Vt)/Dt)-w + Is)/Cm : volt (unless refractory)
		dw/dt = (a*(v-El)-w)/tau_w:ampere
		dGsynI/dt = -GsynI/Ti : siemens
		dGsynE/dt = -GsynE/Te : siemens
		Is:ampere
		Cm:farad
		gl:siemens
		El:volt
		a:siemens
		tau_w:second
		Dt:volt
		Vt:volt
		Ee:volt
		Ei:volt
		Te:second
		Ti:second
		'''

		# Population 1 - FS - Inhibitory

		b_FS = 0.*pA
		G_inh = NeuronGroup(N_inh, eqs, threshold='v > -20*mV', reset='v = -60*mV; w += b_FS', refractory='5*ms', method='heun')
		#init:
		G_inh.v = -65*mV
		G_inh.w = 0.0*pA
		# synaptic
		G_inh.GsynI=0.0*nS
		G_inh.GsynE=0.0*nS
		G_inh.Ee=0.*mV
		G_inh.Ei=-80.*mV
		G_inh.Te=4.*ms
		G_inh.Ti=10.*ms
		# cell parameters
		G_inh.Cm = 200.*pF
		G_inh.gl = 10.*nS
		G_inh.Vt = -50.*mV
		G_inh.Dt = 0.5*mV
		G_inh.tau_w = 200*ms #1.0*ms
		G_inh.Is = 0.0 #[0.0 for i in range(N1)]*nA
		G_inh.El = -65.*mV
		G_inh.a = 0.0*nS


		# Population 2 - RS - Excitatory 
		
		b_RS = 40.*pA
		G_exc = NeuronGroup(N_exc, eqs, threshold='v > -20*mV', reset='v = -60*mV; w += b_RS', refractory='5*ms',  method='heun')
		G_exc.v = -65.*mV
		G_exc.w = 0.0*pA
		# synaptic
		G_exc.GsynI=0.0*nS
		G_exc.GsynE=0.0*nS
		G_exc.Ee=0.*mV
		G_exc.Ei=-80.*mV
		G_exc.Te=4.*ms
		G_exc.Ti=10.*ms
		# cell
		G_exc.Cm = 200.*pF
		G_exc.gl = 10.*nS
		G_exc.Vt = -50.*mV
		G_exc.Dt = 2.*mV
		G_exc.tau_w = 200*ms #1000.*ms
		G_exc.Is = 0.0*nA #2.50*nA #[0.0 for i in range(N2)]*nA
		#G_exc.Is[0]=0.*nA
		G_exc.El = -65.*mV
		G_exc.a = 4.*nS


		# external drive--------------------------------------------------------------------------

		P_ed_inh = PoissonGroup(N_inh, rates=rate_inh*Hz)
		P_ed_exc = PoissonGroup(N_exc, rates=rate_exc*Hz)


		# Network-----------------------------------------------------------------------------

		# quantal increment in synaptic conductances:
		Qe = 1.5*nS
		Qi = 5.*nS

		# probability of connection
		prbC= 0.05

		
		S_edin_ex = Synapses(P_ed_exc, G_inh, on_pre='GsynE_post+=Qe')
		S_edin_ex.connect(p=prbC)

		S_edex_in = Synapses(P_ed_inh, G_exc, on_pre='GsynI_post+=Qi')
		S_edex_in.connect(p=prbC)

		S_edin_in = Synapses(P_ed_inh, G_inh, on_pre='GsynI_post+=Qi')
		S_edin_in.connect(p=prbC)

		S_edex_ex = Synapses(P_ed_exc, G_exc, on_pre='GsynE_post+=Qe')
		S_edex_ex.connect(p=prbC)


		# Recording tools -------------------------------------------------------------------------------

		FRG_inh = PopulationRateMonitor(G_inh)
		FRG_exc = PopulationRateMonitor(G_exc)


		# Run simulation -------------------------------------------------------------------------------

		run(duration)

		# Plots -------------------------------------------------------------------------------



		# prepare firing rate
		def bin_array(array, BIN, time_array):
		    N0 = int(BIN/(time_array[1]-time_array[0]))
		    N1 = int((time_array[-1]-time_array[0])/BIN)
		    return array[:N0*N1].reshape((N1,N0)).mean(axis=1)

		BIN=5
		time_array = arange(int(TotTime/DT))*DT



		LfrG_exc=array(FRG_exc.rate/Hz)
		TimBinned,popRateG_exc=bin_array(time_array, BIN, time_array),bin_array(LfrG_exc, BIN, time_array)

		LfrG_inh=array(FRG_inh.rate/Hz)
		TimBinned,popRateG_inh=bin_array(time_array, BIN, time_array),bin_array(LfrG_inh, BIN, time_array)

		#print(mean(popRateG_inh[150::]))
		#print(mean(popRateG_exc[150::]))
		FRout_inh[i].append(mean(popRateG_inh[150::]))
		FRout_exc[i].append(mean(popRateG_exc[150::]))

		print(f' {j}/{Npts**2}', end='\r')
		j=j+1
	i=i+1


np.save(f'data\\ExpTF_cortex_inh_Nstp{Npts}.npy', FRout_inh)
np.save(f'data\\ExpTF_cortex_exc_Nstp{Npts}.npy', FRout_exc)
plt.figure()
plt.imshow(FRout_inh)
plt.figure()
plt.imshow(FRout_exc)
plt.show()

print()
