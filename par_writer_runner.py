from par_writer import *

import time

#write_miniboone(mdm=0.005,mv=0.4,proddist=["","","proton_brem"],prod_chan=["pi0_decay","eta_decay","V_decay"],partlistfile=["","",""])
rho_decay_switch=False

def miniboone_baryonic_eval(mass_arr,rho_decay_switch=False,partonic_switch=False,brem_switch=False):
    MV=mass_arr[0]
    MX=mass_arr[1]
    t0 = time.time()
    parfile="parameter_run_{0}_{1}.dat".format(str(MV),str(MX))
    #if MX/1000.0<mpi0/2.0:
    
    proddist = []
    prodchan = []
    partlistfile = []
    if MX/1000.0<mpi0/2.0 and MV<600:
        proddist.append("particle_list")
        prodchan.append("pi0_decay_baryonic")
        partlistfile.append("Source/particle_list.dat")
    if MX/1000.0<meta/2.0 and MV<900:
        proddist.append("particle_list")
        prodchan.append("eta_decay_baryonic")
        partlistfile.append("Source/particle_list_k0.dat")
    if MV/1000.0>=mrho and partonic_switch:
	proddist.append("parton_V_baryonic")
	prodchan.append("parton_production_baryonic")
	partlistfile.append("")
    if ((MV<1200) and (MV>=350)) and rho_decay_switch:
        proddist.append("particle_list")
        prodchan.append("omega_decay_baryonic")
        partlistfile.append("Source/particle_list.dat")
    if MV/1000.0>=mrho and MV<=1250 and rho_decay_switch:
        proddist.append("particle_list")
        prodchan.append("phi_decay_baryonic")
        partlistfile.append("Source/particle_list.dat")


    write_miniboone(mdm=MX/1000.0,mv=MV/1000.0,proddist=proddist,prod_chan=prodchan,partlistfile=partlistfile,outfile=parfile,signal_chan="NCE_nucleon_baryonic",eps=0,alpha_D=1e-4,dm_energy_resolution=0.01,sumlog="Events/miniboone_baryonic_paper.dat")
    subp.call(["./Source/main", parfile])
    t1 = time.time()
    print "\ntime={}\n".format(t1-t0) 
    t0 = time.time() 
    #if MX/1000.0<meta/2.0:
    #    write_miniboone(mdm=MX/1000.0,mv=MV/1000.0,proddist="particle_list",prod_chan="eta_decay",outfile=parfile)
    #    subp.call(["./Source/main", parfile])
    #    t1 = time.time()
    #    print "\ntime={}\n".format(t1-t0)
    subp.call(["rm", parfile])

def ship_eval(mass_arr):
    MV=mass_arr[0]
    MX=mass_arr[1]
    t0 = time.time()
    parfile="parameter_run_{0}_{1}.dat".format(str(MV),str(MX))
    #if MX/1000.0<mpi0/2.0:
    
    rho_decay_switch=True
    proton_brem_switch=True
    parton_switch = False;

    proddist = []
    prodchan = []
    partlistfile = []
    if MX/1000.0<mpi0/2.0 and MV<600.0:
        proddist.append("particle_list")
        prodchan.append("pi0_decay")
        partlistfile.append("Source/particle_list_ship.dat")
    if MX/1000.0<meta/2.0 and MV<900.0:
        proddist.append("particle_list")
        prodchan.append("eta_decay")
        partlistfile.append("Source/particle_list_ship.dat")
    if MV/1000.0>=mrho and parton_switch:
	proddist.append("parton_V")
	prodchan.append("parton_production")
	partlistfile.append("")
    if MV/1000.0<=1 and MV/2.0>MX and proton_brem_switch:
	proddist.append("proton_brem")
	prodchan.append("V_decay")
 	partlistfile.append("")
    if ((MV<1200) and (MV>=350)) and rho_decay_switch:
        proddist.append("particle_list")
        prodchan.append("rho_decay")
        partlistfile.append("Source/particle_list_ship.dat")
        proddist.append("particle_list")
        prodchan.append("omega_decay")
        partlistfile.append("Source/particle_list_ship.dat")
    if MV/1000.0>=mrho and MV<=1250:
        proddist.append("particle_list")
        prodchan.append("phi_decay")
        partlistfile.append("Source/particle_list_ship.dat")


    write_ship(mdm=MX/1000.0,mv=MV/1000.0,proddist=proddist,prod_chan=prodchan,partlistfile=partlistfile,outfile=parfile,signal_chan="NCE_electron",samplesize=2000,sumlog="Events/ship.dat")
    subp.call(["./Source/main", parfile])
    t1 = time.time()
    print "\ntime={}\n".format(t1-t0) 
    t0 = time.time() 
    #if MX/1000.0<meta/2.0:
    #    write_miniboone(mdm=MX/1000.0,mv=MV/1000.0,proddist="particle_list",prod_chan="eta_decay",outfile=parfile)
    #    subp.call(["./Source/main", parfile])
    #    t1 = time.time()
    #    print "\ntime={}\n".format(t1-t0)
    subp.call(["rm", parfile])


def miniboone_eval(mass_arr,rho_decay_switch=False,partonic_switch=True,brem_switch=True,signal_channel = "NCE_nucleon"):
    
    MV=mass_arr[0]
    MX=mass_arr[1]
    t0 = time.time()
    parfile="parameter_run_{0}_{1}.dat".format(str(MV),str(MX))
    #if MX/1000.0<mpi0/2.0: 

    proddist = []
    prodchan = []
    partlistfile = []
    if MX/1000.0<mpi0/2.0 and MV<600.0:
        proddist.append("particle_list")
        prodchan.append("pi0_decay")
        partlistfile.append("Source/particle_list.dat")
    if MX/1000.0<meta/2.0 and MV<900.0:
        proddist.append("particle_list")
        prodchan.append("eta_decay")
        partlistfile.append("Source/particle_list_k0.dat")
    if MV/1000.0>=mrho and partonic_switch:
	proddist.append("parton_V")
	prodchan.append("parton_production")
	partlistfile.append("")
    if MV/1000.0<=1 and MV/2.0>MX and brem_switch:
	proddist.append("proton_brem")
	prodchan.append("V_decay")
 	partlistfile.append("")
    if ((MV<1200) and (MV>=350)) and rho_decay_switch:
        proddist.append("particle_list")
        prodchan.append("rho_decay")
        partlistfile.append("Source/particle_list.dat")
        proddist.append("particle_list")
        prodchan.append("omega_decay")
        partlistfile.append("Source/particle_list.dat")
   # if MV/1000.0>=mrho and MV<=1250:
   #     proddist.append("particle_list")
   #     prodchan.append("phi_decay")
   #     partlistfile.append("Source/particle_list.dat")
    if signal_channel=="NCE_nucleon":
        write_miniboone(mdm=MX/1000.0,mv=MV/1000.0,proddist=proddist,prod_chan=prodchan,partlistfile=partlistfile,outfile=parfile,samplesize=1000,signal_chan="NCE_nucleon",sumlog="Events/NCE_nucloen_paper_miniboone.dat")
    elif signal_channel=="Pion_Inelastic":
        write_miniboone(mdm=MX/1000.0,mv=MV/1000.0,proddist=proddist,prod_chan=prodchan,samplesize=1000,min_scatter_energy=0,max_scatter_energy=1000.0,partlistfile=partlistfile,outfile=parfile,signal_chan="Pion_Inelastic",sumlog="Events/pion_inelastic_paper_miniboone.dat")
    subp.call(["./Source/main", parfile])
    t1 = time.time()
    print "\ntime={}\n".format(t1-t0) 
    t0 = time.time() 
    #if MX/1000.0<meta/2.0:
    #    write_miniboone(mdm=MX/1000.0,mv=MV/1000.0,proddist="particle_list",prod_chan="eta_decay",outfile=parfile)
    #    subp.call(["./Source/main", parfile])
    #    t1 = time.time()
    #    print "\ntime={}\n".format(t1-t0)
    subp.call(["rm", parfile])

def t2k_FD_eval(mass_arr,signal_channel="Pion_Inelastic"):
    MV=mass_arr[0]
    MX=mass_arr[1]
    t0 = time.time()
    parfile="parameter_run_{0}_{1}.dat".format(str(MV),str(MX))
    #if MX/1000.0<mpi0/2.0:
    
    rho_decay_switch=False

    proddist = []
    prodchan = []
    partlistfile = []
    if MX/1000.0<mpi0/2.0:
        proddist.append("particle_list")
        prodchan.append("pi0_decay")
        partlistfile.append("Source/particle_list.dat")
    if MX/1000.0<meta/2.0:
        proddist.append("particle_list")
        prodchan.append("eta_decay")
        partlistfile.append("Source/particle_list.dat")
    if MV/1000.0>=mrho:
	proddist.append("parton_V")
	prodchan.append("parton_production")
	partlistfile.append("")
    if MV/1000.0<=2:
	proddist.append("proton_brem")
	prodchan.append("V_decay")
 	partlistfile.append("")
    if ((MV<1200) and (MV>500)) and rho_decay_switch:
        proddist.append("particle_list")
        prodchan.append("rho_decay")
        partlistfile.append("Source/particle_list_k0.dat")
    if signal_channel=="NCE_nucleon":
        write_t2kFD(mdm=MX/1000.0,mv=MV/1000.0,samplesize=1000,proddist=proddist,prod_chan=prodchan,signal_chan="NCE_nucleon",min_scatter_energy=0.0014,max_scatter_energy=1000.0,partlistfile=partlistfile,outfile=parfile)
    elif signal_channel == "Pion_Inelastic":
        write_t2kFD(mdm=MX/1000.0,mv=MV/1000.0,samplesize=1000,proddist=proddist,signal_chan="Pion_Inelastic",min_scatter_energy=0.0,max_scatter_energy=1000.0,prod_chan=prodchan,partlistfile=partlistfile,outfile=parfile)
    elif signal_channel == "NCE_electron": 
        write_t2kFD(mdm=MX/1000.0,mv=MV/1000.0,samplesize=1000,proddist=proddist,signal_chan="NCE_electron",min_scatter_energy=0.0,max_scatter_energy=1000.0,prod_chan=prodchan,partlistfile=partlistfile,outfile=parfile)
    subp.call(["./Source/main", parfile])
    t1 = time.time()
    print "\ntime={}\n".format(t1-t0) 
    t0 = time.time() 
    #if MX/1000.0<meta/2.0:
    #    write_miniboone(mdm=MX/1000.0,mv=MV/1000.0,proddist="particle_list",prod_chan="eta_decay",outfile=parfile)
    #    subp.call(["./Source/main", parfile])
    #    t1 = time.time()
    #    print "\ntime={}\n".format(t1-t0)
    subp.call(["rm", parfile])



def t2k_eval(mass_arr,signal_channel="Pion_Inelastic"):
    MV=mass_arr[0]
    MX=mass_arr[1]
    t0 = time.time()
    parfile="parameter_run_{0}_{1}.dat".format(str(MV),str(MX))
    #if MX/1000.0<mpi0/2.0:
    
    rho_decay_switch=False

    proddist = []
    prodchan = []
    partlistfile = []
    if MX/1000.0<mpi0/2.0:
        proddist.append("particle_list")
        prodchan.append("pi0_decay")
        partlistfile.append("Source/particle_list.dat")
    if MX/1000.0<meta/2.0:
        proddist.append("particle_list")
        prodchan.append("eta_decay")
        partlistfile.append("Source/particle_list.dat")
    if MV/1000.0>=mrho:
	proddist.append("parton_V")
	prodchan.append("parton_production")
	partlistfile.append("")
    if MV/1000.0<=2:
	proddist.append("proton_brem")
	prodchan.append("V_decay")
 	partlistfile.append("")
    if ((MV<1200) and (MV>500)) and rho_decay_switch:
        proddist.append("particle_list")
        prodchan.append("rho_decay")
        partlistfile.append("Source/particle_list_k0.dat")
    if signal_channel=="NCE_nucleon":
        write_nd280(mdm=MX/1000.0,mv=MV/1000.0,proddist=proddist,prod_chan=prodchan,partlistfile=partlistfile,outfile=parfile)
    elif signal_channel == "Pion_Inelastic":
        write_nd280(mdm=MX/1000.0,mv=MV/1000.0,samplesize=1000,proddist=proddist,signal_chan="Pion_Inelastic",min_scatter_energy=0.0,max_scatter_energy=1000.0,prod_chan=prodchan,partlistfile=partlistfile,outfile=parfile)
    subp.call(["./Source/main", parfile])
    t1 = time.time()
    print "\ntime={}\n".format(t1-t0) 
    t0 = time.time() 
    #if MX/1000.0<meta/2.0:
    #    write_miniboone(mdm=MX/1000.0,mv=MV/1000.0,proddist="particle_list",prod_chan="eta_decay",outfile=parfile)
    #    subp.call(["./Source/main", parfile])
    #    t1 = time.time()
    #    print "\ntime={}\n".format(t1-t0)
    subp.call(["rm", parfile])

def lsnd_eval(mass_arr):
    MV=mass_arr[0]
    MX=mass_arr[1]
    t0 = time.time()
    parfile="parameter_run_{0}_{1}.dat".format(str(MV),str(MX))
    #if MX/1000.0<mpi0/2.0:
    
    rho_decay_switch=False

    proddist = []
    prodchan = []
    partlistfile = []
    if MX/1000.0<mpi0/2.0:
        proddist.append("particle_list")
        prodchan.append("pi0_decay")
        partlistfile.append("Source/particle_list_lsnd.dat")
    
    write_lsnd(mdm=MX/1000.0,mv=MV/1000.0,proddist=proddist,prod_chan=prodchan,partlistfile=partlistfile,outfile=parfile,sumlog="Events/lsnd_thesis.dat",p_num_target=8,samplesize=1000)
    subp.call(["./Source/main", parfile])
    t1 = time.time()
    print "\ntime={}\n".format(t1-t0) 
    t0 = time.time() 
    subp.call(["rm", parfile])
'''
#LSND weird
def lsnd_eval(mass_arr):
    MV=mass_arr[0]
    MX=mass_arr[1]
    t0 = time.time()
    parfile="parameter_run_{0}_{1}.dat".format(str(MV),str(MX))
    #if MX/1000.0<mpi0/2.0:
    
    rho_decay_switch=False

    proddist = []
    prodchan = []
    partlistfile = []
    if MX/1000.0<mpi0/2.0:
        proddist.append("particle_list")
        prodchan.append("pi0_decay_baryonic")
        partlistfile.append("Source/particle_list_lsnd.dat")
    
    write_lsnd(mdm=MX/1000.0,mv=MV/1000.0,proddist=proddist,prod_chan=prodchan,partlistfile=partlistfile,outfile=parfile,alpha_D=1e-6,eps=2e-4)
    subp.call(["./Source/main", parfile])
    t1 = time.time()
    print "\ntime={}\n".format(t1-t0) 
    t0 = time.time() 
    #if MX/1000.0<meta/2.0:
    #    write_miniboone(mdm=MX/1000.0,mv=MV/1000.0,proddist="particle_list",prod_chan="eta_decay",outfile=parfile)
    #    subp.call(["./Source/main", parfile])
    #    t1 = time.time()
    #    print "\ntime={}\n".format(t1-t0)
    subp.call(["rm", parfile])
'''

def execute_miniboone_parallel(genlist=True):
    if genlist:
        write_miniboone(prod_chan=["pi0_decay"],proddist=["pi0_sanfordwang"],samplesize=2e6,output_mode="particle_list",partlistfile=["Source/particle_list.dat"])
        subp.call(["./Source/main", "parameter_run.dat"])
        write_miniboone(prod_chan=["eta_decay"],proddist=["k0_sanfordwang"],samplesize=2e6,output_mode="particle_list",partlistfile=["Source/particle_list_k0.dat"])
        subp.call(["./Source/main", "parameter_run.dat"])
        #write_miniboone(prod_chan=["V_decay"],proddist=["proton_brem"],samplesize=4e6,output_mode="particle_list",partlistfile=["Source/particle_list_miniboone_brem.dat"])
    #vmassarr=[1,5]+[10*i for i in xrange(1,14)]+[10*i for i in xrange(15,100,2)]+[775,774,776,777,778,779,781,782,783,785,787]+[1005,1010,1015,1.02,1025,1030,1035,1050,1100,1150,1200,1300]+[i for i in xrange(1000,2100,250)]
    #chimassarr=[1,5]+[10*i for i in xrange(1,27)]
    #vmassarr=[600,700,800,900]
    #vmassarr=[300]
    chimassarr=[10]
    vmassarr=[15,16,17,18,19,21,22,23,24,25]
    #vmassarr=vmassarr+[10*i for i in xrange(1,14)]+[10*i for i in xrange(15,90,2)] 
    vmassarr=vmassarr+[i for i in xrange(20,140,10)]+[i for i in xrange(150,1000,50)]+[770,22,24]+[1000,1100,1150,1200,1300]+[i for i in xrange(1000,4100,250)]
   #+[775,774,776,777,778,779,781,782,783,785,787]
    #chimassarr=[i for i in xrange(10,270,10)]+[132,134,136]+[1,5]
    #chimassarr=[173,175,178]
    massarr=[[MV,MX] for MV in vmassarr for MX in chimassarr]
    for marr in massarr:
        miniboone_eval(marr)
    #pool=Pool(processes=3)
    #pool.map(miniboone_eval,massarr) 

def execute_ship_parallel(genlist=True):
    if genlist:
        #write_ship(prod_chan=["pi0_decay"],proddist=["bmpt"],samplesize=1e6,output_mode="particle_list",partlistfile=["Source/particle_list_ship.dat"])
        write_ship(prod_chan=["parton_production"],proddist=["parton_V"],samplesize=2e5,output_mode="particle_list",partlistfile=["Source/particle_list_ship_parton.dat"])
        #write_ship(prod_chan=["V_decay"],proddist=["proton_brem"],samplesize=4e6,output_mode="particle_list",partlistfile=["Source/particle_list_ship_brem.dat"])
        subp.call(["./Source/main", "parameter_run.dat"])
    vmassarr=[i for i in xrange(100,1000,50)]+[775,774,776,777,778,779,781,782,783,785,787]+[1005,1010,1015,1.02,1025,1030,1035,1050,1100,1150,1200,1300]+[i for i in xrange(1000,4100,250)]+[390,395,405,410]
    #chimassarr=[1,5]+[10*i for i in xrange(1,27)]
    #vmassarr=[600,700,800,900]
    chimassarr=[200]
    massarr=[[MV,MX] for MV in vmassarr for MX in chimassarr]
    #for marr in massarr:
    #    ship_eval(marr)
    #pool=Pool(processes=4)
    #pool.map(ship_eval,massarr) 

def execute_miniboone_baryonic_parallel(genlist=True):
    if genlist:
        write_miniboone(prod_chan=["pi0_decay"],proddist=["pi0_sanfordwang"],samplesize=2e6,output_mode="particle_list",partlistfile=["Source/particle_list.dat"])
        subp.call(["./Source/main", "parameter_run.dat"])
        write_miniboone(prod_chan=["eta_decay"],proddist=["k0_sanfordwang"],samplesize=2e6,output_mode="particle_list",partlistfile=["Source/particle_list_k0.dat"])
        subp.call(["./Source/main", "parameter_run.dat"])
    #vmassarr=[22,25,27]+[10*i for i in xrange(3,14)]+[10*i for i in xrange(15,100,5)]+[775,774,776,777,778,779,781,782,783,785,787]+[1005,1010,1015,1.02,1025,1030,1035,1050,1100,1150,1200,1300]+[i for i in xrange(1000,2100,250)]
    #chimassarr=[1,5]+[10*i for i in xrange(1,27)]
    #vmassarr=[600,700,800,900]
    vmassarr=[10*i for i in xrange(3,14)]+[10*i for i in xrange(15,90,2)]+[21,23,25]
    chimassarr=[10]    
    massarr=[[MV,MX] for MV in vmassarr for MX in chimassarr]
    vmassarr=[300]
    chimassarr=[1,5]+[10*i for i in xrange(1,14)]
    massarr=massarr+[[MV,MX] for MV in vmassarr for MX in chimassarr]
    #for marr in massarr:
    #    miniboone_baryonic_eval(marr)
    pool=Pool(processes=3)
    pool.map(miniboone_baryonic_eval,massarr) 

import numpy as np

def execute_lsnd_parallel(genlist=True):
    if genlist:
        write_lsnd(prod_chan=["pi0_decay"],proddist=["burmansmith"],samplesize=0.55e6,output_mode="particle_list",partlistfile=["Source/particle_list_lsnd_low.dat"],p_num_target=8)
        subp.call(["./Source/main","parameter_run.dat"])
        write_lsnd(prod_chan=["pi0_decay"],proddist=["burmansmith"],samplesize=0.45e6,output_mode="particle_list",partlistfile=["Source/particle_list_lsnd_high.dat"],p_num_target=70)
        subp.call(["./Source/main","parameter_run.dat"])
        arr1 = np.loadtxt("Source/particle_list_lsnd_low.dat")
        arr2 = np.loadtxt("Source/particle_list_lsnd_high.dat")
        arr3 = np.append(arr1,arr2,axis=0)
        np.random.shuffle(arr3)
        np.savetxt("particle_list_lsnd",arr3)
    massarr=[]
    #chimassarr=[5,10]
    #vmassarr=[1,2,3,4,5,6,7,8,10,12,14,16,18,19,20,21,22,25,27]+[i for i in xrange(30,130,10)]+[132,134,136,138]+[i for i in xrange(140,700,20)]
    #massarr=[[MV,MX] for MV in vmassarr for MX in chimassarr]
    #vmassarr=[i for i in xrange(600,910,50)]
    #chimassarr=[10]
    vmassarr=[20,50,100]
    #chimassarr=[i for i in xrange(15,100,10)]+[11,12,13,14,16,18,23,27]
    chimassarr=[5,7,9,49,51,24,26]
    massarr=massarr+[[MV,MX] for MV in vmassarr for MX in chimassarr]
    for marr in massarr:
        lsnd_eval(marr)
    #pool=Pool(processes=4)
    #pool.map(lsnd_eval,massarr) 

def execute_nd280_parallel(genlist=True):
    if genlist:
        write_nd280(prod_chan=["pi0_decay"],proddist=["bmpt"],samplesize=1e6,output_mode="particle_list",partlistfile=["Source/particle_list_nd280.dat"])
        #write_ship(prod_chan=["parton_production"],proddist=["parton_V"],samplesize=2e5,output_mode="particle_list",partlistfile=["Source/particle_list_ship_parton.dat"])
        #write_ship(prod_chan=["V_decay"],proddist=["proton_brem"],samplesize=4e6,output_mode="particle_list",partlistfile=["Source/particle_list_ship_brem.dat"])
        subp.call(["./Source/main", "parameter_run.dat"])
    #vmassarr=[i for i in xrange(20,140,10)]+[i for i in xrange(150,1000,50)]+[770,22,24]+[1000,1100,1150,1200,1300]+[i for i in xrange(1000,4100,250)]+[390,395,405,410]
    #chimassarr=[1,5]+[10*i for i in xrange(1,27)]
    #vmassarr=[600,700,800,900]
    chimassarr=[100]
    vmassarr=[300]
    massarr=[[MV,MX] for MV in vmassarr for MX in chimassarr]
    for marr in massarr:
        #t2k_eval(marr)
        t2k_FD_eval(marr)
    #pool=Pool(processes=4)
    #pool.map(ship_eval,massarr)

execute_nd280_parallel(genlist=True)
#execute_miniboone_baryonic_parallel(genlist=True)
#execute_miniboone_parallel(genlist=False)
#execute_ship_parallel(genlist=True)
#execute_lsnd_parallel(genlist=False)
