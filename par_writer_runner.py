from par_writer import *
import numpy as np
import sys

import time

_pion_decay=1
_eta_decay=2
_rho_decay=3
_parton=4
_brem=5

#write_miniboone(mdm=0.005,mv=0.4,proddist=["","","proton_brem"],prod_chan=["pi0_decay","eta_decay","V_decay"],partlistfile=["","",""])
rho_decay_switch=False

def miniboone_baryonic_eval(mass_arr,rho_decay_switch=False,partonic_switch=True,brem_switch=True,det_switch="miniboone"):
    MV=mass_arr[0]
    MX=mass_arr[1]
    t0 = time.time()
    parfile="parameter_run_{0}_{1}.dat".format(str(MV),str(MX))
    #if MX/1000.0<mpi0/2.0:

    zmin=max(3*MV/1000.0/8.9,0.3)
    zmax=1-zmin

    print("Zmin=",zmin,"Zmax=",zmax)

    proddist = []
    prodchan = []
    partlistfile = []
    if MX/1000.0<mpi0/2.0 and MV<600:
        proddist.append("particle_list")
        prodchan.append("pi0_decay_baryonic")
        partlistfile.append("data/particle_list.dat")
    if MX/1000.0<meta/2.0 and MV<900:
        proddist.append("particle_list")
        prodchan.append("eta_decay_baryonic")
        partlistfile.append("data/particle_list_k0.dat")
    if MV/1000.0>=mrho and partonic_switch:
	   proddist.append("parton_V_baryonic")
	   prodchan.append("parton_production_baryonic")
	   partlistfile.append("")
    if ((MV<1200) and (MV>=350)) and rho_decay_switch:
        proddist.append("particle_list")
        prodchan.append("omega_decay_baryonic")
        partlistfile.append("data/particle_list.dat")
    if zmin<zmax and MV/2.0>MX and brem_switch:
        proddist.append("proton_brem_baryonic")
        prodchan.append("V_decay_baryonic")
        partlistfile.append("")
    #if MV/1000.0>=mrho and MV<=1250:
    #    proddist.append("particle_list")
    #    prodchan.append("phi_decay_baryonic")
    #    partlistfile.append("data/particle_list.dat")
    if det_switch == "sbnd":
        print("SBND run in progress")
        write_sbnd(mdm=MX/1000.0,mv=MV/1000.0,proddist=proddist,samplesize=1000,prod_chan=prodchan,partlistfile=partlistfile,outfile=parfile,signal_chan="NCE_nucleon_baryonic",eps=0,alpha_D=1e-4,dm_energy_resolution=0.01,sumlog="Events/sbnd_y4.dat",output_mode="summary",zmin=zmin,zmax=zmax)
    else:
        #write_miniboone(mdm=MX/1000.0,mv=MV/1000.0,proddist=proddist,prod_chan=prodchan,partlistfile=partlistfile,min_scatter_energy=0.05,max_scatter_energy=2,outfile=parfile,samplesize=1000,eps=0,alpha_D=1e-5,signal_chan="NCE_nucleon_baryonic",sumlog="Events/summary_baryonic_test.dat",output_mode="summary",outlog="Events/mini_baryonic_{}_{}.dat".format(MV,MX),det=miniboone_detector_full,efficiency=1)
        #write_miniboone(mdm=MX/1000.0,mv=MV/1000.0,proddist=proddist,samplesize=1000,prod_chan=prodchan,partlistfile=partlistfile,outfile=parfile,signal_chan="NCE_nucleon_baryonic",eps=0,alpha_D=1e-5,dm_energy_resolution=0.01,sumlog="Events/miniboone_baryonic_y.dat",output_mode="summary",zmin=zmin,zmax=zmax)
        write_miniboone(mdm=MX/1000.0,mv=MV/1000.0,proddist=proddist,\
        prod_chan=prodchan,partlistfile=partlistfile,\
        min_scatter_energy=0.05,max_scatter_energy=2,outfile=parfile,\
        samplesize=40000,eps=0,alpha_D=1e-5,signal_chan="NCE_nucleon_baryonic",\
        sumlog="Events_fit/summary_baryonic.dat",output_mode="comprehensive",\
        outlog="Events_fit/mini_baryonic_{}_{}.dat".format(MV,MX),\
        det=miniboone_detector_full,efficiency=1,zmin=zmin,zmax=zmax)
    subp.call(["./build/main", parfile])
    t1 = time.time()
    print("\ntime={}\n".format(t1-t0))
    t0 = time.time()
    #if MX/1000.0<meta/2.0:
    #    write_miniboone(mdm=MX/1000.0,mv=MV/1000.0,proddist="particle_list",prod_chan="eta_decay",outfile=parfile)
    #    subp.call(["./build/main", parfile])
    #    t1 = time.time()
    #    print "\ntime={}\n".format(t1-t0)
    subp.call(["rm", parfile])

def ship_eval(mass_arr,signal_channel="NCE_electron"):
    MV=mass_arr[0]
    MX=mass_arr[1]
    t0 = time.time()
    parfile="parameter_run_{0}_{1}.dat".format(str(MV),str(MX))
    #if MX/1000.0<mpi0/2.0:


    rho_decay_switch=False
    proton_brem_switch=True
    parton_switch = True;
    alD=0.1

    zmin=max(3*MV/1000.0/8.9,0.3)
    zmax=1-zmin

    proddist = []
    prodchan = []
    partlistfile = []
    if signal_channel=="Inelastic_Nucleon_Scattering_Baryonic":
        if MX/1000.0<mpi0/2.0 and MV<600:
            proddist.append("particle_list")
            prodchan.append("pi0_decay_baryonic")
            partlistfile.append("data/particle_list_ship.dat")
        if MX/1000.0<meta/2.0 and MV<900:
            proddist.append("particle_list")
            prodchan.append("eta_decay_baryonic")
            partlistfile.append("data/particle_list_ship.dat")
        #if ((MV<1200) and (MV>=350)) and rho_decay_switch:
        #    proddist.append("particle_list")
        #    prodchan.append("omega_decay_baryonic")
        #    partlistfile.append("data/particle_list.dat")
        if MV/1000.0>0.9 and MV<1250:
            proddist.append("particle_list")
            prodchan.append("phi_decay_baryonic")
            partlistfile.append("data/particle_list_ship.dat")
        if MV/1000.0<=4:
            proddist.append("proton_brem_baryonic")
            prodchan.append("V_decay_baryonic")
            partlistfile.append("")
        if MV/1000.0>=mrho:
           proddist.append("parton_V_baryonic")
           prodchan.append("parton_production_baryonic")
           partlistfile.append("")
    else:
        if MX/1000.0<mpi0/2.0 and MV<600.0:
            proddist.append("particle_list")
            prodchan.append("pi0_decay")
            partlistfile.append("data/particle_list_ship.dat")
        if MX/1000.0<meta/2.0 and MV<900.0:
            proddist.append("particle_list")
            prodchan.append("eta_decay")
            partlistfile.append("data/particle_list_ship.dat")
        if MV/1000.0>=mrho and parton_switch:
    	   proddist.append("parton_V")
    	   prodchan.append("parton_production")
    	   partlistfile.append("")
        if MV/2.0>MX and proton_brem_switch:
    	   proddist.append("proton_brem")
    	   prodchan.append("V_decay")
     	   partlistfile.append("")
        if ((MV<1200) and (MV>=350)) and rho_decay_switch:
           proddist.append("particle_list")
           prodchan.append("rho_decay")
           partlistfile.append("data/particle_list_ship.dat")
           proddist.append("particle_list")
           prodchan.append("omega_decay")
           partlistfile.append("data/particle_list_ship.dat")
        if MV/1000.0>=mrho and MV<=1250:
            proddist.append("particle_list")
            prodchan.append("phi_decay")
            partlistfile.append("data/particle_list_ship.dat")

    #write_ship(mdm=MX/1000.0,mv=MV/1000.0,proddist=proddist,prod_chan=prodchan,partlistfile=partlistfile,outfile=parfile,signal_chan="NCE_electron",samplesize=2000,sumlog="Events/ship.dat")
    if signal_channel =="NCE_electron":
        write_ship(mdm=MX/1000.0,mv=MV/1000.0,proddist=proddist,prod_chan=prodchan,partlistfile=partlistfile,outfile=parfile,signal_chan="NCE_electron",samplesize=1000,sumlog="Events/ship_y.dat",output_mode="summary",alpha_D=alD)
        #write_ship(mdm=MX/1000.0,mv=MV/1000.0,proddist=proddist,prod_chan=prodchan,partlistfile=partlistfile,outfile=parfile,signal_chan="NCE_electron",samplesize=1e6,sumlog="Events/ship_test.dat",output_mode="comprehensive",outlog="Events/walter_events.dat",alpha_D=alD,det=test_detector,min_scatter_angle=0,max_scatter_angle=10,min_scatter_energy=0,max_scatter_energy=400)
    elif signal_channel =="Inelastic_Nucleon_Scattering_Baryonic":
        write_ship(mdm=MX/1000.0,mv=MV/1000.0,proddist=proddist,prod_chan=prodchan,partlistfile=partlistfile,outfile=parfile,signal_chan="Inelastic_Nucleon_Scattering_Baryonic",samplesize=2000,sumlog="Events/ship_y.dat",output_mode="summary",alpha_D=1e-4,eps=0.0)
    elif signal_channel == "Inelastic_Nucleon_Scattering":
        write_ship(mdm=MX/1000.0,mv=MV/1000.0,proddist=proddist,prod_chan=prodchan,partlistfile=partlistfile,outfile=parfile,signal_chan="Inelastic_Nucleon_Scattering",samplesize=2000,sumlog="Events/ship_y.dat",output_mode="summary",alpha_D=alD)
    elif signal_channel=="test":
        write_ship(mdm=MX/1000.0,mv=MV/1000.0,proddist=proddist,prod_chan=prodchan,partlistfile=partlistfile,outfile=parfile,signal_chan="NCE_electron",samplesize=1000000,sumlog="Events/ship_test.dat",outlog="dm_dist_{}_{}.dat".format(MV,MX),output_mode="dm_detector_distribution",alpha_D=alD,det=test_detector,min_scatter_angle=0,max_scatter_angle=6,min_scatter_energy=0,max_scatter_energy=400)
    subp.call(["./build/main", parfile])
    t1 = time.time()
    print("\ntime={}\n".format(t1-t0))
    t0 = time.time()
    #if MX/1000.0<meta/2.0:
    #    write_miniboone(mdm=MX/1000.0,mv=MV/1000.0,proddist="particle_list",prod_chan="eta_decay",outfile=parfile)
    #    subp.call(["./build/main", parfile])
    #    t1 = time.time()
    #    print "\ntime={}\n".format(t1-t0)
    subp.call(["rm", parfile])


def miniboone_eval(mass_arr,channels={_parton,_brem,_pion_decay,_eta_decay},signal_channel = "NCE_nucleon",det_switch="miniboone",alpha_D=0.1,sumlog="Events/miniboone_testing.dat"):

    MV=mass_arr[0]
    MX=mass_arr[1]
    t0 = time.time()
    parfile="parameter_run_{0}_{1}.dat".format(str(MV),str(MX))
    #if MX/1000.0<mpi0/2.0:

    zmin=max(3*MV/1000.0/8.9,0.3)
    zmax=1-zmin

    proddist = []
    prodchan = []
    partlistfile = [    ]
    executing=False
    if MX/1000.0<mpi0/2.0 and MV<600.0 and _pion_decay in channels:
        proddist.append("particle_list")
        prodchan.append("pi0_decay")
        partlistfile.append("data/particle_list.dat")
        executing=True
    if MX/1000.0<meta/2.0 and MV<900.0 and _eta_decay in channels:
        proddist.append("particle_list")
        prodchan.append("eta_decay")
        partlistfile.append("data/particle_list_k0.dat")
        executing=True
    if MV/1000.0>=mrho and _parton in channels:
	proddist.append("parton_V")
	prodchan.append("parton_production")
	partlistfile.append("")
        executing=True
    if MV/2.0>MX and zmin<zmax and _brem in channels:
        proddist.append("proton_brem")
        prodchan.append("V_decay")
        partlistfile.append("")
        executing=True
    if ((MV<1200) and (MV>=350)) and _rho_decay in channels:
        proddist.append("particle_list")
        prodchan.append("rho_decay")
        partlistfile.append("data/particle_list.dat")
        proddist.append("particle_list")
        prodchan.append("omega_decay")
        partlistfile.append("data/particle_list.dat")
        executing=True
   # if MV/1000.0>=mrho and MV<=1250:
   #     proddist.append("particle_list")
   #     prodchan.append("phi_decay")
   #     partlistfile.append("data/particle_list.dat")
    if not executing:
        return
    if signal_channel=="NCE_nucleon":
        if det_switch == "sbnd":
            print("SBND RUN IN PROGRESS")
            write_sbnd(mdm=MX/1000.0,mv=MV/1000.0,proddist=proddist,prod_chan=prodchan,partlistfile=partlistfile,outfile=parfile,samplesize=2000,signal_chan="NCE_nucleon",sumlog="Events/sbnd_y.dat",zmin=zmin,zmax=zmax,alpha_D=alpha_D)
        else:
            write_miniboone(mdm=MX/1000.0,mv=MV/1000.0,proddist=proddist,prod_chan=prodchan,partlistfile=partlistfile,outfile=parfile,samplesize=2000,signal_chan="NCE_nucleon",sumlog=sumlog,zmin=zmin,zmax=zmax,min_scatter_energy=0,max_scatter_energy=8,alpha_D=alpha_D,eps=1e-3,efficiency=1,det=miniboone_detector,output_mode="summary",outlog="Events/miniboone_events.dat")
            #write_miniboone(mdm=MX/1000.0,mv=MV/1000.0,proddist=proddist,prod_chan=prodchan,partlistfile=partlistfile,outfile=parfile,samplesize=2000,signal_chan="NCE_nucleon",sumlog=sumlog,zmin=zmin,zmax=zmax,min_scatter_energy=0.035,max_scatter_energy=2,alpha_D=alpha_D,eps=1e-3,efficiency=0.35,det=miniboone_detector,output_mode="summary",outlog="Events/miniboone_events.dat")
            #write_miniboone(mdm=MX/1000.0,mv=MV/1000.0,proddist=proddist,prod_chan=prodchan,partlistfile=partlistfile,min_scatter_energy=0.05,max_scatter_energy=2,outfile=parfile,samplesize=40000,signal_chan="NCE_nucleon",sumlog="Events3/summary.dat",efficiency=1.0,output_mode="comprehensive",outlog="Events3/mini_{}_{}.dat".format(MV,MX),det=miniboone_detector_full)
    elif signal_channel=="Pion_Inelastic":
        if det_switch == "sbnd":
            print("SBND RUN IN PROGRESS")
            write_sbnd(mdm=MX/1000.0,mv=MV/1000.0,proddist=proddist,prod_chan=prodchan,samplesize=1000,min_scatter_energy=0,max_scatter_energy=1000.0,partlistfile=partlistfile,outfile=parfile,signal_chan="Pion_Inelastic",sumlog="Events/sbnd_pion_inelastic.dat")
        else:
            write_miniboone(mdm=MX/1000.0,mv=MV/1000.0,proddist=proddist,prod_chan=prodchan,samplesize=10000,min_scatter_energy=0,max_scatter_energy=10.0,partlistfile=partlistfile,outfile=parfile,signal_chan="Pion_Inelastic",sumlog="Events/miniboone_pion_test.dat",alpha_D=alpha_D)
    elif signal_channel=="Inelastic_Delta_to_Gamma":
        if det_switch == "sbnd":
            print("SBND RUN")
            write_sbnd(mdm=MX/1000.0,mv=MV/1000.0,proddist=proddist,prod_chan=prodchan,samplesize=1000,min_scatter_energy=0,max_scatter_energy=1000.0,partlistfile=partlistfile,outfile=parfile,signal_chan="Inelastic_Delta_to_Gamma",sumlog="Events/gamma_inelastic_sbnd.dat")
        else:
            write_miniboone(mdm=MX/1000.0,mv=MV/1000.0,proddist=proddist,prod_chan=prodchan,samplesize=1000,min_scatter_energy=0,max_scatter_energy=10.0,partlistfile=partlistfile,outfile=parfile,signal_chan="Inelastic_Delta_to_Gamma",sumlog="Events/miniboone_gamma_inelastic.dat",alpha_D=alpha_D,output_mode="comprehensive")
    elif signal_channel=="NCE_electron":
        if det_switch == "sbnd":
            write_sbnd(mdm=MX/1000.0,mv=MV/1000.0,proddist=proddist,prod_chan=prodchan,samplesize=1000,min_scatter_energy=0,max_scatter_energy=1.0,max_scatter_angle=0.14,partlistfile=partlistfile,outfile=parfile,signal_chan="NCE_electron",sumlog="Events/sbnd_elec_y.dat",zmin=zmin,zmax=zmax,alpha_D=alpha_D,output_mode="summary")
        else:
            write_miniboone(mdm=MX/1000.0,mv=MV/1000.0,proddist=proddist,prod_chan=prodchan,samplesize=1000,min_scatter_energy=0,max_scatter_energy=1.0,max_scatter_angle=0.14,partlistfile=partlistfile,outfile=parfile,signal_chan="NCE_electron",sumlog="Events/miniboone_electron_10mev.dat",zmin=zmin,zmax=zmax,alpha_D=alpha_D,output_mode="summary")
    subp.call(["./build/main", parfile])
    t1 = time.time()
    print("\ntime={}\n".format(t1-t0))
    t0 = time.time()
    #if MX/1000.0<meta/2.0:
    #    write_miniboone(mdm=MX/1000.0,mv=MV/1000.0,proddist="particle_list",prod_chan="eta_decay",outfile=parfile)
    #    subp.call(["./build/main", parfile])
    #    t1 = time.time()
    #    print "\ntime={}\n".format(t1-t0)
    subp.call(["rm", parfile])

def miniboone_numi_eval(mass_arr,channels={_parton,_brem,_pion_decay,_eta_decay},signal_channel = "NCE_nucleon",det_switch="miniboone_numi",alpha_D=0.5,sumlog="Events/miniboone_numi_500mev_15gev_y.dat"):

    MV=mass_arr[0]
    MX=mass_arr[1]
    t0 = time.time()
    parfile="parameter_run_{0}_{1}.dat".format(str(MV),str(MX))
    #if MX/1000.0<mpi0/2.0:

    BEAM_ENERGY=120
    alD=alpha_D
    zmin=max(3*MV/1000.0/BEAM_ENERGY,0.1)
    zmax=1-zmin

    proddist = []
    prodchan = []
    partlistfile = [    ]
    executing=False
    if MX/1000.0<mpi0/2.0 and MV<600.0 and _pion_decay in channels:

        proddist.append("particle_list")
        prodchan.append("pi0_decay")
        partlistfile.append("data/particle_list.dat")
        executing=True
    if MX/1000.0<meta/2.0 and MV<900.0 and _eta_decay in channels:
        proddist.append("particle_list")
        prodchan.append("eta_decay")
        partlistfile.append("data/particle_list_k0.dat")
        executing=True
    if MV/1000.0>=mrho and _parton in channels:
	proddist.append("parton_V")
	prodchan.append("parton_production")
	partlistfile.append("")
        executing=True
    if MV/2.0>MX and zmin<zmax and _brem in channels:
        proddist.append("proton_brem")
        prodchan.append("V_decay")
        partlistfile.append("")
        executing=True
    if ((MV<1200) and (MV>=350)) and _rho_decay in channels:
        proddist.append("particle_list")
        prodchan.append("rho_decay")
        partlistfile.append("data/particle_list.dat")
        proddist.append("particle_list")
        prodchan.append("omega_decay")
        partlistfile.append("data/particle_list.dat")
        executing=True
   # if MV/1000.0>=mrho and MV<=1250:
   #     proddist.append("particle_list")
   #     prodchan.append("phi_decay")
   #     partlistfile.append("data/particle_list.dat")
    if not executing:
        return
    if signal_channel=="NCE_nucleon":
        write_miniboone_numi(mdm=MX/1000.0,mv=MV/1000.0,proddist=proddist,prod_chan=prodchan,partlistfile=partlistfile,outfile=parfile,samplesize=2000,signal_chan="NCE_nucleon",sumlog=sumlog,zmin=zmin,zmax=zmax,min_scatter_energy=0.5,max_scatter_energy=15,alpha_D=alpha_D,eps=1e-3,efficiency=0.35,det=miniboone_detector,output_mode="summary",outlog="Events/miniboone_events.dat")
    if signal_channel=="Pion_Inelastic":
        write_miniboone_numi(mdm=MX/1000.0,mv=MV/1000.0,proddist=proddist,prod_chan=prodchan,samplesize=2000,min_scatter_energy=0,max_scatter_energy=10.0,partlistfile=partlistfile,outfile=parfile,signal_chan="Pion_Inelastic",sumlog="Events/miniboone_pion_numi.dat",alpha_D=alpha_D)
    if signal_channel=="NCE_electron":
        write_miniboone_numi(mdm=MX/1000.0,mv=MV/1000.0,proddist=proddist,prod_chan=prodchan,samplesize=1000,min_scatter_energy=0.5,max_scatter_energy=3,max_scatter_angle=0.14,partlistfile=partlistfile,outfile=parfile,signal_chan="NCE_electron",sumlog="Events/miniboone_electron_numi_500mev_3gev.dat",zmin=zmin,zmax=zmax,alpha_D=alpha_D,output_mode="summary")
    if signal_channel=="Inelastic_Nucleon_Scattering":
        write_miniboone_numi(mdm=MX/1000.0,mv=MV/1000.0,proddist=proddist,prod_chan=prodchan,partlistfile=partlistfile,outfile=parfile,signal_chan="Inelastic_Nucleon_Scattering",samplesize=2000,sumlog="Events/miniboone_dis.dat",output_mode="summary",alpha_D=alpha_D)
            #write_miniboone(mdm=MX/1000.0,mv=MV/1000.0,proddist=proddist,prod_chan=prodchan,partlistfile=partlistfile,min_scatter_energy=0.05,max_scatter_energy=2,outfile=parfile,samplesize=40000,signal_chan="NCE_nucleon",sumlog="Events3/summary.dat",efficiency=1.0,output_mode="comprehensive",outlog="Events3/mini_{}_{}.dat".format(MV,MX),det=miniboone_detector_full)
    subp.call(["./build/main", parfile])
    t1 = time.time()
    print("\ntime={}\n".format(t1-t0))
    t0 = time.time()
    #if MX/1000.0<meta/2.0:
    #    write_miniboone(mdm=MX/1000.0,mv=MV/1000.0,proddist="particle_list",prod_chan="eta_decay",outfile=parfile)
    #    subp.call(["./build/main", parfile])
    #    t1 = time.time()
    #    print "\ntime={}\n".format(t1-t0)
    subp.call(["rm", parfile])

def t2k_FD_baryonic_eval(mass_arr,signal_channel="NCE_nucleon_baryonic"):
    MV=mass_arr[0]
    MX=mass_arr[1]
    t0 = time.time()
    parfile="parameter_run_{0}_{1}.dat".format(str(MV),str(MX))
    #if MX/1000.0<mpi0/2.0:

    rho_decay_switch=False

    zmin=max(3*MV/1000.0/30.0,0.2)
    zmax=1-zmin

    print("Zmin=",zmin,"Zmax=",zmax)

    proddist = []
    prodchan = []
    partlistfile = []
    if MX/1000.0<mpi0/2.0 and MV<600:
        proddist.append("particle_list")
        prodchan.append("pi0_decay_baryonic")
        partlistfile.append("data/particle_list_t2k.dat")
    if MX/1000.0<meta/2.0 and MV<900:
        proddist.append("particle_list")
        prodchan.append("eta_decay_baryonic")
        partlistfile.append("data/particle_list_t2k.dat")
    #if ((MV<1200) and (MV>=350)) and rho_decay_switch:
    #    proddist.append("particle_list")
    #    prodchan.append("omega_decay_baryonic")
    #    partlistfile.append("data/particle_list.dat")
    if MV/1000.0>0.9 and MV<1250:
        proddist.append("particle_list")
        prodchan.append("phi_decay_baryonic")
        partlistfile.append("data/particle_list_t2k.dat")
    if MV/1000.0<=10:
        proddist.append("proton_brem_baryonic")
        prodchan.append("V_decay_baryonic")
        partlistfile.append("")
    if MV/1000.0>=mrho:
	   proddist.append("parton_V_baryonic")
	   prodchan.append("parton_production_baryonic")
	   partlistfile.append("")
    if signal_channel=="NCE_nucleon_baryonic":
        write_t2kFD(mdm=MX/1000.0,mv=MV/1000.0,samplesize=500,proddist=proddist,prod_chan=prodchan,alpha_D=1e-4,eps=0.0,zmin=zmin,zmax=zmax,signal_chan="NCE_nucleon_baryonic",min_scatter_energy=0.0014,max_scatter_energy=1000.0,partlistfile=partlistfile,outfile=parfile,sumlog="Events/superk_y.dat")
    subp.call(["./build/main", parfile])
    t1 = time.time()
    print("\ntime={}\n".format(t1-t0))
    t0 = time.time()
    #if MX/1000.0<meta/2.0:
    #    write_miniboone(mdm=MX/1000.0,mv=MV/1000.0,proddist="particle_list",prod_chan="eta_decay",outfile=parfile)
    #    subp.call(["./build/main", parfile])
    #    t1 = time.time()
    #    print "\ntime={}\n".format(t1-t0)
    subp.call(["rm", parfile])

def t2k_FD_eval(mass_arr,signal_channel="NCE_nucleon"):
    MV=mass_arr[0]
    MX=mass_arr[1]
    t0 = time.time()
    parfile="parameter_run_{0}_{1}.dat".format(str(MV),str(MX))
    #if MX/1000.0<mpi0/2.0:
    alD=0.5
    rho_decay_switch=False

    zmin=max(3*MV/1000.0/30.0,0.2)
    zmax=1-zmin

    print("Zmin=",zmin,"Zmax=",zmax)

    proddist = []
    prodchan = []
    partlistfile = []
    if MX/1000.0<mpi0/2.0 and MV<600:
        proddist.append("particle_list")
        prodchan.append("pi0_decay")
        partlistfile.append("data/particle_list_t2k.dat")
    if MX/1000.0<meta/2.0 and MV<600:
        proddist.append("particle_list")
        prodchan.append("eta_decay")
        partlistfile.append("data/particle_list_t2k.dat")
    if MV/1000.0>0.9 and MV<1250:
        proddist.append("particle_list")
        prodchan.append("phi_decay")
        partlistfile.append("data/particle_list_t2k.dat")
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
        partlistfile.append("data/particle_list_k0.dat")
    if signal_channel=="NCE_nucleon":
        write_t2kFD(mdm=MX/1000.0,mv=MV/1000.0,zmin=zmin,zmax=zmax,samplesize=500,proddist=proddist,prod_chan=prodchan,signal_chan="NCE_nucleon",min_scatter_energy=0.0014,max_scatter_energy=1000.0,partlistfile=partlistfile,outfile=parfile,sumlog="Events/superk_y.dat",alpha_D=alD)
    elif signal_channel == "Pion_Inelastic":
        write_t2kFD(mdm=MX/1000.0,mv=MV/1000.0,samplesize=1000,proddist=proddist,signal_chan="Pion_Inelastic",min_scatter_energy=0.0,max_scatter_energy=1000.0,prod_chan=prodchan,partlistfile=partlistfile,outfile=parfile)
    elif signal_channel == "NCE_electron":
        write_t2kFD(mdm=MX/1000.0,mv=MV/1000.0,samplesize=1000,proddist=proddist,signal_chan="NCE_electron",min_scatter_energy=0.0,max_scatter_energy=1000.0,prod_chan=prodchan,partlistfile=partlistfile,outfile=parfile)
    subp.call(["./build/main", parfile])
    t1 = time.time()
    print("\ntime={}\n".format(t1-t0))
    t0 = time.time()
    #if MX/1000.0<meta/2.0:
    #    write_miniboone(mdm=MX/1000.0,mv=MV/1000.0,proddist="particle_list",prod_chan="eta_decay",outfile=parfile)
    #    subp.call(["./build/main", parfile])
    #    t1 = time.time()
    #    print "\ntime={}\n".format(t1-t0)
    subp.call(["rm", parfile])

def t2k_eval(mass_arr,signal_channel="Pion_Inelastic"):
    MV=mass_arr[0]
    MX=mass_arr[1]
    t0 = time.time()
    parfile="parameter_run_{0}_{1}.dat".format(str(MV),str(MX))
    #if MX/1000.0<mpi0/2.0:

    zmin=max(3*MV/1000.0/30.0,0.2)
    zmax=1-zmin

    print("Zmin=",zmin,"Zmax=",zmax)

    rho_decay_switch=False
    alD=0.5
    proddist = []
    prodchan = []
    partlistfile = []
    if MX/1000.0<mpi0/2.0:
        proddist.append("particle_list")
        prodchan.append("pi0_decay")
        partlistfile.append("data/particle_list_t2k.dat")
    if MX/1000.0<meta/2.0:
        proddist.append("particle_list")
        prodchan.append("eta_decay")
        partlistfile.append("data/particle_list_t2k.dat")
    if MV/1000.0>=mrho:
        proddist.append("parton_V")
        prodchan.append("parton_production")
        partlistfile.append("")
    if MV/1000.0>0.9 and MV<1250:
        proddist.append("particle_list")
        prodchan.append("phi_decay")
        partlistfile.append("data/particle_list_t2k.dat")
    if MV/1000.0<=2:
	    proddist.append("proton_brem")
	    prodchan.append("V_decay")
 	    partlistfile.append("")
    if ((MV<1200) and (MV>500)) and rho_decay_switch:
        proddist.append("particle_list")
        prodchan.append("rho_decay")
        partlistfile.append("data/particle_list_t2k.dat")
    if signal_channel=="NCE_nucleon":
        write_nd280(mdm=MX/1000.0,mv=MV/1000.0,proddist=proddist,prod_chan=prodchan,partlistfile=partlistfile,outfile=parfile)
    elif signal_channel == "Pion_Inelastic":
        write_nd280(mdm=MX/1000.0,mv=MV/1000.0,samplesize=1000,zmin=zmin,zmax=zmax,proddist=proddist,signal_chan="Pion_Inelastic",min_scatter_energy=0.0,max_scatter_energy=1000.0,prod_chan=prodchan,partlistfile=partlistfile,outfile=parfile,sumlog="Events/t2k_pod_y.dat",alpha_D=0.5)
    subp.call(["./build/main", parfile])
    t1 = time.time()
    print("\ntime={}\n".format(t1-t0))
    t0 = time.time()
    #if MX/1000.0<meta/2.0:
    #    write_miniboone(mdm=MX/1000.0,mv=MV/1000.0,proddist="particle_list",prod_chan="eta_decay",outfile=parfile)
    #    subp.call(["./build/main", parfile])
    #    t1 = time.time()
    #    print "\ntime={}\n".format(t1-t0)
    subp.call(["rm", parfile])

def lsnd_eval(mass_arr,signal_channel="NCE_electron"):
    MV=mass_arr[0]
    MX=mass_arr[1]
    t0 = time.time()
    parfile="parameter_run_{0}_{1}.dat".format(str(MV),str(MX))
    #if MX/1000.0<mpi0/2.0:
    alD=0.1
    rho_decay_switch=False
    proddist = []
    prodchan = []
    partlistfile = []
    if MX/1000.0<mpi0/2.0:
        proddist.append("particle_list")
        prodchan.append("pi0_decay")
        partlistfile.append("data/particle_list_lsnd.dat")
    else:
        return
    if signal_channel=="NCE_electron":
        #write_lsnd(mdm=MX/1000.0,mv=MV/1000.0,proddist=proddist,prod_chan=prodchan,partlistfile=partlistfile,outfile=parfile,sumlog="Events/lsnd_interp.dat",p_num_target=8,samplesize=1000,alpha_D=alD)
        write_lsnd(mdm=MX/1000.0,mv=MV/1000.0,proddist=proddist,prod_chan=prodchan,partlistfile=partlistfile,outfile=parfile,sumlog="Events/lsnd_10mev.dat",p_num_target=8,samplesize=1000,alpha_D=alD)
    elif signal_channel=="NCE_nucleon":
        write_lsnd(mdm=MX/1000.0,mv=MV/1000.0,proddist=proddist,prod_chan=prodchan,partlistfile=partlistfile,outfile=parfile,min_scatter_energy=0.018, max_scatter_energy=1000,sumlog="Events/lsnd_thesis3.dat",p_num_target=8,samplesize=1000,signal_chan="NCE_nucleon")
    subp.call(["./build/main", parfile])
    t1 = time.time()
    print("\ntime={}\n".format(t1-t0))
    t0 = time.time()
    subp.call(["rm", parfile])

def execute_miniboone_parallel(genlist=True):
    if genlist:
        write_miniboone(prod_chan=["pi0_decay"],proddist=["pi0_sanfordwang"],samplesize=2e6,output_mode="particle_list",partlistfile=["data/particle_list.dat"])
        subp.call(["./build/main", "parameter_run.dat"])
        write_miniboone(prod_chan=["eta_decay"],proddist=["k0_sanfordwang"],samplesize=2e6,output_mode="particle_list",partlistfile=["data/particle_list_k0.dat"])
        subp.call(["./build/main", "parameter_run.dat"])
        #write_miniboone(prod_chan=["V_decay"],proddist=["proton_brem"],samplesize=4e6,output_mode="particle_list",partlistfile=["data/particle_list_miniboone_brem.dat"])
    #chimassarr=[10]
    #vmassarr=[1,5]+[10*i for i in range(1,14)]+[10*i for i in range(15,100,2)]+[775,774,776,777,778,779,781,782,783,785,787]+[1005,1010,1015,1.02,1025,1030,1035,1050,1100,1150,1200,1300]+[i for i in range(1000,2100,250)]
    #vmassarr=[19,21,17,15,23]
    #vmassarr=[i for i in range(10,140,10)]+[i for i in range(150,1000,50)]+[770,772,768,762,778]+[1000]+[3,5,7,9]
    #chimassarr=[i for i in range(10,270,10)]+[i for i in range(270,560,20)]+[132,134,136]+[1,5]
    #This is for the 2016 paper.
    chimassarr=[i for i in range(10,270,5)]+[i for i in range(270,501,10)]+[1,2,3,4,5,6,7,8,9,11,12,13,14]+[770/3.0,772/3.0,768/3.0,762/3.0,778/3.0]
    #chimassarr=[i for i in range(525,700,25)]
    massarr=[[3*MX,MX] for MX in chimassarr]
    #chimassarr=[5*i for i in range(13)]+[10*i+65 for i in range(20)]
    #vmassarr=[10*i for i in range(13)]+[20*i+130 for i in range(21)]+[540,135]
    #massarr=np.array([[MV,MX] for MV in vmassarr for MX in chimassarr])
    #massarr=massarr[massarr[:,0]>2*massarr[:,1]]
    #chimassarr=[i for i in range(10,270,10)]+[132,134,136]+[1,5]
    #massarr=[[MV,MV/6.0] for MV in vmassarr]
    #massarr=[[100,11]]
    #for marr in massarr:
    #chimassarr=[7,8,9,10,11,12,13]
    #    miniboone_baryonic_eval(marr,det_switch="miniboone")
    #massarr=[[MV,MX] for MV in vmassarr for MX in chimassarr]
    #massarr2=[]
    #for i in massarr:
    #    if i not in massarr2:
    #        massarr2.append(i)
    #massarr=massarr2
    #massarr=[[80,10],[150,10],[300,40]]
    print(massarr)
    #channs={_parton,_brem,_pion_decay,_eta_decay}
    #pool=Pool(processes=3)
    #pool.map(miniboone_eval,massarr)
    for marr in massarr:
        #miniboone_baryonic_eval(marr,det_switch="miniboone")
        #for chan in channs:
        #    miniboone_eval(marr,signal_channel="NCE_nucleon",det_switch="miniboone",channels=[chan],sumlog="Events/miniboone_split.dat")
        miniboone_eval(marr,signal_channel="Pion_Inelastic",det_switch="sbnd")
        #miniboone_eval(marr,signal_channel="Inelastic_Delta_to_Gamma",det_switch="miniboone")
        miniboone_eval(marr,signal_channel="NCE_electron",det_switch="sbnd")
        #miniboone_eval(marr,signal_channel="NCE_nucleon",det_switch="miniboone")
        #miniboone_baryonic_eval(marr,det_switch="miniboone")
        #if len(sys.argv)>1 and sys.argv[1]=="b":

        #else:
            #miniboone_eval(marr,signal_channel="NCE_nucleon",det_switch="miniboone")
    #print massarr
    #pool=Pool(processes=3)
    #ool.map(miniboone_eval,massarr)


def execute_miniboone_numi_p(genlist=True):
    if genlist:
        write_miniboone_numi(prod_chan=["pi0_decay"],proddist=["bmpt"],samplesize=1e6,output_mode="particle_list",partlistfile=["data/particle_list_numi.dat"])
        subp.call(["./build/main", "parameter_run.dat"])
    vmassarr=[i for i in range(10,140,10)]+[i for i in range(150,1000,50)]+[770,772,768,762,778]+[1000,1100,1150,1200,1300]+[i for i in range(1000,5100,250)]+[1005,1010,1020,1015,1025,1030]+[3,5,7,9]+[30,210,420,600,810,1020,1200,1500,1710]
    massarr=[[MV,MV/3.0] for MV in vmassarr]
    for marr in massarr:
        miniboone_numi_eval(marr,signal_channel="NCE_nucleon")
        miniboone_numi_eval(marr,signal_channel="NCE_electron")
        #iminiboone_numi_eval(marr,signal_channel="Inelastic_Nucleon_Scattering")
        #miniboone_numi_eval(marr,signal_channel="Pion_Inelastic")

def execute_ship_parallel(genlist=True):
    if genlist:
        write_ship(prod_chan=["pi0_decay"],proddist=["bmpt"],samplesize=1e6,output_mode="particle_list",partlistfile=["data/particle_list_ship.dat"])
        #write_ship(prod_chan=["parton_production"],proddist=["parton_V"],samplesize=2e5,output_mode="particle_list",partlistfile=["data/particle_list_ship_parton.dat"])
        #write_ship(prod_chan=["V_decay"],proddist=["proton_brem"],samplesize=4e6,output_mode="particle_list",partlistfile=["data/particle_list_ship_brem.dat"])
        subp.call(["./build/main", "parameter_run.dat"])
    #vmassarr=[i for i in range(100,1000,50)]+[775,774,776,777,778,779,781,782,783,785,787]+[1005,1010,1015,1.02,1025,1030,1035,1050,1100,1150,1200,1300]+[i for i in range(1000,4100,250)]+[390,395,405,410]
    #chimassarr=[1,5]+[10*i for i in range(1,27)]
    #vmassarr=[600,700,800,900]
    #vmassarr=[401,405,425,450,475,772,768,762,778]
    #vmassarr=vmassarr+[771,769,773,767]
    #vmassarr=vmassarr+[400,500,600,700,760,765,770,775,780,800,900,1000,1010,1015,1020,1025,1030,1250,1500,1750,2000,2250,2500,2750,3000,3250,3500,3750,4000]
    #vmassarr=[850,750,825,875,725]
    #vmassarr=[i for i in range(2000,3200,200)]
    vmassarr=[i for i in range(10,140,10)]+[i for i in range(150,1000,50)]+[770,772,768,762,778]+[1000,1100,1150,1200,1300]+[i for i in range(1000,5100,250)]+[1005,1010,1020,1015,1025,1030]+[3,5,7,9]+[30,210,420,600,810,1020,1200,1500,1710]
    #massarr=[[MV,MV/5.0] for MV in vmassarr]
    massarr=[[MV,MV/3.0] for MV in vmassarr]
    #massarr=[[MV,200] for MV in vmassarr]
    for marr in massarr:
        ship_eval(marr,signal_channel="NCE_electron")
        #ship_eval(marr,signal_channel="test")
        ship_eval(marr,signal_channel="Inelastic_Nucleon_Scattering")
        ship_eval(marr,signal_channel="Inelastic_Nucleon_Scattering_Baryonic")
    #pool=Pool(processes=4)
    #pool.map(ship_eval,massarr)

def execute_lsnd_parallel(genlist=True):
    if genlist:
        write_lsnd(prod_chan=["pi0_decay"],proddist=["burmansmith"],samplesize=0.55e6,output_mode="particle_list",partlistfile=["data/particle_list_lsnd_low.dat"],p_num_target=8)
        subp.call(["./build/main","parameter_run.dat"])
        write_lsnd(prod_chan=["pi0_decay"],proddist=["burmansmith"],samplesize=0.45e6,output_mode="particle_list",partlistfile=["data/particle_list_lsnd_high.dat"],p_num_target=70)
        subp.call(["./build/main","parameter_run.dat"])
        arr1 = np.loadtxt("data/particle_list_lsnd_low.dat")
        arr2 = np.loadtxt("data/particle_list_lsnd_high.dat")
        arr3 = np.append(arr1,arr2,axis=0)
        np.random.shuffle(arr3)
        np.savetxt("data/particle_list_lsnd.dat",arr3)
    #massarr=[]
    #chimassarr=[1,2,3,4,5,6,7,8,9,10,11,12,62]+[i for i in range(13,65,3)]
    chimassarr=[10]
    vmassarr=[i for i in range(10,30,1)]+[i for i in range(30,130,5)]+[131,132,133,134,135,136,137,138]+[i for i in range(140,700,20)]
    massarr=[[MV,MX] for MV in vmassarr for MX in chimassarr]
    #vmassarr=[i for i in range(5,400,5)]
    #massarr=[[MV,MV/5] for MV in vmassarr]
    #vmassarr=[i for i in range(600,910,50)]
    #chimassarr=[10]
    #vmassarr=[20,50,100]
    #chimassarr=[i for i in range(15,100,10)]+[11,12,13,14,16,18,23,27]
    #chimassarr=[5,7,9,49,51,24,26]
    #massarr=massarr+[[MV,MX] for MV in vmassarr for MX in chimassarr]
    #massarr=[[MV,MX] for MV in vmassarr for MX in chimassarr if MX<MV/2.0]+[[MV,MV/2.0] for MV in range(10,134,1)]
    #massarr = [[3*MX,MX] for MX in chimassarr]
    #for i in massarr:
    #    if i not in massarr2:
    #        massarr2.append(i)
    #massarr=massarr2
    #print(len(massarr))
    #print(massarr)
    for marr in massarr:
        lsnd_eval(marr,signal_channel="NCE_electron")
    #pool=Pool(processes=4)
    #pool.map(lsnd_eval,massarr)

def execute_t2k_parallel(genlist=True):
    if genlist:
        write_nd280(prod_chan=["pi0_decay"],proddist=["bmpt"],samplesize=1e6,output_mode="particle_list",partlistfile=["data/particle_list_t2k.dat"])
        #write_ship(prod_chan=["parton_production"],proddist=["parton_V"],samplesize=2e5,output_mode="particle_list",partlistfile=["data/particle_list_ship_parton.dat"])
        #write_ship(prod_chan=["V_decay"],proddist=["proton_brem"],samplesize=4e6,output_mode="particle_list",partlistfile=["data/particle_list_ship_brem.dat"])
        subp.call(["./build/main", "parameter_run.dat"])
    #chimassarr=[100]
    #vmassarr=[201,205,225,250,275,772,768,762,778]
    #vmassarr=vmassarr+[300,400,500,600,700,760,765,770,775,780,800,900,1010,1015,1020,1025,1030]+[i for i in range(1000,4100,250)]+[1100,1150,1200,1300]
    #chimassarr=[100]
    #vmassarr=[i for i in range(2250,4100,250)]
    #massarr=[[MV,MX] for MV in vmassarr for MX in chimassarr]
    #vmassarr=[400]
    #chimassarr=[10]
    #massarr=[[MV,MX] for MV in vmassarr for MX in chimassarr]
    vmassarr=[i for i in range(10,140,10)]+[i for i in range(150,1000,50)]+[770,772,768,762,778]+[1000,1100,1150,1200,1300]+[i for i in range(1000,4100,250)]+[1005,1010,1020,1015,1025,1030]+[3,5,7,9]
    massarr=[[MV,MV/3.0] for MV in vmassarr]
    #vmassarr=[400]
    #chimassarr=[10,20,30,40,50,100,150,175,180,185,190,195,199]
    #massarr=[[MV,MX] for MV in vmassarr for MX in chimassarr]
    print(massarr)
    for marr in massarr:
        t2k_eval(marr)
        t2k_FD_baryonic_eval(marr,signal_channel="NCE_nucleon_baryonic")
        t2k_FD_eval(marr,signal_channel="NCE_nucleon")
        #t2k_FD_eval(marr,signal_channel="Pion_Inelastic")
    #chimassarr=[100]
    #vmassarr=[900,1000,1010,1015,1020,1025,1030,1250,1500,1750,2000]
    #massarr=[[MV,MX] for MV in vmassarr for MX in chimassarr]
    #for marr in massarr:
        #t2k_FD_eval(marr,signal_channel="Pion_Inelastic")
    #pool=Pool(processes=4)
    #pool.map(ship_eval,massarr)

#execute_t2k_parallel(genlist=True)
execute_miniboone_parallel(genlist=True)
#execute_miniboone_numi_p(genlist=False)
#execute_ship_parallel(genlist=True)
#execute_lsnd_parallel(genlist=True)
