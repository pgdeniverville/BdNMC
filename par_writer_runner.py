from par_writer import *
import numpy as np
import sys
import math
import copy

import time

_pion_decay=1
_eta_decay=2
_rho_decay=3
_parton=4
_brem=5
_piminus_cap=6
_phi_decay=7

_pion_inelastic_channels={"Pion_Inelastic","Pion_Inelastic_Charged","Inelastic_Delta_to_Gamma"}

#write_miniboone(mdm=0.005,mv=0.4,proddist=["","","proton_brem"],prod_chan=["pi0_decay","eta_decay","V_decay"],partlistfile=["","",""])
rho_decay_switch=False

DET_XPOS = 0.0
_DET_YPOS = 0.0
_DET_ZPOS = 100

mass_delta = 1.1

def ship_detector_modular(f,radius=0.8268,length=3.34,theta=0,phi=0):
    #def ship_detector(f,xpos=0.0,ypos=0,zpos=30.0,radius=0.655,length=2.645,theta=0,phi=0):
    f.write("\ndetector cylinder\n");
    f.write("x-position {0}\ny-position {1}\nz-position {2}\nradius {3}\nlength {4}\ndet-theta {5}\ndet-phi {6}\n".format(str(_DET_XPOS),str(_DET_YPOS),str(_DET_ZPOS),str(radius),str(length),str(theta),str(phi)))
    f.write('\n')
    f.write(Argon_string)

def OPERA_detector(f,xpos=0.0,ypos=0.0,zpos=7.3e6,radius=56.4,length=1,theta=0,phi=0):
    print("Don't use this for event generation!")
    f.write("\ndetector cylinder\n");
    f.write("x-position {0}\ny-position {1}\nz-position {2}\nradius {3}\nlength {4}\ndet-theta {5}\ndet-phi {6}\n".format(str(xpos),str(ypos),str(zpos),str(radius),str(length),str(theta),str(phi)))
    f.write('\n')
    f.write(Argon_string)

def MINOS_detector_modular(f,xpos=0.0,radius=2,length=1.7,theta=0,phi=0):
    f.write("\ndetector cylinder\n");
    f.write("x-position {0}\ny-position {1}\nz-position {2}\nradius {3}\nlength {4}\ndet-theta {5}\ndet-phi {6}\n".format(str(_DET_XPOS),str(_DET_YPOS),str(_DET_ZPOS),str(radius),str(length),str(theta),str(phi)))
    f.write('\n')
    f.write(MINOS_string)
'''
def NOvA_detector_modular(f,xpos=0.0,radius=2,length=14,theta=-NOvA_angle,phi=0):
    f.write("\ndetector cylinder\n");
    f.write("x-position {0}\ny-position {1}\nz-position {2}\nradius {3}\nlength {4}\ndet-theta {5}\ndet-phi {6}\n".format(str(_DET_XPOS),str(_DET_YPOS),str(_DET_ZPOS),str(radius),str(length),str(theta),str(phi)))
    f.write('\n')
    f.write(MINOS_string)
'''
def miniboone_detector_modular(f,radius=5.0):
    f.write("\ndetector sphere\n")
    f.write("x-position {0}\ny-position {1}\nz-position {2}\nradius {3}\n".format(str(_DET_XPOS),str(_DET_YPOS),str(_DET_ZPOS),str(radius)))
    f.write('\n')
    f.write(Hydrogen_string)
    f.write('\n')
    f.write(Carbon_string)

pep2_defaults = {}

def pip2_eval(d_user):
    d=pep2_defaults.copy()
    d.update(d_user)
    MV=d["mv"]
    MX=d["mdm"]
    channels = d["channels"]
    det_switch=d["det_switch"]
    if "alpha_D" in d:
        alD=d["alpha_D"]
    if 'eps' in d:
        eps=d['eps']
    sumlog=d["sumlog"]
    signal_channel = d["signal_chan"]

    t0 = time.time()
    outfile="parameter_run_{0}_{1}.dat".format(str(MV),str(MX))

    proddist = []
    prodchan = []
    partlistfile = [    ]
    executing=False
    if MX/1000.0<mpi0/2.0 and _pion_decay in channels:
        proddist.append("particle_list")
        prodchan.append("pi0_decay")
        partlistfile.append("data/particle_list_pip2.dat")
        executing = True
    if MX/1000.0<0.129/2.0 and MV<600 and MV>2*MX and _piminus_cap in channels:
        proddist.append("")
        prodchan.append("piminus_capture")
        partlistfile.append("")
        executing = True
    if not executing:
        return
    d.update(d_user)
    d.update({"proddist" : proddist, "prod_chan" : prodchan, "partlistfile" : partlistfile,"mv" : MV/1000.0, "mdm" : MX/1000.0, "outfile" : outfile})
    if det_switch =="pip2":
        write_pip2(d=d)
    elif det_switch == "pip2_30":
        write_pip2(d=d,det=pip2_30)
    elif det_switch == "pip2_90":
        write_pip2(d=d,det=pip2_90)
    subp.call(["./build/main", outfile])
    t1 = time.time()
    print("\ntime={}\n".format(t1-t0))
    t0 = time.time()
    subp.call(["rm", outfile])

lanl_defaults = {"mv" : 30, "mdm" : 10, 'channels' : {_pion_decay,_piminus_cap}, "signal_chan" : "NCE_Nucleon", 'det_switch' : 'lanl', 'sumlog' : "Events/lanl.dat", "alpha_D" : 0.5, 'coherent' : True}

def lanl_eval(d_user):
    d=lanl_defaults.copy()
    d.update(d_user)
    MV=d["mv"]
    MX=d["mdm"]
    channels = d["channels"]
    det_switch=d["det_switch"]
    if "alpha_D" in d:
        alD=d["alpha_D"]
    if 'eps' in d:
        eps=d['eps']
    sumlog=d["sumlog"]
    signal_channel = d["signal_chan"]

    t0 = time.time()
    outfile="parameter_run_{0}_{1}.dat".format(str(MV),str(MX))

    proddist = []
    prodchan = []
    partlistfile = [    ]
    particle_list_position = []
    executing=False
    if signal_channel=="NCE_nucleon":
        if MX/1000.0<mpi0/2.0 and _pion_decay in channels:
            proddist.append("particle_list")
            prodchan.append("pi0_decay")
            #partlistfile.append("data/particle_list_lanl.dat")
            partlistfile.append("data/particle_list_lujan.dat")
            particle_list_position.append("true")
            executing = True
        if MX/1000.0<0.129/2.0 and MV<600 and MV>2*MX and _piminus_cap in channels:
            proddist.append("")
            prodchan.append("piminus_capture")
            partlistfile.append("")
            particle_list_position.append("false")
            executing = True
    elif signal_channel=="NCE_nucleon_baryonic":
        if MX/1000.0<mpi0/2.0 and MV<600.0 and _pion_decay in channels:
            proddist.append("particle_list")
            prodchan.append("pi0_decay_baryonic")
            #partlistfile.append("data/particle_list_lanl.dat")
            partlistfile.append("data/particle_list_lujan.dat")
            particle_list_position.append("true")
            executing = True
    elif signal_channel=="Signal_Decay":
        if MV/1000.0<mpi0 and _pion_decay in channels:
            proddist.append("particle_list")
            prodchan.append("pi0_decay")
            #partlistfile.append("data/particle_list_lanl.dat")
            partlistfile.append("data/particle_list_lujan.dat")
            particle_list_position.append("true")
            executing = True
    if not executing:
        return
    user2={}
    if signal_channel=="NCE_nucleon":
        if det_switch == "lanl":
            user2 = {"samplesize" : 1000, "POT" : 17.875e21, "min_scatter_energy" : 10e-6, "max_scatter_energy" : 0.2, "efficiency" : 0.8, "sumlog" : "Events/lanl.dat", "coherent" : "true", "eps" : 1e-3, "burn_max" : 1000}
        elif det_switch == "lanl_far":
            user2 = {"samplesize" : 1000, "POT" : 15.6e21, "min_scatter_energy" : 10e-6, "max_scatter_energy" : 0.2, "efficiency" : 0.8, "sumlog" : "Events/lanl_far.dat", "coherent" : "true", "eps" : 1e-3, "burn_max" : 1000}
    elif signal_channel=="NCE_nucleon_baryonic":
        if det_switch == "lanl":
            user2 = {"samplesize" : 1000, "POT" : 2.2e21, "min_scatter_energy" : 10e-6, "max_scatter_energy" : 0.2, "efficiency" : 0.8, "sumlog" : "Events/lanl_b.dat", "coherent" : "true", "eps" : 0, "burn_max" : 1000}
        elif det_switch == "lanl_far":
            user2 = {"samplesize" : 1000, "POT" : 15.6e21, "min_scatter_energy" : 10e-6, "max_scatter_energy" : 0.2, "efficiency" : 0.8, "sumlog" : "Events/lanl_far_b.dat", "coherent" : "true", "eps" : 0, "burn_max" : 1000}
    d.update(user2)
    d.update(d_user)
    d.update({"proddist" : proddist, "prod_chan" : prodchan, "partlistfile" : partlistfile,"mv" : MV/1000.0, "mdm" : MX/1000.0, "outfile" : outfile})
    if det_switch =="lanl":
        write_lanl(d=d,det=lanl_detector)
    #elif det_switch =="lanl_far":
    #    write_lanl(d=d,det=lanl_detector_far)

    subp.call(["./build/main", outfile])
    t1 = time.time()
    print("\ntime={}\n".format(t1-t0))
    t0 = time.time()
    subp.call(["rm", outfile])

charm_defaults = {"mv" : 30, "mdm" : 10, "mdm1" : 10, "mdm2" : 15, 'channels' : {_brem,_pion_decay,_eta_decay}, 'det_switch' : 'charm'}

def charm_eval(d_user):
    d = charm_defaults.copy()
    d.update(d_user)
    MV=d["mv"]
    MX=d["mdm"]
    MX1=d["mdm1"]
    MX2=d["mdm2"]
    channels = d["channels"]
    det_switch=d["det_switch"]
    if "alpha_D" in d:
        alD=d["alpha_D"]
    if 'eps' in d:
        eps=d['eps']
    sumlog=d["sumlog"]
    signal_channel = d["signal_chan"]

    t0 = time.time()
    outfile="parameter_run_{0}_{1}.dat".format(str(MV),str(MX))

    BEAM_ENERGY=400
    zmin=max(3*MV/1000.0/BEAM_ENERGY,0.1)
    zmax=1-zmin

    proddist = []
    prodchan = []
    partlistfile = []
    executing=False

    if signal_channel=="Inelastic_Nucleon_Scattering_Baryonic" or signal_channel=="Baryonic_Test":
        print("Need to implement Baryonic!")
    else:
        if ((MX/1000.0<mpi0/2.0 and MV<600.0) or (MV/1000.0<mpi0 and (signal_channel=="Signal_Decay"))) and _pion_decay in channels:
            proddist.append("particle_list")
            prodchan.append("pi0_decay")
            partlistfile.append("data/particle_list_charm.dat")
            executing=True
        if ((MX/1000.0<meta/2.0 and MV<900.0) or ((signal_channel=="Signal_Decay") and MV/1000.0<meta)) and _eta_decay in channels:
            proddist.append("particle_list")
            prodchan.append("eta_decay")
            partlistfile.append("data/particle_list_charm.dat")
            executing=True
        if MV/1000.0>=mrho and _parton in channels:
            proddist.append("parton_V")
            prodchan.append("parton_production")
            partlistfile.append("")
            executing=True
        if (MV/2.0>MX or signal_channel=="Signal_Decay") and _brem in channels:
            proddist.append("proton_brem")
            prodchan.append("V_decay")
            partlistfile.append("")
            executing=True
        if MV/1000.0>=mrho and MV<=1250 and _phi_decay in channels:
            proddist.append("particle_list")
            prodchan.append("phi_decay")
            partlistfile.append("data/particle_list_charm.dat")
            executing=True

    if not executing:
        print("No valid channels, skipping!")
        return


    d.update({"proddist" : proddist, "prod_chan" : prodchan, "partlistfile" : partlistfile,"mv" : MV/1000.0, "mdm" : MX/1000.0, "mdm1" : MX1/1000.0, "mdm2" : MX2/1000.0,"zmin" : zmin, "zmax" : zmax, "outfile" : outfile})
    if det_switch == "charm":
        write_charm_decay(d=d)
    subp.call(["./build/main", outfile])
    t1 = time.time()
    print("\ntime={}\n".format(t1-t0))
    t0 = time.time()
    subp.call(["rm", outfile])

ship_defaults = {"mv" : 30, "mdm" : 10, "mdm1" : 10, "mdm2" : 10, 'channels' : {_parton,_brem,_pion_decay,_eta_decay,_phi_decay}, "signal_chan" : "NCE_electron", 'det_switch' : 'ship', 'sumlog' : "Events/ship.dat"}

def ship_eval(d_user):
    d = ship_defaults.copy()
    d.update(d_user)
    MV=d["mv"]
    MX1=d["mdm1"]
    MX2=d["mdm2"]
    MX=d["mdm"]
    channels = d["channels"]
    det_switch=d["det_switch"]
    if "alpha_D" in d:
        alD=d["alpha_D"]
    if 'eps' in d:
        eps=d['eps']
    sumlog=d["sumlog"]
    signal_channel = d["signal_chan"]
    model=d["model"]

    t0 = time.time()
    outfile="parameter_run_{0}_{1}.dat".format(str(MV),str(MX))

    BEAM_ENERGY=400
    zmin=max(3*MV/1000.0/BEAM_ENERGY,0.1)
    zmax=1-zmin

    proddist = []
    prodchan = []
    partlistfile = []
    executing=False

    if model=="Inelastic_Dark_Matter":
        if MX1+MX2<mpi0*1000 and _pion_decay in channels:
            proddist.append("particle_list")
            prodchan.append("pi0_decay")
            partlistfile.append("data/particle_list_ship.dat")
            executing=True
        if MX1+MX2<meta*1000 and _eta_decay in channels:
            proddist.append("particle_list")
            prodchan.append("eta_decay")
            partlistfile.append("data/particle_list_ship.dat")
            executing=True
        if (MV>MX1+MX2) and _brem in channels:
            proddist.append("proton_brem")
            prodchan.append("V_decay")
            partlistfile.append("")
            executing=True
    else:
        if ((MX/1000.0<mpi0/2.0 and MV<600.0) or (MV/1000.0<mpi0 and signal_channel=="Signal_Decay")) and _pion_decay in channels:
            proddist.append("particle_list")
            prodchan.append("pi0_decay")
            partlistfile.append("data/particle_list_ship.dat")
            executing=True
        if ((MX/1000.0<meta/2.0 and MV<900.0) or (signal_channel=="Signal_Decay" and MV/1000.0<meta)) and _eta_decay in channels:
            proddist.append("particle_list")
            prodchan.append("eta_decay")
            partlistfile.append("data/particle_list_ship.dat")
            executing=True
        if MV/1000.0>=mrho and _parton in channels:
            proddist.append("parton_V")
            prodchan.append("parton_production")
            partlistfile.append("")
            executing=True
        if (MV/2.0>MX or signal_channel=="Signal_Decay") and _brem in channels:
            proddist.append("proton_brem")
            prodchan.append("V_decay")
            partlistfile.append("")
            executing=True
        if MV/1000.0>=mrho and MV<=1250 and _phi_decay in channels:
            proddist.append("particle_list")
            prodchan.append("phi_decay")
            partlistfile.append("data/particle_list_ship.dat")
            executing=True
    if not executing:
        print("No valid channels, skipping!")
        return

    d.update({"proddist" : proddist, "prod_chan" : prodchan, "partlistfile" : partlistfile,"mv" : MV/1000.0, "mdm" : MX/1000.0, "mdm1" : MX1/1000.0, "mdm2" : MX2/1000.0, "zmin" : zmin, "zmax" : zmax, "outfile" : outfile})
    if det_switch == "ship":
        write_ship(d=d)
    elif det_switch == "test_sphere":
        write_ship(d=d,det=test_sphere)
    elif det_switch == "test_cylinder":
        write_ship(d=d,det=test_cylinder)
    elif det_switch == "test_cuboid":
        write_ship(d=d,det=test_cuboid)
    elif det_switch == "charm2":
        write_charm2(d=d)
    elif det_switch == "na62":
        write_ship(d=d,det=NA62_decay_vol)
    subp.call(["./build/main", outfile])
    t1 = time.time()
    print("\ntime={}\n".format(t1-t0))
    t0 = time.time()
    subp.call(["rm", outfile])

charm2_defaults = {"mv" : 30, "mdm" : 10, 'channels' : {_parton,_brem,_pion_decay,_eta_decay,_phi_decay}, "signal_chan" : "NCE_electron", 'det_switch' : 'charm', 'sumlog' : "Events/charm2.dat"}

def charm2_eval(d_user):
    d = charm2_defaults.copy()
    d.update(d_user)
    MV=d["mv"]
    MX=d["mdm"]
    channels = d["channels"]
    det_switch=d["det_switch"]
    if "alpha_D" in d:
        alD=d["alpha_D"]
    if 'eps' in d:
        eps=d['eps']
    sumlog=d["sumlog"]
    signal_channel = d["signal_chan"]

    t0 = time.time()
    outfile="parameter_run_{0}_{1}.dat".format(str(MV),str(MX))

    BEAM_ENERGY=450
    zmin=max(3*MV/1000.0/BEAM_ENERGY,0.1)
    zmax=1-zmin

    proddist = []
    prodchan = []
    partlistfile = []
    executing=False

    if signal_channel=="Inelastic_Nucleon_Scattering_Baryonic" or signal_channel=="Baryonic_Test":
        print("Need to implement Baryonic!")
    else:
        if ((MX/1000.0<mpi0/2.0 and MV<600.0) or (MV/1000.0<mpi0 and signal_channel=="Signal_Decay")) and _pion_decay in channels:
            proddist.append("particle_list")
            prodchan.append("pi0_decay")
            partlistfile.append("data/particle_list_charm2.dat")
            executing=True
        if ((MX/1000.0<meta/2.0 and MV<900.0) or (signal_channel=="Signal_Decay" and MV/1000.0<meta)) and _eta_decay in channels:
            proddist.append("particle_list")
            prodchan.append("eta_decay")
            partlistfile.append("data/particle_list_charm2.dat")
            executing=True
        if MV/1000.0>=mrho and _parton in channels:
            proddist.append("parton_V")
            prodchan.append("parton_production")
            partlistfile.append("")
            executing=True
        if (MV/2.0>MX or signal_channel=="Signal_Decay") and _brem in channels:
            proddist.append("proton_brem")
            prodchan.append("V_decay")
            partlistfile.append("")
            executing=True
        if MV/1000.0>=mrho and MV<=1250 and _phi_decay in channels:
            proddist.append("particle_list")
            prodchan.append("phi_decay")
            partlistfile.append("data/particle_list_charm2.dat")
            executing=True

    if not executing:
        print("No valid channels, skipping!")
        return

    d.update({"proddist" : proddist, "prod_chan" : prodchan, "partlistfile" : partlistfile,"mv" : MV/1000.0, "mdm" : MX/1000.0, "zmin" : zmin, "zmax" : zmax, "outfile" : outfile})
    if det_switch == "charm2":
        write_charm2(d=d)
    subp.call(["./build/main", outfile])
    t1 = time.time()
    print("\ntime={}\n".format(t1-t0))
    t0 = time.time()
    subp.call(["rm", outfile])

nucal_defaults = {"mv" : 30, "mdm" : 10, 'channels' : { _brem, _pion_decay, _eta_decay}, "signal_chan" : "Signal_Decay", "det_switch" : 'nucal', "alpha_D" : 0.5, "mdm1" : 10, "mdm2" : 10}

def nucal_eval(d_user):
    d = nucal_defaults.copy()
    d.update(d_user)
    MV=d["mv"]
    MX=d["mdm"]
    MX1=d["mdm1"]
    MX2=d["mdm2"]
    channels = d["channels"]
    det_switch=d["det_switch"]
    if "alpha_D" in d:
        alD=d["alpha_D"]
    if 'eps' in d:
        eps=d['eps']
    sumlog=d["sumlog"]
    signal_channel = d["signal_chan"]
    model=d["model"]

    t0 = time.time()
    outfile="parameter_run_{0}_{1}_{2}.dat".format(str(MV),str(MX),str(eps))
    BEAM_ENERGY=70
    zmin=max(3*MV/1000.0/BEAM_ENERGY,0.1)
    zmax=1-zmin

    proddist = []
    prodchan = []
    partlistfile = [    ]
    executing=False
    if signal_channel=="Inelastic_Nucleon_Scattering_Baryonic" or signal_channel=="Baryonic_Test":
        if MX/1000.0<mpi0/2.0 and MV<1000 and _pion_decay in channels:
            proddist.append("particle_list")
            prodchan.append("pi0_decay_baryonic")
            partlistfile.append("data/particle_list_nucal.dat")
            executing=True
        if MX/1000.0<meta/2.0 and MV<1000 and _eta_decay in channels:
            proddist.append("particle_list")
            prodchan.append("eta_decay_baryonic")
            partlistfile.append("data/particle_list_nucal.dat")
            executing=True
        #if ((MV<1200) and (MV>=350)) and rho_decay_switch:
        #    proddist.append("particle_list")
        #    prodchan.append("omega_decay_baryonic")
        #    partlistfile.append("data/particle_list.dat")
        if MV/1000.0>0.9 and MV<1250 and _phi_decay in channels:
            proddist.append("particle_list")
            prodchan.append("phi_decay_baryonic")
            partlistfile.append("data/particle_list_nucal.dat")
            executing=True
        if MV>2.0*MX and _brem in channels:
            proddist.append("proton_brem_baryonic")
            prodchan.append("V_decay_baryonic")
            partlistfile.append("")
            executing=True
        if MV/1000.0>=mrho and _parton in channels:
            proddist.append("parton_V_baryonic")
            prodchan.append("parton_production_baryonic")
            partlistfile.append("")
            executing=True
    elif model=="Inelastic_Dark_Matter":
        if MX1+MX2<mpi0*1000 and _pion_decay in channels:
            proddist.append("particle_list")
            prodchan.append("pi0_decay")
            partlistfile.append("data/particle_list_nucal.dat")
            executing=True
        if MX1+MX2<meta*1000 and _eta_decay in channels:
            proddist.append("particle_list")
            prodchan.append("eta_decay")
            partlistfile.append("data/particle_list_nucal.dat")
            executing=True
        if (MV>MX1+MX2) and _brem in channels:
            proddist.append("proton_brem")
            prodchan.append("V_decay")
            partlistfile.append("")
            executing=True
    else:
        if ((MX/1000.0<mpi0/2.0 and MV<600.0 and signal_channel!="Signal_Decay") or (MV/1000.0<mpi0 and signal_channel=="Signal_Decay")) and _pion_decay in channels:
            proddist.append("particle_list")
            prodchan.append("pi0_decay")
            partlistfile.append("data/particle_list_nucal.dat")
            executing=True
        if ((signal_channel!="Signal_Decay" and MX/1000.0<meta/2.0 and MV<900.0) or (signal_channel=="Signal_Decay" and MV/1000.0<meta)) and _eta_decay in channels:
            proddist.append("particle_list")
            prodchan.append("eta_decay")
            partlistfile.append("data/particle_list_nucal.dat")
            executing=True
        if MV/1000.0>=mrho and _parton in channels:
            proddist.append("parton_V")
            prodchan.append("parton_production")
            partlistfile.append("")
            executing=True
        if (MV/2.0>MX or signal_channel=="Signal_Decay") and _brem in channels:
            proddist.append("proton_brem")
            prodchan.append("V_decay")
            partlistfile.append("")
            executing=True
        if MV/1000.0>=mrho and MV<=1250 and _phi_decay in channels:
            proddist.append("particle_list")
            prodchan.append("phi_decay")
            partlistfile.append("data/particle_list_nucal.dat")
            executing=True
   # if MV/1000.0>=mrho and MV<=1250:
   #     proddist.append("particle_list")
   #     prodchan.append("phi_decay")
   #     partlistfile.append("data/particle_list.dat")
    if not executing:
        return
    d.update(d_user)
    d.update({"proddist" : proddist, "prod_chan" : prodchan, "partlistfile" : partlistfile,"mv" : MV/1000.0, "mdm" : MX/1000.0, "mdm1" : MX1/1000.0, "mdm2" : MX2/1000.0, "zmin" : zmin, "zmax" : zmax, "outfile" : outfile})
    if det_switch == "nucal":
        write_nucal(d=d)
    subp.call(["./build/main", outfile])
    t1 = time.time()
    print("\ntime={}\n".format(t1-t0))
    t0 = time.time()
    subp.call(["rm", outfile])

miniboone_defaults = {"mv" : 30, "mdm" : 10, 'channels' : {_parton,_brem,_pion_decay,_eta_decay}, "signal_chan" : "NCE_nucleon", 'det_switch' : 'miniboone', 'sumlog' : "Events/miniboone_y.dat", "model" : "Dark_Photon", 'eps' : 1e-3, 'alpha_D' : 0.5}

def miniboone_eval(d_user):
    d = miniboone_defaults.copy()
    d.update(d_user)
    MV=d["mv"]
    MX=d["mdm"]
    channels = d["channels"]
    det_switch=d["det_switch"]
    if "alpha_D" in d:
        alD=d["alpha_D"]
    if 'eps' in d:
        eps=d['eps']
    sumlog=d["sumlog"]
    signal_channel = d["signal_chan"]
    model = d["model"]

    t0 = time.time()
    outfile="parameter_run_{0}_{1}.dat".format(str(MV),str(MX))

    zmin=max(3*MV/1000.0/8.9,0.3)
    zmax=1-zmin

    proddist = []
    prodchan = []
    partlistfile = [    ]
    executing=False
    if model == "Dark_Photon_Baryonic":
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
    elif model == "Dark_Photon" or model == "Dark_Photon2" or model == "Axion_Dark_Photon":
        if ((MX/1000.0<mpi0/2.0) or (signal_channel=="Signal_Decay" and MV/1000.0<mpi0)) and _pion_decay in channels:
            proddist.append("particle_list")
            prodchan.append("pi0_decay")
            if use_miniboone_sample:
                partlistfile.append("data/offTarget_pi0Avg.txt")
            else:
                partlistfile.append("data/particle_list.dat")
            executing=True
        if ((MX/1000.0<meta/2.0) or (signal_channel=="Signal_Decay" and MV/1000.0<meta)) and _eta_decay in channels:
            proddist.append("particle_list")
            prodchan.append("eta_decay")
            if use_miniboone_sample:
                partlistfile.append("data/offTarget_etaAvg.txt")
            else:
                partlistfile.append("data/particle_list_k0.dat")
            executing=True
        if MV/1000.0>=mrho and MV>2*MX and _parton in channels:
	    proddist.append("parton_V")
	    prodchan.append("parton_production")
	    partlistfile.append("")
            executing=True
        if (MV/2.0>MX or signal_channel=="Signal_Decay") and zmin<zmax and _brem in channels:
            proddist.append("proton_brem")
            prodchan.append("V_decay")
            partlistfile.append("")
            executing=True
        if MV/1000.0>=mrho and MV<=1250 and _phi_decay in channels:
            proddist.append("particle_list")
            prodchan.append("phi_decay")
            partlistfile.append("data/particle_list.dat")
    if not executing:
        return

    if signal_channel=="NCE_nucleon":
        #user2 = {"samplesize" : 5000, "min_scatter_energy" : 0.035, "max_scatter_energy" : 2, "efficiency" : 0.35, "outlog" : "Events/miniboone_nucleon_events_{}_{}.dat".format(str(MV),str(eps))}
        user2 = {"samplesize" : 1000, "min_scatter_energy" : 0.035, "max_scatter_energy" : 1, "efficiency" : 0.35, "outlog" : "Events/miniboone_nucleon_events_{}_{}.dat".format(str(MV),str(MX))}
    elif signal_channel in _pion_inelastic_channels:
        user2 = {"outlog" : "Events/miniboone_pion_events_{}_{}.dat".format(str(MV),str(MX)), "efficiency" : 0.35, "min_scatter_energy" : 0.0, "max_scatter_energy" : 10.0}
    elif signal_channel == "NCE_electron":
        user2 = {"min_scatter_energy" : 0.0, "max_scatter_energy" : 1.0, "max_scatter_angle" : 0.14, "outlog" : "Events/miniboone_electron_events_{}_{}.dat".format(str(MV),str(MX))}
    elif signal_channel == "Signal_Decay":
        outfile="parameter_run_{0}_{1}.dat".format(str(MV),str(eps))
        user2 = {"min_scatter_energy" : 0.0, "max_scatter_energy" : 10.0, "outlog" : "Decay_Events/miniboone_decay_events_{}_{}.dat".format(str(MV),str(eps)), "efficiency" : 1, "outfile" : outfile}
        if det_switch == "sbnd":
            user2["outlog"] = "Events/sbnd_decay_events_{}_{}".format(str(MV),str(eps))

    d.update(user2)
    d.update(d_user)
    d.update({"proddist" : proddist, "prod_chan" : prodchan, "partlistfile" : partlistfile,"mv" : MV/1000.0, "mdm" : MX/1000.0, "zmin" : zmin, "zmax" : zmax, "outfile" : outfile})
    if det_switch == "sbnd":
        write_miniboone(d=d,det=SBND_detector)
    elif det_switch == "miniboone_full":
        write_miniboone(d=d,det=miniboone_detector_full)
    elif det_switch == "icarus":
        write_miniboone(d=d,det=ICARUS_detector)
    else:
        write_miniboone(d=d)

    subp.call(["./build/main", outfile])
    t1 = time.time()
    print("\ntime={}\n".format(t1-t0))
    t0 = time.time()
    subp.call(["rm", outfile])
#END OF miniboone_eval

t2k_defaults = {"mv" : 30, "mdm" : 10, "signal_chan" : "Pion_Inelastic","channels" : {_parton,_brem,_pion_decay,_eta_decay,_phi_decay}, 'det_switch' : "nd280", 'sumlog' : "Events/nova_decay.dat"}

def t2k_eval(d_user):
    d = t2k_defaults.copy()
    d.update(d_user)
    MV=d["mv"]
    MX=d["mdm"]
    channels = d["channels"]
    det_switch=d["det_switch"]
    if "alpha_D" in d:
        alD=d["alpha_D"]
    if 'eps' in d:
        eps=d['eps']
    sumlog=d["sumlog"]
    signal_channel = d["signal_chan"]
    model = d["model"]

    t0 = time.time()
    outfile="parameter_run_{0}_{1}.dat".format(str(MV),str(MX))

    zmin=max(3*MV/1000.0/30.0,0.2)
    zmax=1-zmin

    print("Zmin=",zmin,"Zmax=",zmax)

    proddist = []
    prodchan = []
    partlistfile = []
    executing = False
    if model == "Dark_Photon":
        if MX/1000.0<mpi0/2.0 and MV<600.0 and _pion_decay in channels:
            proddist.append("particle_list")
            prodchan.append("pi0_decay")
            partlistfile.append("data/particle_list_t2k.dat")
            executing=True
        if MX/1000.0<meta/2.0 and MV<900.0 and _eta_decay in channels:
            proddist.append("particle_list")
            prodchan.append("eta_decay")
            partlistfile.append("data/particle_list_t2k.dat")
            executing=True
        if MV/1000.0>=mrho and _parton in channels:
            proddist.append("parton_V")
    	    prodchan.append("parton_production")
    	    partlistfile.append("")
            executing=True
        if MV/2.0>MX and _brem in channels:
    	    proddist.append("proton_brem")
    	    prodchan.append("V_decay")
     	    partlistfile.append("")
            executing=True
        if MV/1000.0>=mrho and MV<=1250 and _phi_decay in channels:
            proddist.append("particle_list")
            prodchan.append("phi_decay")
            partlistfile.append("data/particle_list_t2k.dat")
            executing=True
    if not executing:
        return

    if signal_channel=="NCE_nucleon":
        if det_switch == "superk":
            user2 = {"samplesize" : 500, "min_scatter_energy" : 0.0014, "max_scatter_energy" : 1000.0, "efficiency" : 0.66/1000.0, "sumlog" : "Events/superk_y,dat","timing" : 5e-8}
        elif det_switch == "P0D":
            pass
    elif signal_channel in _pion_inelastic_channels:
        if det_switch == "superk":
            user2 = {"efficiency" : 0.66/1000.0, "samplesize" : 1000, "min_scatter_energy" : 0.0, "max_scatter_energy" : 1000.0, "timing" : 5e-8}
        elif det_switch == "P0D":
            #Warning, efficiency is arbitrary right now!
            user2 = {"samplesize" : 1000, "min_scatter_energy" : 0.0, "max_scatter_energy" : 1000.0, "efficiency" : 0.5, "sumlog" : "Events/t2k_pod_y.dat"}
    elif signal_channel == "NCE_electron":
        if det_switch == "superk":
            pass
        elif det_switch == "P0D":
            pass
    d.update(user2)
    d.update(d_user)
    d.update({"proddist" : proddist, "prod_chan" : prodchan, "partlistfile" : partlistfile,"mv" : MV/1000.0, "mdm" : MX/1000.0, "zmin" : zmin, "zmax" : zmax, "outfile" : outfile})

    if det_switch == "superk":
        write_t2k(d, det=t2k_superK1000)
    elif det_switch == "P0D":
        write_t2k(d)

    subp.call(["./build/main", outfile])
    t1 = time.time()
    print("\ntime={}\n".format(t1-t0))
    t0 = time.time()
    subp.call(["rm", outfile])
#END OF t2k_eval

numi_defaults = {"mv" : 30, "mdm" : 10, 'channels' : {_parton,_brem,_pion_decay,_eta_decay,_phi_decay}, "signal_chan" : "Signal_Decay", 'det_switch' : 'nova', 'sumlog' : "Events/nova_decay.dat", 'mdm1' : 10, 'mdm2' : 10}

def numi_eval(d_user):
    d = numi_defaults.copy()

    d.update(d_user)
    MV=d["mv"]
    MX=d["mdm"]
    MX1=d["mdm1"]
    MX2=d["mdm2"]
    channels = d["channels"]
    det_switch=d["det_switch"]
    if "alpha_D" in d:
        alD=d["alpha_D"]
    if 'eps' in d:
        eps=d['eps']
    sumlog=d["sumlog"]
    signal_channel = d["signal_chan"]
    model=d["model"]

    t0 = time.time()
    outfile="parameter_run_{0}_{1}.dat".format(str(MV),str(MX))
    #if MX/1000.0<mpi0/2.0:
    if "DIST" in d and "ANGLE" in d:
        _DIST = d["DIST"]
        _ANGLE = d["ANGLE"]

        global _DET_XPOS, _DET_YPOS, _DET_ZPOS
        if det_switch == "nova" or det_switch == "nova_abs":
            _DET_XPOS = 0
            _DET_YPOS = -_DIST*math.sin(math.radians(_ANGLE))
            _DET_ZPOS = _DIST*math.cos(math.radians(_ANGLE))
        else:
            _DET_XPOS = _DIST*math.sin(math.radians(_ANGLE))
            _DET_YPOS = 0
            _DET_ZPOS = _DIST*math.cos(math.radians(_ANGLE))

    BEAM_ENERGY=120
    zmin=max(3*MV/1000.0/BEAM_ENERGY,0.1)
    zmax=1-zmin

    proddist = []
    prodchan = []
    partlistfile = [    ]
    executing=False
    if signal_channel=="Inelastic_Nucleon_Scattering_Baryonic" or signal_channel=="Baryonic_Test":
        if MX/1000.0<mpi0/2.0 and MV<1000 and _pion_decay in channels:
            proddist.append("particle_list")
            prodchan.append("pi0_decay_baryonic")
            partlistfile.append("data/particle_list_numi_pion.dat")
            executing=True
        if MX/1000.0<meta/2.0 and MV<1000 and _eta_decay in channels:
            proddist.append("particle_list")
            prodchan.append("eta_decay_baryonic")
            partlistfile.append("data/particle_list_numi_eta.dat")
            executing=True
        #if ((MV<1200) and (MV>=350)) and rho_decay_switch:
        #    proddist.append("particle_list")
        #    prodchan.append("omega_decay_baryonic")
        #    partlistfile.append("data/particle_list.dat")
        if MV/1000.0>0.9 and MV<1250 and _phi_decay in channels:
            proddist.append("particle_list")
            prodchan.append("phi_decay_baryonic")
            partlistfile.append("data/particle_list_numi.dat")
            executing=True
        if MV>2.0*MX and _brem in channels:
            proddist.append("proton_brem_baryonic")
            prodchan.append("V_decay_baryonic")
            partlistfile.append("")
            executing=True
        if MV/1000.0>=mrho and _parton in channels:
            proddist.append("parton_V_baryonic")
            prodchan.append("parton_production_baryonic")
            partlistfile.append("")
            executing=True
    elif model=="Inelastic_Dark_Matter":
        if MX1+MX2<mpi0*1000 and _pion_decay in channels:
            proddist.append("particle_list")
            prodchan.append("pi0_decay")
            partlistfile.append("data/particle_list_numi_pion.dat")
            executing=True
        if MX1+MX2<meta*1000 and _eta_decay in channels:
            proddist.append("particle_list")
            prodchan.append("eta_decay")
            partlistfile.append("data/particle_list_numi_eta.dat")
            executing=True
        if (MV>MX1+MX2) and _brem in channels:
            proddist.append("proton_brem")
            prodchan.append("V_decay")
            partlistfile.append("")
            executing=True
    else:
        if ((MX/1000.0<mpi0/2.0 and MV<600.0) or (MV/1000.0<mpi0 and signal_channel=="Signal_Decay")) and _pion_decay in channels:
            proddist.append("particle_list")
            prodchan.append("pi0_decay")
            partlistfile.append("data/particle_list_numi_pion.dat")
            executing=True
        if ((MX/1000.0<meta/2.0 and MV<900.0) or (signal_channel=="Signal_Decay" and MV/1000.0<meta)) and _eta_decay in channels:
            proddist.append("particle_list")
            prodchan.append("eta_decay")
            partlistfile.append("data/particle_list_numi_eta.dat")
            executing=True
        if MV/1000.0>=mrho and _parton in channels:
            proddist.append("parton_V")
    	    prodchan.append("parton_production")
    	    partlistfile.append("")
            executing=True
        if (MV/2.0>MX or signal_channel=="Signal_Decay") and _brem in channels:
    	    proddist.append("proton_brem")
    	    prodchan.append("V_decay")
     	    partlistfile.append("")
            executing=True
        if MV/1000.0>=mrho and MV<=1250 and _phi_decay in channels:
            proddist.append("particle_list")
            prodchan.append("phi_decay")
            partlistfile.append("data/particle_list_numi_eta.dat")
            executing=True
   # if MV/1000.0>=mrho and MV<=1250:
   #     proddist.append("particle_list")
   #     prodchan.append("phi_decay")
   #     partlistfile.append("data/particle_list.dat")
    if not executing:
        return
    if signal_channel == "Signal_Decay":
        outfile="parameter_run_{0}_{1}_{2}.dat".format(str(MV),str(eps),str(det_switch))
        user2 = {"min_scatter_energy" : 0.0, "max_scatter_energy" : 120.0,  "efficiency" : 1, "outfile" : outfile}
        if det_switch == "nova":
            user2["sumlog"] = "Decay_Events/nova_decay.dat"
            user2["outlog"] = "Decay_Events/nova_decay_events_{}_{}.dat".format(str(MV),str(eps))
        elif det_switch == "nova_absorber":
            user2["sumlog"] = "Decay_Events/nova_abs_decay.dat"
            user2["outlog"] = "Decay_Events/nova_abs_decay_events_{}_{}.dat".format(str(MV),str(eps))
        elif det_switch == "minos":
            user2["outlog"] = "Decay_Events/minos_decay_events_{}_{}.dat".format(str(MV),str(eps))
            user2["sumlog"] = "Decay_Events/minos_decay.dat"
        elif det_switch == "minos_absorber":
            user2["outlog"] = "Decay_Events/minos_abs_decay_events_{}_{}.dat".format(str(MV),str(eps))
            user2["sumlog"] = "Decay_Events/minos_abs_decay.dat"
        elif det_switch == "miniboone_numi":
            user2["outlog"] = "Decay_Events/mini_numi_decay_events_{}_{}.dat".format(str(MV),str(eps))
            user2["sumlog"] = "Decay_Events/mini_numi_decay.dat"
        elif det_switch == "test_sphere" or det_switch == "test_cuboid" or det_switch == "test_cylinder":
            user2["outlog"] = det_switch
            user2["sumlog"] = "test.dat"
        d.update(user2)
    d.update(d_user)
    d.update({"proddist" : proddist, "prod_chan" : prodchan, "partlistfile" : partlistfile,"mv" : MV/1000.0, "mdm" : MX/1000.0, "mdm1" : MX1/1000.0, "mdm2" : MX2/1000.0, "zmin" : zmin, "zmax" : zmax, "outfile" : outfile})
    if det_switch == "nova":
        write_numi(d=d)
    elif det_switch == "nova_absorber":
        write_numi(d=d,det=NOvA_absorber_detector)
    elif det_switch == "minos":
        write_numi(d=d,det=MINOS_detector)
    elif det_switch == "minos_absorber":
        write_numi(d=d,det=MINOS_absorber_detector)
    elif det_switch == "miniboone_numi":
        write_numi(d=d,det=miniboone_detector_numi)
    elif det_switch == "icarus":
        write_numi(d=d,det=ICARUS_detector_NuMI)
    elif det_switch=="dune_mpd":
        write_numi(d=d,det=dune_mpd_detector)
    subp.call(["./build/main", outfile])
    t1 = time.time()
    print("\ntime={}\n".format(t1-t0))
    t0 = time.time()
    subp.call(["rm", outfile])
#####END OF numi_eval############

def coherent_eval(d_user):
    d = numi_defaults.copy()
    d.update(d_user)
    MV=d["mv"]
    MX=d["mdm"]
    channels = d["channels"]
    det_switch=d["det_switch"]
    if "alpha_D" in d:
        alD=d["alpha_D"]
    if 'eps' in d:
        eps=d['eps']
    sumlog=d["sumlog"]
    signal_channel = d["signal_chan"]

    t0 = time.time()
    outfile="parameter_run_{0}_{1}.dat".format(str(MV),str(MX))

    proddist = []
    prodchan = []
    partlistfile = [    ]
    executing=False
    if signal_channel=="NCE_nucleon_baryonic":
        if MX/1000.0<mpi0/2.0 and MV<600.0 and _pion_decay in channels:
            proddist.append("particle_list")
            prodchan.append("pi0_decay_baryonic")
            partlistfile.append("data/particle_list_coherent.dat")
            executing = True
    else:
        if MX/1000.0<mpi0/2.0 and _pion_decay in channels:
            proddist.append("particle_list")
            prodchan.append("pi0_decay")
            partlistfile.append("data/particle_list_coherent.dat")
            executing = True
        if MX/1000.0<0.129/2.0 and MV<600 and MV>2*MX and _piminus_cap in channels:
            proddist.append("")
            prodchan.append("piminus_capture")
            partlistfile.append("")
            executing = True
    if not executing:
        return
    _captain_dets = ["captain","captain_off"]
    _coherent_dets = ["csi",'csi1T']
    if signal_channel=="NCE_nucleon":
        if det_switch in _captain_dets:
            user2 = {"samplesize" : 1000, "min_scatter_energy" : 0.03, "max_scatter_energy" : 0.2, "efficency" : 0.5, "coherent" : "false", "burn_max" : 1000, "sumlog" : "Events/captain.dat"}
        else:
                user2 = {"samplesize" : 1000, "min_scatter_energy" : 5e-6, "max_scatter_energy" : 0.2, "efficiency" : 0.5, "sumlog" : "Events/coherent.dat", "coherent" : "true", "eps" : 1e-3, "burn_max" : 1000}
    elif signal_channel == "NCE_nucleon_baryonic":
        user2 = {"samplesize" : 1000, "min_scatter_energy" : 5e-6, "max_scatter_energy" : 0.2, "efficiency" : 0.5, "sumlog" : "Events/coherent_baryonic.dat", "eps" : 0.0, "coherent" : "true", "alpha_D" : 1e-4, "burn_max" : 1000}
    if signal_channel=="NCE_electron":
        if det_switch in _captain_dets:
            user2 = {"samplesize" : 1000, "min_scatter_energy" : 0.01, "max_scatter_energy" : 0.2, "efficency" : 0.5, "coherent" : "false", "burn_max" : 1000, "sumlog" : "Events/captain.dat"}
    if signal_channel in _pion_inelastic_channels:
        if det_switch in _captain_dets:
            user2 = {"samplesize" : 1000, "burn_max" : 1000, "min_scatter_energy" : 0.0, "max_scatter_energy" : 1, "efficency" : 0.5, "coherent" : "false", "burn_max" : 1000, "sumlog" : "Events/captain.dat"}
    d.update(user2)
    d.update(d_user)
    d.update({"proddist" : proddist, "prod_chan" : prodchan, "partlistfile" : partlistfile,"mv" : MV/1000.0, "mdm" : MX/1000.0, "outfile" : outfile})
    if det_switch == "csi":
        write_coherent(d=d,det=coherent_detector_CsI)
    elif det_switch == 'csi1T':
        write_coherent(d=d,det=coherent_detector_CsI_1T)
    elif det_switch == "captain":
        write_coherent(d=d,det=captain_detector)
    elif det_switch == "captain_off_axis":
        write_coherent(d=d,det=captain_detector_off)
    elif det_switch == "LAr29":
        write_coherent(d=d,det=coherent_detector_LAr_29kg)
    else:
        write_coherent(d=d)

    subp.call(["./build/main", outfile])
    t1 = time.time()
    print("\ntime={}\n".format(t1-t0))
    t0 = time.time()
    subp.call(["rm", outfile])
#END OF coherent_eval

use_miniboone_sample=True
def execute_miniboone_parallel(genlist = True):
    if genlist and not use_miniboone_sample:
        d={"prod_chan" : ["pi0_decay"],"proddist" : ["pi0_sanfordwang"],"samplesize" : 2e6,"output_mode" : "particle_list","partlistfile" : ["data/particle_list.dat"]}
        write_miniboone(d=d)
        subp.call(["./build/main", "parameter_run.dat"])
        d={"prod_chan" : ["eta_decay"],"proddist" : ["k0_sanfordwang"],"samplesize" : 2e6,"output_mode" : "particle_list","partlistfile" : ["data/particle_list_k0.dat"]}
        write_miniboone(d=d)
        subp.call(["./build/main", "parameter_run.dat"])
    channs={_brem,_pion_decay,_eta_decay}
    #vmassarr=[2,2.5,3,4,5,10,15,17,19,19.5,20.5,21,23,25,30,40,50,60,70,80,90,95,100,105,110,115,120,125,130,132,134,135,137,140,150,160,180,200,220,250,300,350,400,500,530,540,550,700,760,770,773,777,779,790,800,900,1000,1500,2000,3000,5000,7500,10000]
    #vmarr={10,100,300,500,770}
    #epsarr={1e-6,3e-7,1e-7,3e-7,1e-8}
    epsarr=[1e-3]
    #massarr=[[mv,10,eps] for mv in vmassarr for eps in epsarr]
    massarr=[[10,10/3.0,1e-3],[30,10,1e-3],[300,100,1e-3]]
    for marr in massarr:
        #d={"mv" : marr[0],"mdm" : marr[1], "eps" : marr[2], "signal_chan" : "Signal_Decay", "output_mode" : "comprehensive", "det_switch" : "miniboone_full", "sumlog" : "Decay_Events/miniboone_decay.dat", "samplesize" : 1000}
        #d={"mv" : marr[0],"mdm" : marr[1], "eps" : marr[2], "signal_chan" : "NCE_electron", "output_mode" : "summary", "det_switch" : "sbnd", "alpha_D" : 0.1, "channels" : channs, "model" : "Dark_Photon", "sumlog" : "Events/sbnd_elec.dat","POT" : 1e21, "samplesize" : 2000}
        #miniboone_eval(d)
        #d.update({"alpha_D" : 0.5})
        #miniboone_eval(d)
        #d={"mv" : marr[0],"mdm" : marr[1], "eps" : marr[2], "signal_chan" : "Pion_Inelastic", "output_mode" : "summary", "det_switch" : "sbnd", "alpha_D" : 0.1, "channels" : channs, "model" : "Dark_Photon", "sumlog" : "Events/sbnd_pi0.dat", "min_scatter_energy" : 0, "max_scatter_energy" : 100, "min_scatter_angle" : 0, "max_scatter_angle" : 10, "POT" : 1e21, "samplesize" : 2000}
        #miniboone_eval(d)
        #d.update({"alpha_D" : 0.5})
        #miniboone_eval(d)
        #d.update({"mv" : marr[0],"mdm" : marr[1], "eps" : marr[2], "signal_chan" : "NCE_nucleon", 'min_scatter_energy' : 0.035, "output_mode" : "comprehensive", 'sumlog' : "Claudia/sbnd.dat", 'outlog' :  "Claudia/sbnd_nucleon_{}_{}.dat".format(marr[0],marr[1])})
        #miniboone_eval(d)
        #d.update({"mv" : marr[0],"mdm" : marr[1], "eps" : marr[2], "signal_chan" : "NCE_nucleon", 'min_scatter_energy' : 0, "output_mode" : "comprehensive", 'sumlog' : "Claudia/sbnd.dat", 'outlog' :  "Claudia/sbnd_nucleon_no_cut_{}_{}.dat".format(marr[0],marr[1])})
        #miniboone_eval(d)
        #d.update({"mv" : marr[0],"mdm" : marr[1], "eps" : marr[2], "signal_chan" : "NCE_nucleon", "output_mode" : "dm_detector_distribution", "channels" : channs, 'sumlog' : "Claudia/sbnd.dat", 'outlog' :  "Claudia/sbnd_dm_{}_{}.dat".format(marr[0],marr[1])})
        #miniboone_eval(d)
        d={"mv" : marr[0],"mdm" : marr[1], "eps" : marr[2], "signal_chan" : "NCE_electron", "output_mode" : "summary", "det_switch" : "icarus", "alpha_D" : 0.1, "channels" : channs, "model" : "Dark_Photon", "sumlog" : "Events/icarus.dat","POT" : 1e21, "samplesize" : 20000, "pi0_per_POT" : 2.42}
        miniboone_eval(d)

pythia_120=False
def execute_numi(genlist=True):
    #pi0_per_POT_numi=4.176
    #meson_per_pi0_numi=meson_per_pi0_ship = {'pi0_decay' : '1.0', 'eta_decay' : str(0.474/4.176), 'rho_decay' : str(0.11), 'omega_decay' : '0.11', 'phi_decay' : str(0.02)}
    if genlist and not pythia_120:
        d={"prod_chan" : ["pi0_decay"],"proddist" : ["bmpt"],"samplesize" : 2e6,"output_mode" : "particle_list","partlistfile" : ["data/particle_list_numi_pion.dat"]}
        write_numi(d=d)
        subp.call(["./build/main", "parameter_run.dat"])
        d={"prod_chan" : ["pi0_decay"],"proddist" : ["bmpt"],"samplesize" : 2e6,"output_mode" : "particle_list","partlistfile" : ["data/particle_list_numi_eta.dat"]}
        write_numi(d=d)
        subp.call(["./build/main", "parameter_run.dat"])
    elif genlist and pythia_120:
        arr1=np.loadtxt("data/pion_pythia_120GeV.dat")
        np.random.shuffle(arr1)
        np.savetxt("data/particle_list_numi_pion.dat", arr1)
        arr1=np.loadtxt("data/eta_pythia_120GeV.dat")
        np.random.shuffle(arr1)
        np.savetxt("data/particle_list_numi_eta.dat", arr1)
        meson_per_pi0_numi=meson_per_pi0_ship = {'pi0_decay' : '1.0', 'eta_decay' : str(0.474/4.176), 'rho_decay' : str(0.11), 'omega_decay' : '0.11', 'phi_decay' : str(0.02)}
    #vmarr=[10,25,50,75,100,125,150,175,200,250,300,350,400]
    #vmarr=[10,20,30,40,60,80,100,150,200,250,300,400,500,550,600,700,800,900,1000]
    #vmarr=[10,30,60,100,150,200,250,300,400,500,600,700,800,900,1000]
    #epsarr=[10**n for n in range(-8,-3)]+[3*10**n for n in range(-9,-4)]
    #vmarr=[1.15,1,1.2,1.3,1.5,2,3,4,5,10,15,20,30,40,60,80,100,130,140,150,200,250,300,400,500,540,550,560,600,700,725,750,760,765,767,770,772,775,780,790,800,900,1000,1100,1200,1300,1400,1500,1750,2000,2250,2500,2750,3000,3500,4000,4500,5000]
    #vmarr=[1.15,1.2,1.3,1.5,2,3,4,5,10,15,20,30,40,60,80,100,130,140,150,200,250,300,400,500,540,550,560,600,700,725,750,760,765,767,770,772,775,780,790,800,850,900,950,1000]
    vmarr=[1310,1320,1330,1340,1350,1360,1370,1380,1390]
    massarr=[[mv,mv] for mv in vmarr]
    #vmarr=[50]
    #epsarr=[10**-7]
    #massarr=[[mv,mv/3.0,eps] for mv in vmarr for eps in epsarr]
    #massarr=[[10,10/3.0,1e-3],[30,10,1e-3],[300,100,1e-3]]
    #massarr=[[mv,mv,eps] for mv in vmarr for eps in epsarr]
    #massarr=[[mv,10,eps] for mv in vmarr for eps in epsarr]
    #massarr=[[1000,300,1e-3]]
    #d=({"signal_chan" : "Signal_Decay", "output_mode" : "summary", "samplesize" : 1000, "model" : "Dark_Photon","min_scatter_energy" : 5, "min_scatter_angle" : 0});
    #d=({"signal_chan" : "NCE_electron", "output_mode" : "summary", "samplesize" : 1000, "min_scatter_energy" : 5, "max_scatter_energy" : 35, "efficiency" : 0.5, "alpha_D" : 0.5, "POT" : 6e20});
    d_list=[]
    #d=({"channels" : [_pion_decay,_eta_decay,_brem], "signal_chan" : "DM2_Signal_Decay", "output_mode" : "comprehensive", "samplesize" : 50000,  "sumlog" : "IDM_Events/dune_decay.dat", "model" : "Inelastic_Dark_Matter", "min_scatter_energy" : 0, "max_scatter_energy" : 120, "min_scatter_angle" : 0, "max_scatter_angle" : 6, "det_switch" : "dune_mpd", "eps" : 1e-4, "weighted" : 'true'});
    d=({"channels" : [_pion_decay,_eta_decay,_brem], "signal_chan" : "Signal_Decay", "output_mode" : "comprehensive", "samplesize" : 50000,  "sumlog" : "Visible_Dark_Photon/dune_decay.dat", "model" : "Dark_Photon", "min_scatter_energy" : 0, "max_scatter_energy" : 120, "min_scatter_angle" : 0, "max_scatter_angle" : 6, "det_switch" : "dune_mpd", "eps" : 1e-8, "weighted" : 'true','POT' : 1.47e22});
    for marr in massarr:
        print(marr[0],marr[1])
        d.update({"mv" : marr[0], "mdm" : marr[1], "outlog" : "Visible_Dark_Photon/dune/dune_{}.dat".format(marr[0]),"alpha_D" : 0})
        #d.update({"mv" : marr[0], "mdm1" : marr[1], "mdm2" : 1.05*marr[1], "outlog" : "IDM_Events/dune_aD0.5_Delta0.05/dune_{}_{}.dat".format(marr[0],marr[1]),"alpha_D" : 0.5})
        #numi_eval(d)
        #d.update({"mv" : marr[0],"mdm1" : marr[1], "mdm2" : 1.1*marr[1], "outlog" : "IDM_Events/dune_aD0.1_Delta0.1/dune_{}_{}.dat".format(marr[0],marr[1]),"alpha_D" : 0.1})
        #numi_eval(d)
        #d.update({"mv" : marr[0],"mdm1" : marr[1], "mdm2" : 1.4*marr[1], "outlog" : "IDM_Events/dune_aD0.1_Delta0.4/dune_{}_{}.dat".format(marr[0],marr[1]),"alpha_D" : 0.1})
        numi_eval(d)
    #for marr in massarr:
    #    d={"mv" : marr[0],"mdm" : marr[1], "eps" : marr[2], "signal_chan" : "NCE_electron", "output_mode" : "summary", "det_switch" : "icarus", "alpha_D" : 0.1, "channels" : channs, "model" : "Dark_Photon", "sumlog" : "Events/icarus_numi.dat","POT" : 1e21, "samplesize" : 20000, "pi0_per_POT" : pi0_per_POT_numi, "min_scatter_energy" : 0.075, "max_scatter_energy" : 1, "max_scatter_angle" :  0.14, "min_scatter_angle" : 0}
        #d.update({"mv" : marr[0],"mdm" : marr[1], "eps" : marr[2]})
        #d.update({"det_switch" : "nova","channels" : [_pion_decay,_eta_decay,_brem,_parton], "sumlog" : "Events/nova_electron_high.dat", "outlog" : "Events/nova_electron_{}_{}.dat".format(marr[0],marr[1])})
        #d.update({"det_switch" : "nova","channels" : [_pion_decay,_eta_decay,_brem], "sumlog" : "Events/nova_dec.dat"})
        #numi_eval(d)
        #d_list.append(copy.deepcopy(d))
        #d.update({"det_switch" : "nova_absorber","channels" : [_pion_decay,_eta_decay,_brem], "sumlog" : "Events/nova_dec_abs.dat"})
        #d.update({"det_switch" : "nova_absorber","channels" : [_pion_decay,_eta_decay,_brem,_parton], "sumlog" : "Events/nova_electron_abs_high.dat"})
        #numi_eval(d)
        #d_list.append(copy.deepcopy(d))
    #pool = Pool(processes=3)
    #pool.map(numi_eval,d_list)
    '''
    for marr in massarr:
        d.update({"mv" : marr[0],"mdm" : marr[1], "eps" : marr[2]})
        numi_eval(d)
        #d_list.append(copy.deepcopy(d))
        numi_eval(d)
    '''
    '''
    d_list=[]
    for marr in massarr:
        d.update({"mv" : marr[0],"mdm" : marr[1], "eps" : marr[2]})
        d.update({"det_switch" : "minos_absorber","channels" : [_pion_decay,_eta_decay,_brem]})
        #numi_eval(d)
        d_list.append(copy.deepcopy(d))
        d.update({"det_switch" : "minos","channels" : [_pion_decay,_eta_decay,_brem]})
        d_list.append(copy.deepcopy(d))
        #numi_eval(d)
        #d.update({"det_switch" : "nova","channels" : [_pion_decay,_eta_decay,brem]})
        #d_list.append(d)
        #numi_eval(d)
        d.update({"det_switch" : "nova_absorber","channels" : [_pion_decay,_eta_decay,_brem]})
        d_list.append(copy.deepcopy(d))
        #numi_eval(d)
        #d={"model" : "Dark_Photon", "mv" : marr[0],"mdm" : marr[1], "eps" : marr[2], "signal_chan" : "Signal_Decay", "output_mode" : "comprehensive", "det_switch" : "miniboone_numi", "samplesize" : 1000}
        #numi_eval(d)
    pool = Pool(processes=4)
    pool.map(numi_eval,d_list)
   '''


def execute_t2k(genlist=True):
    if genlist:
        d={"prod_chan" : ["pi0_decay"],"proddist" : ["bmpt"],"samplesize" : 2e6,"output_mode" : "particle_list","partlistfile" : ["data/particle_list_numi.dat"]}
        write_t2k(d=d)
        subp.call(["./build/main", "parameter_run.dat"])
    vmarr=[10,50,100,150,200,300,400,500,600,700,760,770,780,800,900,1000,1250,1500,1750,2000]
    vmarr=vmarr+[1100,1200,1300,1400,2250,2500,2750,3000,3250,3500,3750,4000]
    epsarr={1e-3}
    chans={_pion_decay,_eta_decay,_brem,_parton}
    massarr=[[mv,mv/3.0,eps,chan] for mv in vmarr for eps in epsarr for chan in chans]
    for marr in massarr:
        d={"mv" : marr[0],"mdm" : marr[1], "eps" : marr[2], "channels" : [marr[3]], "signal_chan" : "NCE_nucleon", "output_mode" : "comprehensive", "det_switch" : "superk", "samplesize" : 500, "model" : "Dark_Photon", "sumlog" : "Events/t2k_superk_chan2.dat", "ptmax" : 2}
        t2k_eval(d)

def execute_coherent(genlist=True):
    if genlist:
        d={"prod_chan" : ["pi0_decay"], "proddist" : ["burmansmith"], "samplesize" : 1e6, "output_mode" : "particle_list", "partlistfile" : ["data/particle_list_coherent.dat"], "p_num_target" : 80}
        write_coherent(d=d)
        subp.call(["./build/main","parameter_run.dat"])
    #vmassarr=[i for i in range(11,30,2)]+[i for i in range(30,130,10)]+[129,131,132,134,136,138,140,145,150,155,160,165,170,180,190,195,200]+[3,5,6,9]
    #vmassarr=[9]
    vmassarr=[1,2,5,10,15,17,19,21,23,25,30,40,50,60,70,80,90,95,100,105,110,115,120,125,130,132,134,135,137,140,150,160,180,200,220,250,300,350,400,500,700,900,1000,1500,2000,3000,5000,7500,10000,100000]
    #vmassarr=[10]
    massarr=[[MV,10] for MV in vmassarr]
    #massarr=[[MV,MV/3.0] for MV in vmassarr]
    for marr in massarr:
        #d={"mv" : marr[0], "alpha_D" : 0.5, "mdm" : marr[1], "channels" : [_pion_decay], "signal_chan" : "NCE_nucleon", "det_switch" : "LAr29", "samplesize" : 1000, "sumlog" : "IDM_Events/coherent_Ar_29kg_10mev.dat", "outlog" : "Events_coherent/coherent_Ar_{}_{}.dat".format(str(marr[0]),str(marr[1])), "efficiency" : 1, "min_scatter_energy" : 2e-5, "max_scatter_energy" : 0.05, "burn_max" : 100, "output_mode" : "summary", "coherent" : "true", "POT" : 4.2e22, "model" : "Inelastic_Dark_Matter", "kinetic_energy_cut" : "true"}
        #d={"mv" : marr[0], "alpha_D" : 0.5, "mdm" : marr[1], "channels" : [_pion_decay], "signal_chan" : "NCE_nucleon", "det_switch" : "LAr29", "samplesize" : 2000, "sumlog" : "Events/coherent_Ar_29kg.dat", "outlog" : "Events_coherent/coherent_Ar_{}_{}.dat".format(str(marr[0]),str(marr[1])), "efficiency" : 1, "min_scatter_energy" : 8e-5, "max_scatter_energy" : 0.05, "burn_max" : 1000, "output_mode" : "summary", "coherent" : "true", "POT" : 4.2e22, "model" : "Dark_Photon"}
        #coherent_eval(d)
        #d.update({"model" : "Inelastic_Dark_Matter", "sumlog" : "IDM_Events/coherent_Ar_29kg.dat"})
        #d={"mv" : marr[0],  "mdm" : marr[1], "channels" : [_pion_decay], "signal_chan" : "NCE_nucleon_baryonic", "det_switch" : "csi1T", "samplesize" : 500, "sumlog" : "Events/coherent_CsI_1T.dat"}
        d={"mv" : marr[0], "alpha_D" : 0.5, "mdm" : marr[1], "channels" : [_pion_decay], "signal_chan" : "NCE_nucleon", "det_switch" : "LAr29", "samplesize" : 1000, "sumlog" : "IDM_Events/coherent_Ar_29kg.dat", "outlog" : "Events_coherent/coherent_Ar_{}_{}.dat".format(str(marr[0]),str(marr[1])), "efficiency" : 1, "min_scatter_energy" : 8e-5, "max_scatter_energy" : 0.05, "burn_max" : 300, "output_mode" : "summary", "coherent" : "true", "POT" : 4.2e22, "model" : "Inelastic_Dark_Matter", "kinetic_energy_cut" : "true"}
        coherent_eval(d)
        d.update({"alpha_D" : "0.1"})
        coherent_eval(d)


def execute_lsnd(genlist=True):
    if genlist:
        d={"prod_chan" : ["pi0_decay"], "proddist" : ["burmansmith"], "samplesize" : 0.55e6, "output_mode" : "particle_list", "partlistfile" : ["data/particle_list_lsnd_low.dat"], "p_num_target" : 6}
        write_lsnd(d=d)
        subp.call(["./build/main","parameter_run.dat"])
        d={"prod_chan" : ["pi0_decay"], "proddist" : ["burmansmith"], "samplesize" : 0.45e6, "output_mode" : "particle_list", "partlistfile" : ["data/particle_list_lsnd_high.dat"], "p_num_target" : 70}
        write_lsnd(d=d)
        subp.call(["./build/main","parameter_run.dat"])
        arr1 = np.loadtxt("data/particle_list_lsnd_low.dat")
        arr2 = np.loadtxt("data/particle_list_lsnd_high.dat")
        arr3 = np.append(arr1,arr2,axis=0)
        np.random.shuffle(arr3)
        np.savetxt("data/particle_list_lsnd.dat",arr3)
    #vmassarr = [1,2,3,4,5,7,9,10,12,15,17]+[x for x in range(20,135,5)]+[131,132,133,134,135,136,137,139]+[x for x in range(140,200,10)]+[x for x in range(200,700,25)]
    #dmarr =[1,2,3,4,5,7,9,10,12,15,17]+[x for x in range(20,135,5)]+[64,66]
    #vmarr = [1,2,3,4,5,6,7,8,9,10]
    #dmarr = [0.5,1,1.5,2,2.5,3,3.5,4,4.5,5]
    #massarr = [[v,x] for v in vmarr for x in dmarr if x<=v/2.0 and x<134/2.0]
    vmarr=[3,90]
    massarr = [[v,v/3.0] for v in vmarr]
    d_list=[]
    #massarr = [[4,1]]
    for marr in massarr:
        #d={"mv" : marr[0], "alpha_D" : 0.1, "mdm" : marr[1], "channels" : [_pion_decay], "signal_chan" : "NCE_electron", "det_switch" : "lsnd", "samplesize" : 10000, "sumlog" : "Claudia/lsnd.dat", 'outlog' : 'Claudia/lsnd_dm_{}_{}.dat'.format(marr[0],marr[1]), 'output_mode' : 'dm_detector_distribution'}
        #lsnd_eval(d)
        d={"mv" : marr[0], "alpha_D" : 0.5, "eps" : 1, "mdm" : marr[1], "channels" : [_pion_decay], "signal_chan" : "NCE_electron", "det_switch" : "lsnd", "samplesize" : 10000, "sumlog" : "lsnd.dat", 'outlog' : 'lsnd_elec_events.dat'.format(marr[0],marr[1]), 'output_mode' : 'summary', "max_scatter_energy" : 0.0528}
        lsnd_eval(d)
        #d_list.append(copy.deepcopy(d))
    #pool = Pool(processes=4)
    #pool.map(lsnd_eval,d_list)

def execute_pip2(genlist=True):
    if genlist:
        d={"prod_chan" : ["pi0_decay"], "proddist" : ["burmansmith"], "samplesize" : 1e6, "output_mode" : "particle_list", "partlistfile" : ["data/particle_list_pip2.dat"], "p_num_target" : 74}
        write_pip2(d=d)
        subp.call(["./build/main","parameter_run.dat"])
    vmassarr=[i for i in range(11,30,2)]+[i for i in range(30,130,10)]+[129,131,132,134,136,138,140,145,150,155,160,165,170,180,190,195,200]+[3,5,6,9]
    massarr=[[MV,MV/3.0] for MV in vmassarr]
    for marr in massarr:
        pip2_angle=0
        #d={"POT" : 2e22, "mv" : marr[0], "alpha_D" : 0.5, "mdm" : marr[1], "channels" : [_pion_decay], "signal_chan" : "NCE_nucleon", "det_switch" : "pip2", "samplesize" : 3000, "sumlog" : "Events/pip2_0.dat", "model" : "Dark_Photon", "output_mode" : "summary", "kinetic_energy_cut" : "true", "coherent" : "true", "min_scatter_energy" : 1e-5, "max_scatter_energy" : 0.05, "burn_max" : 1000, "max_trials" : -1}
        d={"POT" : 2e22, "mv" : marr[0], "alpha_D" : 0.5, "mdm" : marr[1], "channels" : [_pion_decay], "signal_chan" : "Pion_Inelastic", "det_switch" : "pip2", "samplesize" : 3000, "sumlog" : "Events/pip2_0_pi0.dat", "model" : "Dark_Photon", "output_mode" : "summary", "kinetic_energy_cut" : "false", "coherent" : "false", "min_scatter_energy" : 0, "max_scatter_energy" : 1, "burn_max" : 1000, "max_trials" : -1}
        pip2_eval(d)
        d.update({"sumlog" : "Events/pip2_30_pi0.dat", "det_switch" : "pip2_30"})
        pip2_eval(d)
        d.update({"sumlog" : "Events/pip2_90_pi0.dat","det_switch" : "pip2_90"})
        pip2_eval(d)

def execute_lanl(genlist=True):
    #if genlist:
        #d={"prod_chan" : ["pi0_decay"], "proddist" : ["burmansmith"], "samplesize" : 1e6, "output_mode" : "particle_list", "partlistfile" : ["data/particle_list_lanl.dat"], "p_num_target" : 74}
        #write_lanl(d=d)
        #subp.call(["./build/main","parameter_run.dat"])
    #vmassarr=[i for i in range(11,30,2)]+[i for i in range(30,130,10)]+[129,131,132,134,136,138,140,145,150,155,160]+[3,5,6,9]
    #vmassarr=[i for i in range(11,30,2)]+[i for i in range(30,130,10)]+[129,131,132,134,136,138,140,145,150,155,160,165,170,180,190,195,200]+[3,5,6,9]
    #vmassarr=[2,2.5,3,4,5,10,15,17,19,19.5,20.5,21,23,25,30,40,50,60,70,80,90,95,100,105,110,115,120,125,130,132,134,135,137,140,150,160,180,200,220,250,300,350,400,500,700,900,1000,1500,2000,3000,5000,7500,10000,100000]
    #massarr=[[MV,MX] for MV in vmassarr for MX in chimassarr]
    #massarr=[[MV,MV/3.0] for MV in vmassarr]
    vmassarr=[5,10,20,50]
    massarr=[[MV,MV] for MV in vmassarr]
    channel="NCE_nucleon"
    for marr in massarr:
        #arr1 = np.loadtxt("data/particle_list_lujan.dat")
        #np.random.shuffle(arr1)
        #np.savetxt("data/particle_list_lujan.dat",arr1)
        #d={"POT" : 2e22,"mv" : marr[0], "alpha_D" : 0.5, "mdm" : marr[1], "channels" : [_pion_decay], "signal_chan" : channel, "det_switch" : "lanl", "samplesize" : 2000, "sumlog" : "Events/CCM_high.dat","model" : "Dark_Photon","output_mode" : "summary", "kinetic_energy_cut" : "true", "coherent" : "true", "min_scatter_energy" : 5e-5, "max_scatter_energy" : 0.05, "burn_max" : 1000, "output_mode" : "summary", "pi0_per_POT" : 0.1145}
        #lanl_eval(d)
        #d={"POT" : 2e22,"mv" : marr[0], "alpha_D" : 0.5, "mdm" : marr[1], "channels" : [_pion_decay], "signal_chan" : channel, "det_switch" : "lanl", "samplesize" : 2000, "sumlog" : "Events/CCM.dat","model" : "Dark_Photon","output_mode" : "summary", "kinetic_energy_cut" : "true", "coherent" : "true", "min_scatter_energy" : 2e-5, "max_scatter_energy" : 0.05, "burn_max" : 1000, "output_mode" : "summary", "pi0_per_POT" : 0.1145}
        d={"POT" : 1e21,"mv" : marr[0], "alpha_D" : 0, "mdm" : marr[1], "channels" : [_pion_decay], "signal_chan" : channel, "det_switch" : "lanl", "samplesize" : 20000, "sumlog" : "Events/CCM_DP.dat","outlog" : "Events/CCM_DP_{}mev.dat".format(str(marr[0])),"model" : "Dark_Photon","output_mode" : "dm_detector_distribution", "kinetic_energy_cut" : "false", "coherent" : "false", "min_scatter_energy" : 0, "max_scatter_energy" : 10, "burn_max" : 1000, "max_trials" : -1, "signal_chan" : "Signal_Decay"}
        lanl_eval(d)
        #d.update({"alpha_D" : 0.1})
        #lanl_eval(d)
        #d={"mv" : marr[0], "alpha_D" : 1e-3, "mdm" : marr[1], "channels" : [_pion_decay], "signal_chan" : channelb, "det_switch" : "lanl", "samplesize" : 2000, "sumlog" : "Events/lanl20_b.dat"}
        #lanl_eval(d)
        #d={"mv" : marr[0], "alpha_D" : 1e-3, "mdm" : marr[1], "channels" : [_pion_decay], "signal_chan" : channelb, "det_switch" : "lanl_far", "samplesize" : 2000, "sumlog" : "Events/lanl40_b.dat"}
        #lanl_eval(d)

def execute_nucal(genlist=True):
    d_list=[]
    if genlist:
        d={"prod_chan" : ["pi0_decay"],"proddist" : ["bmpt"],"samplesize" : 2e6,"output_mode" : "particle_list","partlistfile" : ["data/particle_list_nucal.dat"]}
        write_nucal(d=d)
        subp.call(["./build/main", "parameter_run.dat"])
    #vmarr=[1.15,1,1.2,1.3,1.5,2,3,4,5,10,15,20,30,40,60,80,100,130,140,150,200,250,300,400,500,540,550,560,600,700,725,750,760,765,767,770,772,775,780,790,800,900,1000,1100,1200,1300,1400,1500,1750,2000,2250,2500,2750,3000,3500,4000,4500,5000]
    vmarr=[10,15,20,25,30,35,40,45,50,55,60,65,70,75,80,90,100,115,130,140,150,200,250,300,400,500,540,550,560,600,700,725,750,760,765,767,770,772,775,780,790,800,900,1000,1100,1200,1300,1400,1500,1750,2000,2100,2250,2400,2500]
    #,2750,3000,3500,4000,4500,5000,5500,6000]
    #vmarr=[1712,1714,1716,1718]
    #vmarr=[71,72,73,74,4550,4600,4650,4700,4750,4800,4850,4900,4950]
    #vmarr=[4500+x for x in range(50,500,50)]
    #vmarr=[1.1,1.075,1.05,1.04,1.03,1.15,1.2,1.3,1.4,1.5,1.7,1.9,2,3,4,5,10,15,20,30,40,60,80,100,130,140,150, 175,200,225,250,275,300,325,350,375,400,405,410,415,420,700,770,800,405,410,415,420,460,470,480,490,500,550,600, 650,700,770,800,725,750,760,765,767,772,775,790,810,820,830,840,850,860,870,880,890,900,910,920,930,940,950,960, 970,980,990,1000]
    #massarr=[[mv,mv/3.0,1.05*mv/3.0] for mv in vmarr]
    massarr=[[mv,mv/2.5,1.4*mv/2.5] for mv in vmarr]
    #massarr=[[mv,mv,1e-7] for mv in vmarr]
    #d=({"signal_chan" : "Signal_Decay", "output_mode" : "comprehensive", "samplesize" : 50000, "sumlog" : "Visible_Dark_Photon/nucal_events.dat", "model" : "Dark_Photon", "min_scatter_energy" : 3, "max_scatter_energy" : 70, "min_scatter_angle" : 0, "max_scatter_angle" : 6, "det_switch" : "nucal", "weighted" : "true"});
    d=({"channels" : [_pion_decay,_eta_decay,_brem], "signal_chan" : "DM2_Signal_Decay", "output_mode" : "comprehensive", "samplesize" : 50000,  "sumlog" : "IDM_Events2.5/nucal_decay.dat", "model" : "Inelastic_Dark_Matter", "min_scatter_energy" : 3, "max_scatter_energy" : 70, "min_scatter_angle" : 0, "max_scatter_angle" : 6, "det_switch" : "nucal", "eps" : 1e-4, "weighted" : 'true'});
    for marr in massarr:
        #d.update({"mv" : marr[0],"mdm1" : marr[1], "mdm2" : marr[2], "outlog" : "IDM_Events/nucal_aD0.5_Delta_0.05/nucal_{}_{}.dat".format(marr[0],marr[2]), "alpha_D" : 0.5})
        #nucal_eval(d)
        #d.update({"mv" : marr[0],"mdm1" : marr[1], "mdm2" : marr[1]*1.1, "outlog" : "IDM_Events/nucal_aD0.1_Delta_0.1/nucal_{}_{}.dat".format(marr[0],marr[2]), "alpha_D" : 0.1})
        #nucal_eval(d)
        d.update({"mv" : marr[0],"mdm1" : marr[1], "mdm2" : marr[2], "outlog" : "IDM_Events2.5/nucal_aD0.1_Delta0.4/nucal_{}_{}.dat".format(marr[0],marr[2]),"alpha_D" : 0.1})
        #d.update({"mv" : marr[0],"mdm" : marr[1], "eps" : marr[2], "outlog" : "Visible_Dark_Photon/nucal/nucal_decays_{}_{}.dat".format(marr[0],marr[2])})
        nucal_eval(d)
    #massarr=[[mv,mv/2.5,1.4*mv/2.5] for mv in vmarr]
    #for marr in massarr:
    #    d.update({"mv" : marr[0],"mdm1" : marr[1], "mdm2" : marr[2], "outlog" : "IDM_Events2.5/nucal_aD0.1_Delta0.4/nucal_{}_{}.dat".format(marr[0],marr[2]),"alpha_D" : 0.1})
    #    nucal_eval(d)

def execute_charm(genlist=True):
    if genlist:
        d={"prod_chan" : ["pi0_decay"],"proddist" : ["bmpt"],"samplesize" : 2e6,"output_mode" : "particle_list","partlistfile" : ["data/particle_list_charm.dat"]}
        write_charm_decay(d=d)
        subp.call(["./build/main", "parameter_run.dat"])
    vmarr=[1,5,8,10,15,20,30,40,50,60,65,70,75,80,85,90,95,100,115,130,140,150,200,250,300,400,500,540,550,560,600,700,725,750,760,765,767,770,772,775,780,790,800,900,1000,1100,1200,1300,1400,1500,1750,2000,2250,2500,2750,3000,3500,4000]+[x for x in range(4100,6500,100)]
    #vmarr=[2,3,4,5,7,9,10,15,20,30,40,60,80,100,130,140,150,200,250,300,350,400,450,500,540,550,560,600,650,700,725,750,760,765,767,770,772,775,780,790,800,900,1000,1100,1200,1300,1400,1500,1750,2000,2250,2500,2750,3000,3500,4000]
    vmarr=[2050,2100,2150,2200]
    #epsarr=[2e-7]
    #massarr=[[mv,mv,eps] for mv in vmarr for eps in epsarr]
    #massarr=[[mv,mv/3.0,1.05*mv/3.0] for mv in vmarr]
    massarr=[[mv,mv/2.5,1.4*mv/2.5] for mv in vmarr]    
    #d=({"channels" : [_pion_decay,_eta_decay,_brem], "det_switch" : "charm", "signal_chan" : "Signal_Decay", "output_mode" : "dm_detector_distribution", "samplesize" : 50000, "min_scatter_energy" : 1, "max_scatter_energy" : 1e5, "efficiency" : 1, "alpha_D" : 0.5, "POT" : 2.4e18, "sumlog" : "YuDai_project/charm/charm_events.dat"});
    #d=({"channels" : [_pion_decay,_eta_decay,_brem], "det_switch" : "charm", "signal_chan" : "Signal_Decay", "output_mode" : "comprehensive", "samplesize" : 50000, "min_scatter_energy" : 3, "max_scatter_energy" : 1e5, "efficiency" : 1, "alpha_D" : 0.5, "POT" : 2.4e18, "sumlog" : "YuDai_project/charm/charm_decay_events.dat"});
    #d=({"channels" : [_pion_decayx,_eta_decay,_brem], "det_switch" : "charm", "signal_chan" : "Signal_Decay", "output_mode" : "comprehensive", "samplesize" :    50000, "min_scatter_energy" : 3, "max_scatter_energy" : 1e5, "efficiency" : 1,     "alpha_D" : 0.5, "POT" : 2.4e18, "sumlog" : "Visible_Dark_Photon/charm_decay_events.dat", "weighted" : "true"});
    d=({"channels" : [_pion_decay,_eta_decay,_brem], "det_switch" : "charm", "signal_chan" : "DM2_Signal_Decay", "output_mode" : "comprehensive", "samplesize" : 20000, "min_scatter_energy" : 3, "max_scatter_energy" : 1e5, "efficiency" : 1, "alpha_D" : 0.1, "POT" : 2.4e18, "sumlog" : "IDM_Events2.5/charm_decay_events.dat", "eps" : 1e-4, "model" : "Inelastic_Dark_Matter", "weighted" : 'true'})
    #d=({"channels" : [_pion_decay,_eta_decay,_brem], "det_switch" : "charm", "signal_chan" : "DM2_Signal_Decay", "output_mode" : "dm_detector_distribution", "samplesize" : 10000, "min_scatter_energy" : 3, "max_scatter_energy" : 1e5, "efficiency" : 1, "alpha_D" : 0.5, "POT" : 2.4e18, "sumlog" : "IDM_Events/charm_decay_events.dat", "eps" : 1e-3, "model" : "Inelastic_Dark_Matter", "weighted" : 'true'});
    for marr in massarr:
        d.update({"mv" : marr[0], "mdm1" : marr[1], "mdm2" : marr[2]})
        #d.update({"mv" : marr[0], "mdm" : marr[1], "eps" : marr[2]})
        #d.update({"outlog" : "Visible_Dark_Photon/CHARM/charm_decays_{}_{}.dat".format(marr[0],marr[2])})
        d.update({"outlog" : "IDM_Events2.5/charm_aD0.1_Delta0.4/charm_decays_{}_{}.dat".format(marr[0],marr[2])})
        #d.update({"outlog" : "IDM_Events/charm_aD0.5_Delta0.1/charm_decays_{}_{}.dat".format(marr[0],marr[2])})
        charm_eval(d)

def execute_charm2(genlist=True):
    if genlist:
        d={"prod_chan" : ["pi0_decay"],"proddist" : ["bmpt"],"samplesize" : 2e6,"output_mode" : "particle_list","partlistfile" : ["data/particle_list_charm2.dat"]}
        write_charm2(d=d)
        subp.call(["./build/main", "parameter_run.dat"])
    vmarr=[1,5,10,15,20,30,40,60,80,100,130,140,150,200,250,300,400,500,540,550,560,600,700,725,750,760,765,767,770,772,775,780,790,800,900,1000,1100,1200,1300,1400,1500]
    #massarr=[[mv,mv,1e-7] for mv in vmarr]
    #vmarr=[30]
    massarr=[[mv,mv/3.0,1e-3] for mv in vmarr]
    #d=({"channels" : [_pion_decay,_eta_decay,_brem], "det_switch" : "charm2", "signal_chan" : "Signal_Decay", "output_mode" : "dm_detector_distribution", "samplesize" : 50000, "min_scatter_energy" : 0, "max_scatter_energy" : 1e5, "efficiency" : 1, "alpha_D" : 0.5, "POT" : 0.25e20, "sumlog" : "YuDai_project/charm2/charm2events.dat"});
    d=({"channels" : [_pion_decay,_eta_decay,_brem,_parton], "signal_chan" : "NCE_electron", "det_switch" : "charm2", "output_mode" : "comprehensive", "samplesize" : 20000, "alpha_D" : 0.5, "sumlog" : "Claudia/charm2.dat", "efficiency" : 1, "model" : "Dark_Photon"})
    for marr in massarr:
        d.update({"mv" : marr[0],"mdm" : marr[1], "eps" : marr[2]})
        d.update({"outlog" : "Claudia/charm2_{}_{}.dat".format(marr[0],round(marr[1],3))})
        #d.update({"outlog" : "YuDai_project/charm2/charm2_dp_{}_{}.dat".format(marr[0],marr[2])})
        charm2_eval(d)

def execute_ship(genlist=True):
    if genlist:
        d={"prod_chan" : ["pi0_decay"],"proddist" : ["bmpt"],"samplesize" : 2e6,"output_mode" : "particle_list","partlistfile" : ["data/particle_list_ship.dat"]}
        write_ship(d=d)
        subp.call(["./build/main", "parameter_run.dat"])
    #vmarr=[1,5,10,15,20,30,40,60,80,100,130,140,150,200,250,300,400,500,540,550,560,600,700,725,750,760,765,767,770,772,775,780,790,800,900,1000,1100,1200,1300,1400,1500]
    #vmarr=[7,8,9,10,15,20,30,40,50,60,65,70,75,80,90,100,115,130,140,150,200,250,300,400,500,540,550,560,600,700,725,750,760,765,767,770,772,775,780,790,800,900,1000,1100,1200,1300,1400,1500,1750,2000,2250,2500,2750,3000,3500,4000,4500,5000,5500,6000,6500,7000]
    vmarr=[2800,2850,2900,2950]
    #vmarr=[1.2,1.3,1.5,1.7,2,2.5,3,3.5,4,5,7,10,15,20]+[x for x in range(30,150,10)]+[x for x in range(200,1400,50)]+[134,541,760,765,767,770,772,775,780,790]
    #vmarr=[x/10.0 for x in range (26,30)]
    #vmarr+=[x for x in range(4600,7500,100)]
    #vmarr=[51,52,53,54]
    #massarr=[[mv,mv,1e-6] for mv in vmarr]
    #vmarr=[6600,6700,6800,6900]
    mass_arr=[[mv,mv/2.5,1.4*mv/2.5] for mv in vmarr]
    #mass_arr=[[mv,mv/3.0,1.05*mv/3.0] for mv in vmarr]
    #d=({"signal_chan" : "NCE_nucleon", "output_mode" : "dm_detector_distribution", "samplesize" : 50000, "min_scatter_energy" : 0, "max_scatter_energy" : 1e5, "efficiency" : 0.5, "alpha_D" : 0.5, "POT" : 6e20});
    #d=({"channels" : [_pion_decay,_eta_decay,_brem,_parton] ,"signal_chan" : "NCE_electron", "output_mode" : "comprehensive", "samplesize" : 10000, "alpha_D" : 0.5, "sumlog" : "Claudia/ship.dat", "model" : "Dark_Photon"})
    d=({"channels" : [_pion_decay,_eta_decay,_brem], "signal_chan" : "DM2_Signal_Decay", "output_mode" : "comprehensive", "samplesize" : 20000,  "sumlog" : "IDM_Events2.5/na62_decay.dat", "min_scatter_energy" : 3, "max_scatter_energy" : 400, "min_scatter_angle" : 0, "max_scatter_angle" : 6, "eps" : 1e-3, "weighted" : 'true', 'model' : 'Inelastic_Dark_Matter', 'efficiency' : 1, 'det_switch' : "na62", 'POT' : 1e18});
    #d=({"channels" : [_pion_decay,_eta_decay,_brem], "signal_chan" : "Signal_Decay", "output_mode" : "comprehensive", "samplesize" : 100000,  "sumlog" : "Visible_Dark_Photon/na62_decay.dat", "model" : "Dark_Photon", "min_scatter_energy" : 3, "max_scatter_energy" : 400, "min_scatter_angle" : 0, "max_scatter_angle" : 6, "eps" : 1e-3, "weighted" : 'true', 'efficiency' : 1, 'det_switch' : "na62", 'POT' : 1e18});
    #d=({"channels" : [_pion_decay,_eta_decay,_brem,_parton], "signal_chan" : "NCE_electron", "det_switch" : "charm2", "output_mode" : "comprehensive", "samplesize" : 10000, "alpha_D" : 0.5, "sumlog" : "Claudia/charm2.dat", "efficiency" : 1, "model" : "Dark_Photon"})
    for marr in mass_arr:
        d.update({"mv" : marr[0],"mdm1" : marr[1], "mdm2" : marr[2], "outlog" : "IDM_Events2.5/na62_aD0.1_Delta0.4/na62_{}_{}.dat".format(marr[0],marr[2]), "alpha_D" : 0.1})
        ship_eval(d)
        #d.update({"mv" : marr[0],"mdm" : marr[0], "eps" : marr[2], "outlog" : "Visible_Dark_Photon/NA62/na62_{}_{}.dat".format(marr[0],marr[2]), "alpha_D" : 0})
        #ship_eval(d)
        #d.update({"mv" : marr[0],"mdm1" : marr[1], "mdm2" : marr[1]*1.4, "outlog" : "IDM_Events/na62_aD0.1_Delta0.4/na62_{}_{}.dat".format(marr[0],marr[2]),"alpha_D" : 0.1})
        #ship_eval(d)
