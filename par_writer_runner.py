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
    elif model == "Dark_Photon" or model == "Axion_Dark_Photon":
        if ((MX/1000.0<mpi0/2.0 and MV<600.0) or (signal_channel=="Signal_Decay" and MV/1000.0<mpi0)) and _pion_decay in channels:
            proddist.append("particle_list")
            prodchan.append("pi0_decay")
            partlistfile.append("data/particle_list.dat")
            executing=True
        if ((MX/1000.0<meta/2.0 and MV<900.0) or (signal_channel=="Signal_Decay" and MV/1000.0<meta)) and _eta_decay in channels:
            proddist.append("particle_list")
            prodchan.append("eta_decay")
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
        write_sbnd(d=d)
    elif det_switch == "miniboone_full":
        write_miniboone(d=d,det=miniboone_detector_full)
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

ship_defaults = {"mv" : 30, "mdm" : 10, 'channels' : {_parton,_brem,_pion_decay,_eta_decay,_phi_decay}, "signal_chan" : "NCE_electron", 'det_switch' : 'ship', 'sumlog' : "Events/ship.dat"}

def ship_eval(d_user):
    d = ship_defaults.copy()
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

    BEAM_ENERGY=400
    zmin=max(3*MV/1000.0/BEAM_ENERGY,0.1)
    zmax=1-zmin

    proddist = []
    prodchan = []
    partlistfile = [    ]
    executing=False

    if signal_channel=="Inelastic_Nucleon_Scattering_Baryonic" or signal_channel=="Baryonic_Test":
        print("Need to implement Baryonic!")
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

    d.update({"proddist" : proddist, "prod_chan" : prodchan, "partlistfile" : partlistfile,"mv" : MV/1000.0, "mdm" : MX/1000.0, "zmin" : zmin, "zmax" : zmax, "outfile" : outfile})
    if det_switch == "ship":
        write_ship(d=d)
    elif det_switch == "test_sphere":
        write_ship(d=d,det=test_sphere)
    elif det_switch == "test_cylinder":
        write_ship(d=d,det=test_cylinder)
    elif det_switch == "test_cuboid":
        write_ship(d=d,det=test_cuboid)
    subp.call(["./build/main", outfile])
    t1 = time.time()
    print("\ntime={}\n".format(t1-t0))
    t0 = time.time()
    subp.call(["rm", outfile])


numi_defaults = {"mv" : 30, "mdm" : 10, 'channels' : {_parton,_brem,_pion_decay,_eta_decay,_phi_decay}, "signal_chan" : "Signal_Decay", 'det_switch' : 'nova', 'sumlog' : "Events/nova_decay.dat", "alpha_D" : 0.5}

def numi_eval(d_user):
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
            partlistfile.append("data/particle_list_numi.dat")
            executing=True
        if MX/1000.0<meta/2.0 and MV<1000 and _eta_decay in channels:
            proddist.append("particle_list")
            prodchan.append("eta_decay_baryonic")
            partlistfile.append("data/particle_list_numi.dat")
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
    else:
        if ((MX/1000.0<mpi0/2.0 and MV<600.0) or (MV/1000.0<mpi0 and signal_channel=="Signal_Decay")) and _pion_decay in channels:
            proddist.append("particle_list")
            prodchan.append("pi0_decay")
            partlistfile.append("data/particle_list_numi.dat")
            executing=True
        if ((MX/1000.0<meta/2.0 and MV<900.0) or (signal_channel=="Signal_Decay" and MV/1000.0<meta)) and _eta_decay in channels:
            proddist.append("particle_list")
            prodchan.append("eta_decay")
            partlistfile.append("data/particle_list_numi.dat")
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
            partlistfile.append("data/particle_list_numi.dat")
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
    d.update({"proddist" : proddist, "prod_chan" : prodchan, "partlistfile" : partlistfile,"mv" : MV/1000.0, "mdm" : MX/1000.0, "zmin" : zmin, "zmax" : zmax, "outfile" : outfile})
    if det_switch == "nova":
        write_numi(d=d)
    elif det_switch == "nova_absorber":
        write_numi(d=d,det=NOvA_absorber_detector)
    elif det_switch == "minos":
        print("EXECUTING MINOS")
        write_numi(d=d,det=MINOS_detector)
    elif det_switch == "minos_absorber":
        write_numi(d=d,det=MINOS_absorber_detector)
    elif det_switch == "miniboone_numi":
        write_numi(d=d,det=miniboone_detector_numi)
    elif det_switch == "test_sphere":
        write_numi(d=d,det=test_sphere)
    elif det_switch == "test_cylinder":
        write_numi(d=d,det=test_cylinder)
    elif det_switch == "test_cuboid":
        write_numi(d=d,det=test_cuboid)
    subp.call(["./build/main", outfile])
    t1 = time.time()
    print("\ntime={}\n".format(t1-t0))
    t0 = time.time()
    subp.call(["rm", outfile])
#####END OF numi_eval############

coherent_defaults = {"mv" : 30, "mdm" : 10, 'channels' : {_pion_decay,_piminus_cap}, "signal_chan" : "NCE_Nucleon", 'det_switch' : 'csi', 'sumlog' : "Events/coherent.dat", "alpha_D" : 0.5}

def coherent_eval(d_user):
    d = coherent_defaults.copy()
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
        if MX/1000.0<mpi0/2.0 and MV<600.0 and _pion_decay in channels:
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
    if det_switch == 'csi1T':
        write_coherent(d=d,det=coherent_detector_CsI_1T)
    if det_switch == "captain":
        write_coherent(d=d,det=captain_detector)
    if det_switch == "captain_off_axis":
        write_coherent(d=d,det=captain_detector_off)
    else:
        write_coherent(d=d)

    subp.call(["./build/main", outfile])
    t1 = time.time()
    print("\ntime={}\n".format(t1-t0))
    t0 = time.time()
    subp.call(["rm", outfile])
#END OF coherent_eval

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
    executing=False
    if signal_channel=="NCE_nucleon":
        if MX/1000.0<mpi0/2.0 and MV<600.0 and _pion_decay in channels:
            proddist.append("particle_list")
            prodchan.append("pi0_decay")
            partlistfile.append("data/particle_list_lanl.dat")
            executing = True
        if MX/1000.0<0.129/2.0 and MV<600 and MV>2*MX and _piminus_cap in channels:
            proddist.append("")
            prodchan.append("piminus_capture")
            partlistfile.append("")
            executing = True
    elif signal_channel=="NCE_nucleon_baryonic":
        if MX/1000.0<mpi0/2.0 and MV<600.0 and _pion_decay in channels:
            proddist.append("particle_list")
            prodchan.append("pi0_decay_baryonic")
            partlistfile.append("data/particle_list_lanl.dat")
            executing = True
    if not executing:
        return
    if signal_channel=="NCE_nucleon":
        if det_switch == "lanl":
            user2 = {"samplesize" : 1000, "POT" : 2.2e21, "min_scatter_energy" : 10e-6, "max_scatter_energy" : 0.2, "efficiency" : 0.8, "sumlog" : "Events/lanl.dat", "coherent" : "true", "eps" : 1e-3, "burn_max" : 1000}
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
    elif det_switch =="lanl_far":
        write_lanl(d=d,det=lanl_detector_far)

    subp.call(["./build/main", outfile])
    t1 = time.time()
    print("\ntime={}\n".format(t1-t0))
    t0 = time.time()
    subp.call(["rm", outfile])

def lsnd_eval(d_user):
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

    zmin =0
    zmax =1

    t0 = time.time()
    outfile="parameter_run_{0}_{1}.dat".format(str(MV),str(MX))
    proddist = []
    prodchan = []
    partlistfile = [    ]
    executing=False
    if signal_channel=="Inelastic_Nucleon_Scattering_Baryonic" or signal_channel=="Baryonic_Test":
        pass
    else:
        print("trying this")
        if ((MX/1000.0<mpi0/2.0 and MV<600.0) or (MV/1000.0<mpi0 and signal_channel=="Signal_Decay")) and _pion_decay in channels:
            proddist.append("particle_list")
            prodchan.append("pi0_decay")
            partlistfile.append("data/particle_list_lsnd.dat")
            executing=True
    if not executing:
        return

    d.update(d_user)
    d.update({"proddist" : proddist, "prod_chan" : prodchan, "partlistfile" : partlistfile,"mv" : MV/1000.0, "mdm" : MX/1000.0, "zmin" : zmin, "zmax" : zmax, "outfile" : outfile})
    if det_switch == "lsnd":
        write_lsnd(d=d)
    subp.call(["./build/main", outfile])
    t1 = time.time()
    print("\ntime={}\n".format(t1-t0))
    t0 = time.time()
    subp.call(["rm", outfile])
#####END OF lsnd_eval############



def execute_miniboone_parallel(genlist = True):
    if genlist:
        d={"prod_chan" : ["pi0_decay"],"proddist" : ["pi0_sanfordwang"],"samplesize" : 2e6,"output_mode" : "particle_list","partlistfile" : ["data/particle_list.dat"]}
        write_miniboone(d=d)
        subp.call(["./build/main", "parameter_run.dat"])
        d={"prod_chan" : ["eta_decay"],"proddist" : ["k0_sanfordwang"],"samplesize" : 2e6,"output_mode" : "particle_list","partlistfile" : ["data/particle_list_k0.dat"]}
        write_miniboone(d=d)
        subp.call(["./build/main", "parameter_run.dat"])
    channs={_brem,_pion_decay,_eta_decay}
    #vmarr={10,100,300,500,770}
    #epsarr={1e-6,3e-7,1e-7,3e-7,1e-8}
    #massarr=[[mv,mv,eps] for mv in vmarr for eps in epsarr]
    #massarr = [[80,10,0.0005],[300,10,0.001],[600,10,0.009],[800,10,0.006],[800,10,0.001]]
    massarr=[[90,30,1e-3],[600,200,1e-3]]
    for marr in massarr:
        #d={"mv" : marr[0],"mdm" : marr[1], "eps" : marr[2], "signal_chan" : "Signal_Decay", "output_mode" : "comprehensive", "det_switch" : "miniboone_full", "sumlog" : "Decay_Events/miniboone_decay.dat", "samplesize" : 1000}
        #d={"mv" : marr[0],"mdm" : marr[1], "eps" : marr[2], "signal_chan" : "NCE_nucleon", "output_mode" : "summary", "det_switch" : "miniboone", "alpha_D" : 0.1, "channels" : channs}
        d={"mv" : marr[0],"mdm" : marr[1], "eps" : marr[2], "signal_chan" : "NCE_electron", "output_mode" : "summary", "det_switch" : "miniboone", "alpha_D" : 0.5, "channels" : channs}
        miniboone_eval(d)

def execute_ship(genlist=True):
    if genlist:
        d={"prod_chan" : ["pi0_decay"],"proddist" : ["bmpt"],"samplesize" : 2e6,"output_mode" : "particle_list","partlistfile" : ["data/particle_list_ship.dat"]}
        write_ship(d=d)
        subp.call(["./build/main", "parameter_run.dat"])
    massarr=[[600,200,1e-3]]
    d=({"signal_chan" : "NCE_nucleon", "output_mode" : "dm_detector_distribution", "samplesize" : 50000, "min_scatter_energy" : 0, "max_scatter_energy" : 1e5, "efficiency" : 0.5, "alpha_D" : 0.5, "POT" : 6e20});
    for marr in massarr:
        d.update({"mv" : marr[0],"mdm" : marr[1], "eps" : marr[2]})
        d.update({"det_switch" : "test_sphere"})
        ship_eval(d)


def execute_numi_test(genlist=True):
    if genlist:
        d={"prod_chan" : ["pi0_decay"],"proddist" : ["bmpt"],"samplesize" : 2e6,"output_mode" : "particle_list","partlistfile" : ["data/particle_list_numi.dat"]}
        write_numi(d=d)
        subp.call(["./build/main", "parameter_run.dat"])
    massarr=[[600,200,1e-3]]
    d=({"signal_chan" : "NCE_nucleon", "output_mode" : "dm_detector_distribution", "samplesize" : 50000, "min_scatter_energy" : 0, "max_scatter_energy" : 1e5, "efficiency" : 0.5, "alpha_D" : 0.5, "POT" : 6e20});
    for marr in massarr:
        d.update({"mv" : marr[0],"mdm" : marr[1], "eps" : marr[2]})
        d.update({"det_switch" : "test_sphere"})
        numi_eval(d)

def execute_numi(genlist=True):
    d_list=[]
    if genlist:
        d={"prod_chan" : ["pi0_decay"],"proddist" : ["bmpt"],"samplesize" : 2e6,"output_mode" : "particle_list","partlistfile" : ["data/particle_list_numi.dat"]}
        write_numi(d=d)
        subp.call(["./build/main", "parameter_run.dat"])
    #vmarr=[1,5,10,15,20,30,40,60,80,100,130,140,150,200,250,300,400,500,540,550,560,600,700,725,750,760,765,767,770,772,775,780,790,800,900,1000,1100,1200,1300,1400,1500]
    #epsarr=[10**n for n in range(-8,-3)]+[3*10**n for n in range(-9,-4)]
    vmarr=[300]
    epsarr=[1e-3]
    massarr=[[mv,mv/3.0,eps] for mv in vmarr for eps in epsarr]
    #massarr=[[mv,10,eps] for mv in vmarr for eps in epsarr]
    #massarr=[[90,30,1e-3],[600,200,1e-3]]
    d=({"signal_chan" : "NCE_nucleon", "output_mode" : "dm_detector_distribution", "samplesize" : 10000, "sumlog" : "Scatter_Events/miniboone_numi_high.dat", "model" : "Dark_Photon","min_scatter_energy" : 0.05, "max_scatter_energy" : 100, "min_scatter_angle" : 0, "ptmax" : 5});
    #d_low=({"signal_chan" : "NCE_electron", "output_mode" : "summary", "samplesize" : 2000,  "sumlog" : "Scatter_Events/miniboone_numi_low.dat", "model" : "Dark_Photon","min_scatter_energy" : 0.2, "max_scatter_energy" : 0.9, "min_scatter_angle" : 0, "ptmax" : 5});
    #d=({"signal_chan" : "NCE_nucleon", "output_mode" : "summary", "samplesize" : 1000, "min_scatter_energy" : 5, "max_scatter_energy" : 35, "efficiency" : 0.5, "alpha_D" : 0.5, "POT" : 6e20});
    d_low=({"signal_chan" : "NCE_nucleon", "output_mode" : "dm_detector_distribution", "samplesize" : 10000, "min_scatter_energy" : 0.05, "max_scatter_energy" : 100, "efficiency" : 0.5, "alpha_D" : 0.5, "POT" : 6e20});
    #d_miniboone=({"signal_chan" : "NCE_electron", "det_switch" : "miniboone_numi", "output_mode" : "summary", "samplesize" : 1000, "min_scatter_energy" : 0.05, "max_scatter_energy" : 10, "efficiency" : 0.35, "alpha_D" : 0.5, "POT" : 6e20,"channels" : [_pion_decay,_eta_decay,_brem,_parton], "sumlog" : "Events/mini_numi_electron.dat", "ptmax" : 5, "burn_max" : 100});
    #d_miniboone_nucleon=({"signal_chan" : "NCE_nucleon", "det_switch" : "miniboone_numi", "output_mode" : "summary", "samplesize" : 1000, "min_scatter_energy" : 0.1, "max_scatter_energy" : 10, "efficiency" : 0.35, "alpha_D" : 0.5, "POT" : 6e20,"channels" : [_pion_decay,_eta_decay,_brem,_parton], "sumlog" : "Events/mini_numi_nucleon.dat","ptmax" : 5, "burn_max" : 100});
    #d = ({"signal_chan" : "NCE_electron", "samplesize" : 1000, "min_scatter_energy" : 0.5, "max_scatter_energy" : 100, "output_mode" : "comprehensive", "efficiency" : 1})
    '''
    for marr in massarr:
        d_miniboone.update({"mv" : marr[0],"mdm" : marr[1], "eps" : marr[2]})
        d_miniboone_nucleon.update({"mv" : marr[0],"mdm" : marr[1], "eps" : marr[2]})
        d_list.append(copy.deepcopy(d_miniboone))
        d_list.append(copy.deepcopy(d_miniboone_nucleon))
    '''
    for marr in massarr:
        d.update({"mv" : marr[0],"mdm" : marr[1], "eps" : marr[2]})
        d_low.update({"mv" : marr[0],"mdm" : marr[1], "eps" : marr[2]})
        d.update({"det_switch" : "miniboone_numi", "channels" : [_pion_decay,_eta_decay,_brem,_parton],  "outlog" : "Scatter_Events/miniboone_numi_nucleon_{}_{}.dat".format(marr[0],marr[1])})
        d_low.update({"det_switch" : "nova", "channels" : [_pion_decay,_eta_decay,_brem,_parton], "outlog" : "Scatter_Events/nova_nucleon_low_{}_{}.dat".format(marr[0],round(marr[1],2))})
        #numi_eval(d_low)
        #d_low.update({"det_switch" : "nova","channels" : [_pion_decay,_eta_decay,_brem,_parton], "sumlog" : "Scatter_Events/nova_electron_low.dat", "outlog" : "Scatter_Events/nova_electron_low_{}_{}.dat".format(marr[0],round(marr[1],2))})
        numi_eval(d_low)
        #d_list.append(copy.deepcopy(d))
        #d_list.append(copy.deepcopy(d_low))
        #d.update({"det_switch" : "nova_absorber","channels" : [_pion_decay,_eta_decay,_brem,_parton], "sumlog" : "Events/nova_electron_abs_high.dat"})
        #d_low.update({"det_switch" : "nova_absorber","channels" : [_pion_decay,_eta_decay,_brem,_parton], "sumlog" : "Events/nova_electron_abs_low.dat"})
        numi_eval(d)
        #d_list.append(copy.deepcopy(d))
        #d_list.append(copy.deepcopy(d_low))
    '''
    #vmarr=[100]
    #epsarr=[2.6*10**-5]
    vmarr=[10,30,60,100,150,200,250,300,400,500,600,700,800,900,1000]
    epsarr=[10**n for n in range(-8,-3)]+[3*10**n for n in range(-9,-4)]
    massarr=[[mv,mv,eps] for mv in vmarr for eps in epsarr]
    #d=({"signal_chan" : "Signal_Decay", "output_mode" : "comprehensive", "samplesize" : 5000, "model" : "Dark_Photon","min_scatter_energy" : 5, "max_scatter_energy" : 35, "min_scatter_angle" : 0});
    d_low=({"signal_chan" : "Signal_Decay", "output_mode" : "comprehensive", "samplesize" : 1000, "model" : "Dark_Photon","min_scatter_energy" : 0.5, "max_scatter_energy" : 5, "min_scatter_angle" : 0});
    for marr in massarr:
        #d.update({"mv" : marr[0],"mdm" : marr[1], "eps" : marr[2]})
        d_low.update({"mv" : marr[0],"mdm" : marr[1], "eps" : marr[2]})
        #d.update({"det_switch" : "nova","channels" : [_pion_decay,_eta_decay,_brem], "sumlog" : "Events/nova_dec_high.dat"})
        #d_low.update({"det_switch" : "nova","channels" : [_pion_decay,_eta_decay,_brem], "sumlog" : "Decay_Events/nova_dec_low.dat", "outlog" : "Decay_Events/nova_dec_low_TEST.eps"})
        #d_list.append(copy.deepcopy(d))
        #d_list.append(copy.deepcopy(d_low))
        #d.update({"det_switch" : "nova_absorber","channels" : [_pion_decay,_eta_decay,_brem], "sumlog" : "Events/nova_dec_abs_high.dat", "outlog" : "Decay_Events/nova_dec_abs_high_{}_{}.dat".format(marr[0],marr[2]) })
        d_low.update({"det_switch" : "nova_absorber","channels" : [_pion_decay,_eta_decay,_brem], "sumlog" : "Events/nova_dec_abs_low.dat", "outlog" : "Decay_Events/nova_dec_abs_low_{}_{}.dat".format(marr[0],marr[2])})
        #d_list.append(copy.deepcopy(d))
        d_list.append(copy.deepcopy(d_low))
    for marr in massarr:
        d.update({"mv" : marr[0],"mdm" : marr[1], "eps" : marr[2]})
        numi_eval(d)
        #d_list.append(copy.deepcopy(d))
        numi_eval(d)

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
    '''
    #pool = Pool(processes=3)
    #pool.map(numi_eval,d_list)

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
    #vmassarr=[i for i in range(11,30,2)]+[i for i in range(30,130,10)]+[129,131,132,134,136,138,140,145,150,155,160]+[3,5,6,9]
    vmassarr=[170,180,190,200,225,250]
    #massarr=[[MV,MX] for MV in vmassarr for MX in chimassarr]
    massarr=[[MV,MV/3.0] for MV in vmassarr]
    for marr in massarr:
        d={"mv" : marr[0], "alpha_D" : 0.5, "mdm" : marr[1], "channels" : [_pion_decay, _piminus_cap], "signal_chan" : "NCE_nucleon", "det_switch" : "csi1T", "samplesize" : 2000, "sumlog" : "Events/coherent_CsI_1T.dat"}
        coherent_eval(d)
        #d={"mv" : marr[0],  "mdm" : marr[1], "channels" : [_pion_decay], "signal_chan" : "NCE_nucleon_baryonic", "det_switch" : "csi1T", "samplesize" : 500, "sumlog" : "Events/coherent_CsI_1T.dat"}
        #coherent_eval(d)

def execute_lsnd(genlist=True):
    if genlist:
        d={"prod_chan" : ["pi0_decay"], "proddist" : ["burmansmith"], "samplesize" : 0.55e6, "output_mode" : "particle_list", "partlistfile" : ["data/particle_list_lsnd_low.dat"], "p_num_target" : 8}
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
    vmarr = [1,2,3,4,5,6,7,8,9,10]
    dmarr = [0.5,1,1.5,2,2.5,3,3.5,4,4.5,5]
    massarr = [[v,x] for v in vmarr for x in dmarr if x<=v/2.0 and x<134/2.0]
    #massarr = [[v,v/3.0] for v in vmarr]
    d_list=[]
    #massarr = [[4,1]]
    for marr in massarr:
        d={"mv" : marr[0], "alpha_D" : 0.1, "mdm" : marr[1], "channels" : [_pion_decay], "signal_chan" : "NCE_electron", "det_switch" : "lsnd", "samplesize" : 1000, "sumlog" : "Events/lsnd_test.dat"}
        lsnd_eval(d)
        #d_list.append(copy.deepcopy(d))
    #pool = Pool(processes=4)
    #pool.map(lsnd_eval,d_list)

def execute_lanl(genlist=True):
    if genlist:
        d={"prod_chan" : ["pi0_decay"], "proddist" : ["burmansmith"], "samplesize" : 1e6, "output_mode" : "particle_list", "partlistfile" : ["data/particle_list_lanl.dat"], "p_num_target" : 74}
        write_coherent(d=d)
        subp.call(["./build/main","parameter_run.dat"])
    vmassarr=[i for i in range(11,30,2)]+[i for i in range(30,130,10)]+[129,131,132,134,136,138,140,145,150,155,160]+[3,5,6,9]
    #vmassarr=[170,180,190,200,225,250]
    #massarr=[[MV,MX] for MV in vmassarr for MX in chimassarr]
    massarr=[[MV,MV/3.0] for MV in vmassarr]
    channelb="NCE_nucleon_baryonic"
    channel="NCE_nucleon"
    for marr in massarr:
        d={"mv" : marr[0], "alpha_D" : 0.5, "mdm" : marr[1], "channels" : [_pion_decay, _piminus_cap], "signal_chan" : channel, "det_switch" : "lanl", "samplesize" : 2000, "sumlog" : "Events/lanl20.dat"}
        lanl_eval(d)
        d={"mv" : marr[0], "alpha_D" : 0.5, "mdm" : marr[1], "channels" : [_pion_decay, _piminus_cap], "signal_chan" : channel, "det_switch" : "lanl_far", "samplesize" : 2000, "sumlog" : "Events/lanl40.dat"}
        lanl_eval(d)
        d={"mv" : marr[0], "alpha_D" : 1e-3, "mdm" : marr[1], "channels" : [_pion_decay], "signal_chan" : channelb, "det_switch" : "lanl", "samplesize" : 2000, "sumlog" : "Events/lanl20_b.dat"}
        lanl_eval(d)
        d={"mv" : marr[0], "alpha_D" : 1e-3, "mdm" : marr[1], "channels" : [_pion_decay], "signal_chan" : channelb, "det_switch" : "lanl_far", "samplesize" : 2000, "sumlog" : "Events/lanl40_b.dat"}
        lanl_eval(d)

#execute_lanl(genlist=True)
#execute_lsnd(genlist=False)
execute_numi(genlist=False)
#execute_ship(genlist=True)
#execute_miniboone_parallel(genlist=True)
#execute_t2k(genlist=False)
#execute_coherent(genlist=False)
