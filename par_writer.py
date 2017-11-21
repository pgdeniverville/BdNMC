import subprocess as subp
from contextlib import contextmanager
import os
import pandas as pd
import numpy as np
from multiprocessing import Pool

import Parsing_inelastic
import math

@contextmanager
def cd(newdir):
    prevdir = os.getcwd()
    os.chdir(os.path.expanduser(newdir))
    try:
        yield
    finally:
        os.chdir(prevdir)

mb=1e-31
m_proton=0.938
mpi0=0.134
meta=0.547862
mrho=0.77549
momega=0.782
mphi=1.020
pi=math.pi

meson_per_pi0_miniboone = {'pi0_decay' : '1.0', 'eta_decay' : str(1.0/30.0), 'rho_decay' : str(1.0/20.0), 'omega_decay' : '0.046', 'phi_decay' : str(1.0/150.0), 'pi0_decay_baryonic' : '1.0', 'eta_decay_baryonic' : str(1.0/30.0), 'rho_decay_baryonic' : str(1.0/20.0), 'omega_decay_baryonic' : '0.046', 'phi_decay_baryonic' : str(1.0/150.0)}

meson_per_pi0_lsnd = {'pi0_decay' : '1.0'}

meson_per_pi0_coherent = {'pi0_decay' : 1.0, 'piminus_capture' : '0.63'}

meson_per_pi0_ship = {'pi0_decay' : '1.0', 'eta_decay' : str(0.078), 'rho_decay' : str(0.11), 'omega_decay' : '0.11', 'phi_decay' : str(0.02)}

Hydrogen_string = "material Hydrogen\nnumber_density 7.26942e22\nproton_number 1\nneutron_number 0\nelectron_number 1\nmass 0.945778\n"

Water_string = "material Oxygen\nnumber_density 3.34184e22\nproton_number 8\nneutron_number 8\nelectron_number 8\nmass 0.94578\nmaterial Hydrogen\nnumber_density 6.68368e22\nproton_number 1\nneutron_number 0\nelectron_number 1\nmass 0.945778\n"

Carbon_string = "material Carbon\nnumber_density 3.63471e22\nproton_number 6\nneutron_number 6\nelectron_number 6\nmass 11.2593\n"

Argon_string = "material Argon\nnumber_density 2.11e22\nproton_number 18\nneutron_number 22\nelectron_number 18\nmass {0}\n".format(str(39.948*0.938))

ND280_string = "material nd280stuff\nnumber_density 3.7e23\nproton_number 1\nneutron_number 1\nelectron_number 1\nmass 0.945778\n"

Sodium_Iodide_string = "material Sodium\nnumber_density 1.58e22\nproton_number 11\n neutron_number 23\nelectron_number 11\nmass 21.61\nmaterial Iodine\nnumber_density 1.58e22\nproton_number 53\nneutron_number 72\nelectron_number 53\nmass 119.03\n"

Cesium_Iodide_string = "material Cesium\nnumber_density 1.04e22\nproton_number 55\n neutron_number 78\nelectron_number 55\nmass 132.9\nmaterial Iodine\nnumber_density 1.04e22\nproton_number 53\nneutron_number 72\nelectron_number 53\nmass 119.03\n"

#Don't know what the atomic makeup of the MINOS detector is. Not using this for event generation, so it should be okay.
MINOS_string = "material Steel\nnumber_density 5e24\nproton_number 1\nneutron_number 1\nelectron_number 1\nmass 0\n"

NOvA_string = "material Liquid_Scintillator\nnumber_density 5.16e22\nproton_number 8\nneutron_number 6\nelectron_number 8\nmass 14.011"

defaults = {"eps" : 1e-3, "mdm" : 0.03, "mv" : 0.1, "alpha_D" : 0.1, "prod_chan" : ["pi0_decay"], "signal_chan" : "NCE_nucleon", "outfile" : "parameter_run.dat", "proddist" : [""], "partlistfile" : ["Source/particle_list.dat"], "sumlog" : "Events/miniboone.dat", "outlog" : "Events/miniboone_events.dat", "output_mode" :"summary", "samplesize" : 5000, "min_scatter_energy" : 0.035, "max_scatter_energy" : 1.0, "dm_energy_resolution" : 0.01, "efficiency" : 0.35, "beam_energy" : 8.9, "n_num_target" :
        4, "p_num_target" : 4, "max_trials" : 80e6, "ptmax" : 0.2, "zmin" : 0.3, "zmax" : 0.7, "run" : -1, "POT" : 2e20, "pi0_per_POT" : 0.9, "p_cross" : 25*mb, "meson_per_pi0" : meson_per_pi0_miniboone, "min_scatter_angle" : 0.0, "max_scatter_angle" : 2.1*pi, "repeat" : 1, "timing" : 0.0, "burn_max" : -1,"inelastic_dist" : "data/DIS.dat", "coherent" : 'false', "model" : "Dark_Photon_DM", "gagg" : 0, "gagpg" : 0, "gagpgp" : 0}

def write_experiment(write_detector,user):
    context = defaults.copy()
    context.update(user)
    print("Writing this!")
    print(context)
    prod_chan = context["prod_chan"]; signal_chan = context["signal_chan"]; proddist = context["proddist"]; partlistfile = context["partlistfile"]
    eps = context["eps"]; mdm = context["mdm"]; mv = context["mv"]; alpha_D = context["alpha_D"]; outfile = context["outfile"];
    sumlog = context["sumlog"]; outlog = context["outlog"]; output_mode = context["output_mode"]; samplesize = context["samplesize"];
    min_scatter_energy = context["min_scatter_energy"]; max_scatter_energy=context["max_scatter_energy"]; dm_energy_resolution = context["dm_energy_resolution"];
    efficiency = context["efficiency"]; beam_energy = context["beam_energy"]; n_num_target = context["n_num_target"]; p_num_target = context["p_num_target"];
    max_trials = context["max_trials"]; ptmax = context["ptmax"]; zmin = context["zmin"]; zmax = context["zmax"]; run = context["run"];
    pi0_per_POT=context["pi0_per_POT"]; p_cross = context["p_cross"]; meson_per_pi0 = context["meson_per_pi0"]; min_scatter_angle=context["min_scatter_angle"];
    max_scatter_angle=context["max_scatter_angle"]; repeat = context["repeat"]; timing = context["timing"]; burn_max = context["burn_max"];
    inelastic_dist = context["inelastic_dist"]; coherent = context["coherent"]; model = context["model"]; gagg= context["gagg"]; gagpg = context["gagpg"];
    gagpgp = context["gagpgp"]; POT = context["POT"]
    with open(context["outfile"],'w') as f:
        if run>=0:
            f.write('run {}\n'.format(context["run"]))
        for i in range(len(prod_chan)):
	    f.write('production_channel {}\n'.format(prod_chan[i]))
	    if(prod_chan[i]=="parton_production"):
            	prepare_parton(mA=mv,energy=beam_energy,file_path=v_parton_kinetic)
                f.write('parton_V_neutron_file {}\n'.format("data/parton_V_n.dat"))
                f.write('parton_V_proton_file {}\n'.format("data/parton_V_p.dat"))
            elif(prod_chan[i]=="parton_production_baryonic"):
            	prepare_parton(mA=mv,energy=beam_energy,file_path=v_parton_baryonic)
                f.write('parton_V_neutron_file {}\n'.format("data/parton_V_n.dat"))
                f.write('parton_V_proton_file {}\n'.format("data/parton_V_p.dat"))
            if(proddist[i]!=""):
            	f.write("production_distribution {}\n".format(proddist[i]))
	    if(partlistfile[i]!=""):
		f.write("particle_list_file {}\n".format(partlistfile[i]))
	    if prod_chan[i] in meson_per_pi0:
                f.write('meson_per_pi0 {}\n'.format(str(meson_per_pi0[prod_chan[i]])))
	    if proddist[i]=='proton_brem' or proddist[i]=='proton_brem_baryonic':
                f.write('zmax {}\n'.format(str(zmax)))
                f.write('zmin {}\n'.format(str(zmin)))
                f.write('ptmax {}\n'.format(str(ptmax)))
                f.write('\n')
        if repeat!=1:
            f.write("repeat {}\n".format(str(repeat)))
        if timing!=0.0:
            f.write("timing_cut {}\n".format(str(timing)))
        if burn_max!=-1:
            f.write("burn_max {}\n".format(str(burn_max)))
        f.write('proton_target_cross_section {}\n'.format(str(p_cross)))
        f.write('max_trials {}\n'.format(str(max_trials)))
        f.write("efficiency {}\n".format(str(efficiency)))
        f.write('min_scatter_energy {}\n'.format(str(min_scatter_energy)))
        f.write('max_scatter_energy {}\n'.format(str(max_scatter_energy)))
        f.write('min_scatter_angle {}\n'.format(str(min_scatter_angle)))
        f.write('max_scatter_angle {}\n'.format(str(max_scatter_angle)))
        f.write('dm_energy_resolution {}\n'.format(str(dm_energy_resolution)))
        f.write('epsilon {}\n'.format(str(eps)))
        f.write('n_num_target {}\n'.format(str(n_num_target)))
        f.write('p_num_target {}\n'.format(str(p_num_target)))
        f.write('beam_energy {}\n'.format(str(beam_energy)))
        f.write('dark_matter_mass {}\n'.format(str(mdm)))
        f.write('dark_photon_mass {}\n'.format(str(mv)))
        f.write('alpha_D {}\n'.format(str(alpha_D)))
        if(model=="Axion_Dark_Photon"):
            f.write('gagg {}\n'.format(str(gagg)))
            f.write('gagpg {}\n'.format(str(gagpg)))
            f.write('gagpgp {}\n'.format(str(gagpgp)))
        f.write('model {}\n'.format(model));
        f.write('POT {}\n'.format(str(POT)))
        f.write('signal_channel {}\n'.format(signal_chan))
        if signal_chan=="Inelastic_Nucleon_Scattering" or signal_chan == "Inelastic_Nucleon_Scattering_Baryonic":
            Parsing_inelastic.gen_cross_section(mv,eps=eps,alpha_d=alpha_D,model=signal_chan,output=inelastic_dist)
            f.write('scatter_dist_filename {}\n'.format(inelastic_dist))
        f.write('output_file {}\n'.format(outlog))
        f.write('summary_file {}\n'.format(sumlog))
        f.write('output_mode {}\n'.format(output_mode))
        f.write('samplesize {}\n'.format(str(samplesize)))
        f.write('coherent {}\n'.format(coherent))
        f.write('pi0_per_POT {}\n\n'.format(str(pi0_per_POT)))
        write_detector(f)
        f.close()

##################
#DETECTOR PRESETS#
##################
def miniboone_detector(f,xpos=0.0,ypos=-1.9,zpos=491.0,radius=5.0):
    f.write("\ndetector sphere\n")
    f.write("x-position {0}\ny-position {1}\nz-position {2}\nradius {3}\n".format(str(xpos),str(ypos),str(zpos),str(radius)))
    f.write('\n')
    f.write(Hydrogen_string)
    f.write('\n')
    f.write(Carbon_string)

def miniboone_detector_numi(f,xpos=0.0,ypos=0.0,zpos=100.0,radius=5.0):
    f.write("\ndetector sphere\n")
    f.write("x-position {0}\ny-position {1}\nz-position {2}\nradius {3}\n".format(str(xpos),str(ypos),str(zpos),str(radius)))
    f.write('\n')
    f.write(Hydrogen_string)
    f.write('\n')
    f.write(Carbon_string)

def miniboone_detector_full(f,xpos=0.0,ypos=-1.9,zpos=491.0,radius=6.106):
    f.write("\ndetector sphere\n")
    f.write("x-position {0}\ny-position {1}\nz-position {2}\nradius {3}\n".format(str(xpos),str(ypos),str(zpos),str(radius)))
    f.write('\n')
    f.write(Hydrogen_string)
    f.write('\n')
    f.write(Carbon_string)

def test_sphere(f,xpos=0.0,ypos=0.0,zpos=0.0,radius=1.0):
    f.write("\ndetector sphere\n")
    f.write("x-position {0}\ny-position {1}\nz-position {2}\nradius {3}\n".format(str(xpos),str(ypos),str(zpos),str(radius)))
    f.write('\n')
    f.write(Hydrogen_string)
    f.write('\n')
    f.write(Carbon_string)

def test_cylinder(f,xpos=0.0,ypos=0.0,zpos=0.0,radius=1.0, length=1.0, theta=0, phi=0):
    f.write("\ndetector cylinder\n");
    f.write("x-position {0}\ny-position {1}\nz-position {2}\nradius {3}\nlength {4}\ndet-theta {5}\ndet-phi {6}\n".format(str(xpos),str(ypos),str(zpos),str(radius),str(length),str(theta),str(phi)))
    f.write('\n')
    f.write(Hydrogen_string)
    f.write('\n')
    f.write(Carbon_string)

def test_cuboid(f,xpos=0.0,ypos=0.0,zpos=0.0,height=1.0,length=1.0,width=1.0,theta=0.0,phi=0,psi=0):
    f.write("\ndetector cuboid\n");
    f.write("x-position {0}\ny-position {1}\nz-position {2}\nwidth {3}\nlength {4}\nheight {5}\ndet-phi {6}\ndet-theta {7}\ndet-psi {8}".format(str(xpos),str(ypos),str(zpos),str(width),str(length),str(height),str(phi),str(theta),str(psi)))
    f.write('\n')
    f.write(Hydrogen_string)
    f.write('\n')
    f.write(Carbon_string)

MINOS_absorber_z=270
MINOS_target_z=950
#Don't use this for actual event generation!
def MINOS_detector(f,xpos=0.0,ypos=0.0,zpos=MINOS_target_z,radius=2.2,length=1.7,theta=0,phi=0):
    print("This detector should not be used for event generation!")
    f.write("\ndetector cylinder\n");
    f.write("x-position {0}\ny-position {1}\nz-position {2}\nradius {3}\nlength {4}\ndet-theta {5}\ndet-phi {6}\n".format(str(xpos),str(ypos),str(zpos),str(radius),str(length),str(theta),str(phi)))
    f.write('\n')
    f.write(MINOS_string)

def MINOS_absorber_detector(f):
    MINOS_detector(f,zpos=MINOS_absorber_z)

#NOvA_absorber_d=240
#NOvA_target_d=920
NOvA_absorber_d=447
NOvA_target_d=990
#NOvA_angle=0.0575959#3.3 degree
NOvA_Target_Angle=0.0122
NOvA_Absorber_Angle=0.024

def NOvA_detector(f,xpos=0.0,ypos=NOvA_target_d*math.sin(NOvA_Target_Angle),zpos=NOvA_target_d*math.cos(NOvA_Target_Angle),height=4.2,length=14.3,width=2.9,theta=-NOvA_Target_Angle,phi=0,psi=0):
    print("This NOvA detector is prelimary!")
    #print("x-position {0}\ny-position {1}\nz-position {2}\nwidth {3}\nlength {4}\nheight {5}\ndet-phi {6}\ndet-theta {7}\ndet-psi {8}".format(str(xpos),str(ypos),str(zpos),str(width),str(length),str(height),str(phi),str(theta),str(psi)))
    f.write("\ndetector cuboid\n");
    f.write("x-position {0}\ny-position {1}\nz-position {2}\nwidth {3}\nlength {4}\nheight {5}\ndet-phi {6}\ndet-theta {7}\ndet-psi {8}".format(str(xpos),str(ypos),str(zpos),str(width),str(length),str(height),str(phi),str(theta),str(psi)))
    f.write('\n')
    f.write(NOvA_string)

def NOvA_absorber_detector(f):
    NOvA_detector(f,xpos=0.0,ypos=NOvA_absorber_d*math.sin(NOvA_Absorber_Angle),zpos=NOvA_absorber_d*math.cos(NOvA_Absorber_Angle))

def SBND_detector(f,xpos=0.0,ypos=0,zpos=112.0,radius=2.38,length=4.76,theta=0,phi=0):
    f.write("\ndetector cylinder\n");
    f.write("x-position {0}\ny-position {1}\nz-position {2}\nradius {3}\nlength {4}\ndet-theta {5}\ndet-phi {6}\n".format(str(xpos),str(ypos),str(zpos),str(radius),str(length),str(theta),str(phi)))
    f.write('\n')
    f.write(Argon_string)

def SBND_detector_old(f,xpos=0.0,ypos=0,zpos=62.0,width=2.38,length=4.76,theta=0,phi=0):
    f.write("\ndetector cylinder\n");
    f.write("x-position {0}\ny-position {1}\nz-position {2}\nradius {3}\nlength {4}\ndet-theta {5}\ndet-phi {6}\n".format(str(xpos),str(ypos),str(zpos),str(radius),str(length),str(theta),str(phi)))
    f.write('\n')
    f.write(Argon_string)

SHIP_Side=1.92856

def ship_detector(f,xpos=0.0,ypos=0,zpos=100.0,width=SHIP_Side,length=SHIP_Side,height=SHIP_Side,phi=0,theta=0,psi=0):
    #def ship_detector(f,xpos=0.0,ypos=0,zpos=30.0,radius=0.655,length=2.645,theta=0,phi=0):
    f.write("\ndetector cuboid\n");
    f.write("x-position {0}\ny-position {1}\nz-position {2}\nwidth {3}\nlength {4}\nheight {5} \ndet-theta {6}\ndet-phi {7}\ndet-psi {8}".format(str(xpos),str(ypos),str(zpos),str(width),str(length),str(height),str(phi),str(theta),str(psi)))
    f.write('\n')
    #Need to figure out what it's made of!
    f.write(Argon_string)

#temp detector until I implement proper geometry handling
#This is actually the P0D, with only the proper number of neutrons and protons. NO ATOMS IMPLEMENTED
#Double check the fiducial mass on the pod. Is it 3 tons of water? 13 tons of stuff?
#Need to update with cuboid shape!
def t2k_ND280(f):
    xpos=11;ypos=0;zpos=280;detphi=0;radius=0.9413;dettheta=0.0436332;length=1.7;
    f.write("\ndetector cylinder\n")
    f.write("x-position {0}\ny-position {1}\nz-position {2}\nradius {3}\ndet-theta {4}\ndet-phi {5}\nlength {6}\n".format(str(xpos),str(ypos),str(zpos),str(radius),str(dettheta),str(detphi),str(length),str(length)))
    f.write('\n')
    f.write(ND280_string)

def t2k_superK1000(f):
    xpos=12867.7;ypos=0;zpos=294719;detphi=0;radius=190.5;dettheta=1.5708;length=410;
    f.write("\ndetector cylinder\n")
    f.write("x-position {0}\ny-position {1}\nz-position {2}\nradius {3}\ndet-theta {4}\ndet-phi {5}\nlength {6}\n".format(str(xpos),str(ypos),str(zpos),str(radius),str(dettheta),str(detphi),str(length),str(length)))
    print("WARNING: This detector is enlarged 1000 fold. Its efficiency must be suppressed in a corresponding manner.")
    f.write('\n')
    f.write(Water_string)

def lsnd_detector(f,xpos=0.0,ypos=-4.65,zpos=29.8,radius=5.7/2.0-0.35,length=8.3,theta=0,phi=0):
    f.write("\ndetector cylinder\n");
    f.write("x-position {0}\ny-position {1}\nz-position {2}\nradius {3}\nlength {4}\ndet-theta {5}\ndet-phi {6}\n".format(str(xpos),str(ypos),str(zpos),str(radius),str(length),str(theta),str(phi)))
    f.write('\n')
    f.write(Carbon_string)
    f.write('\n')
    f.write(Hydrogen_string)

def coherent_detector_LAr(f,xpos=20.0,ypos=0.0,zpos=0.0,radius=0.48,length=0.96,theta=pi/2.0,phi=pi/2.0):
    f.write("\ndetector cylinder\n");
    f.write("x-position {0}\ny-position {1}\nz-position {2}\nradius {3}\nlength {4}\ndet-theta {5}\ndet-phi {6}\n".format(str(xpos),str(ypos),str(zpos),str(radius),str(length),str(theta),str(phi)))
    f.write('\n')
    f.write(Argon_string)

def coherent_detector_NaI(f,xpos=20.0,ypos=0.0,zpos=0.0,radius=0.601,length=2*0.601,theta=pi/2.0,phi=pi/2.0):
    f.write("\ndetector cylinder\n");
    f.write("x-position {0}\ny-position {1}\nz-position {2}\nradius {3}\nlength {4}\ndet-theta {5}\ndet-phi {6}\n".format(str(xpos),str(ypos),str(zpos),str(radius),str(length),str(theta),str(phi)))
    f.write('\n')
    f.write(Sodium_Iodide_string)
0.08
def coherent_detector_CsI(f,xpos=19.6,ypos=0.0,zpos=0.0,radius=0.08,length=2*0.08,theta=pi/2.0,phi=pi/2.0):
    f.write("\ndetector cylinder\n");
    f.write("x-position {0}\ny-position {1}\nz-position {2}\nradius {3}\nlength {4}\ndet-theta {5}\ndet-phi {6}\n".format(str(xpos),str(ypos),str(zpos),str(radius),str(length),str(theta),str(phi)))
    f.write('\n')
    f.write(Cesium_Iodide_string)

def coherent_detector_CsI_1T(f,xpos=20,ypos=0.0,zpos=0.0,radius=0.328,length=2*0.328,theta=pi/2.0,phi=pi/2.0):
    f.write("\ndetector cylinder\n");
    f.write("x-position {0}\ny-position {1}\nz-position {2}\nradius {3}\nlength {4}\ndet-theta {5}\ndet-phi {6}\n".format(str(xpos),str(ypos),str(zpos),str(radius),str(length),str(theta),str(phi)))
    f.write('\n')
    f.write(Cesium_Iodide_string)

def captain_detector(f,xpos=0.0,ypos=0.0,zpos=30.0,radius=1,length=1.15,theta=pi/2.0,phi=0.0):
    f.write("\ndetector cylinder\n");
    f.write("x-position {0}\ny-position {1}\nz-position {2}\nradius {3}\nlength {4}\ndet-theta {5}\ndet-phi {6}\n".format(str(xpos),str(ypos),str(zpos),str(radius),str(length),str(theta),str(phi)))
    f.write('\n')
    f.write(Argon_string)

def captain_detector_off(f,xpos=30.0,ypos=0.0,zpos=0.0,radius=1,length=1.15,theta=pi/2.0,phi=0.0):
    f.write("\ndetector cylinder\n");
    f.write("x-position {0}\ny-position {1}\nz-position {2}\nradius {3}\nlength {4}\ndet-theta {5}\ndet-phi {6}\n".format(str(xpos),str(ypos),str(zpos),str(radius),str(length),str(theta),str(phi)))
    f.write('\n')
    f.write(Argon_string)

#############
#EXPERIMENTS#
#############

miniboone_default = {}

def write_miniboone(d={},det=miniboone_detector):
    context = miniboone_default.copy()
    context.update(d)
    write_experiment(det,context)

t2k_default = {"POT" : 8e21, "proddist" : ["bmpt"], "partlistfile" : ["data/particle_list_t2k.dat"], "sumlog" : "Events/t2k.dat", "outlog" : "Events/t2k_events.dat", "beam_energy" : 30, "n_num_target" : 6, "p_num_target" : 6, "ptmax" : 1, "zmin" : 0.2, "zmax" : 0.8, "pi0_per_POT" : 1.0, "meson_per_pi0" : meson_per_pi0_miniboone, "min_scatter_energy" : 0.035, "max_scatter_energy" : 2, "p_cross" : 15*mb}

def write_t2k(d={}, det=t2k_ND280):
    context = t2k_default.copy()
    context.update(d)
    write_experiment(det,context)

SHIP_Energy = 400

ship_default = {"proddist" : ["bmpt"], "partlistfile" : ["data/particle_list_ship.dat"], "sumlog" : "Events/ship.dat", "outlog" : "Events/ship_events.dat", "beam_energy" : SHIP_Energy, "n_num_target" : 54, "p_num_target" : 42, "ptmax" : 1, "zmin" : 0.1, "zmax" : 0.9, "signal_chan" : "NCE_electron", "min_scatter_energy" : 2, "max_scatter_energy" : 20, "min_scatter_angle" : 0.01, "max_scatter_angle" : 0.02, "samplesize" : 1000, "efficiency" : 0.5, "POT" : 2e20, "pi0_per_POT" : 1.8, "p_cross" :
        11*mb, "meson_per_pi0" : meson_per_pi0_ship}

def write_ship(d={}, det=ship_detector):
    context = ship_default.copy()
    context.update(d)
    write_experiment(det,context)

lsnd_default = {"proddist" : ["burmansmith"], "partlistfile" : ["data/particle_list_lsnd.dat"], "signal_chan" : "NCE_electron", "sumlog" : "Events/lsnd.dat", "outlog" : "Events/lsnd_events.dat", "samplesize" : 5000, "min_scatter_energy" : 0.018, "max_scatter_energy" : 0.05, "efficiency" : 0.19*0.793, "beam_energy" : 0.8, "n_num_target" : 0, "p_num_target" : 1, "POT" : 1.8e23, "pi0_per_POT" : 0.06, "p_cross" : 30*mb}

def write_lsnd(d={}, det=lsnd_detector):
    context = lsnd_default.copy()
    context.update(d)
    write_experiment(det,context)

coherent_default = {"coherent" : "true", "proddist" : ["burmansmith"], "partlistfile" : ["data/particle_list_coherent.dat"], "sumlog" : "Events/coherent.dat", "outlog" : "Events/coherent_events.dat", "min_scatter_energy" : 0.018, "max_scatter_energy" : 0.05, "dm_energy_resolution" : 0.001, "efficiency" : 0.5, "beam_energy" : 1.0, "p_num_target" : 80, "n_num_target" : 0, "POT" : 1e23, "pi0_per_POT" : 0.1, "p_cross" : 30*mb, "burn_max" : 100, "meson_per_pi0" :
        meson_per_pi0_coherent, "outlog" : "Events/coherent_events.dat"}

def write_coherent(d={}, det=coherent_detector_LAr):
    context = coherent_default.copy()
    context.update(d)
    write_experiment(det,context)

numi_energy = 120

numi_default = {"proddist" : ["bmpt"], "partlistfile" : ["data/particle_list_numi.dat"], "sumlog" : "Events/nova.dat", "outlog" : "Events/nova_events.dat", "output_mode" : "summary", "samplesize" : 5000, "min_scatter_energy" : 0.05, "max_scatter_energy" : 3.0, "dm_energy_resolution" : 0.01, "efficiency" : 0.35, "beam_energy" : 120, "n_num_target" : 8, "p_num_target" : 8, "ptmax" : 2, "zmin" : 0.1, "zmax" : 0.9, "POT" : 1e21,
        "pi0_per_POT" : 1.0, "p_cross" : 15*mb}

def write_numi(d={}, det=NOvA_detector):
    context = numi_default.copy()
    context.update(d)
    write_experiment(det,context)

v_parton_kinetic = "~/Code/DMcode/"
v_parton_baryonic = "~/Code/DMcodeBaryon/"

miniboone_Energy= 8.9

local_directory=os.getcwd()

def prepare_parton(mA=1,energy=miniboone_Energy,file_path=v_parton_kinetic):
    with cd(file_path):
        subp.call(["./scalarDM",str(mA),"1",str(energy)])
        subp.call(["mv","dsigdpV.out","{}/data/parton_V_p.dat".format(local_directory)])
        subp.call(["./scalarDM",str(mA),"2",str(energy)])
        subp.call(["mv","dsigdpV.out","{}/data/parton_V_n.dat".format(local_directory)])
    textarr=["data/parton_V_p.dat","data/parton_V_n.dat"]
    for sfile in textarr:
        df = pd.read_csv(sfile, sep=" ", engine='python',skipinitialspace=True,header=None)
        arr= df.values.copy()
        arr2 = [[i, np.interp(i,arr[:,0],arr[:,1])] for i in np.arange(arr[0,0],arr[-1,0],0.1)]
        np.savetxt(sfile,arr2,delimiter=" ")

