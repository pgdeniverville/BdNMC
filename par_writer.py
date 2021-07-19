import subprocess as subp
from contextlib import contextmanager
import os
import pandas as pd
import numpy as np
from multiprocessing import Pool

from par_writer_constants import *
from par_writer_detectors import *

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

defaults = {"eps" : 1e-3, "mdm" : 0.03, "mv" : 0.1, "alpha_D" : 0.1, "prod_chan" : ["pi0_decay"], "signal_chan" : "NCE_nucleon", "outfile" : "parameter_run.dat", "proddist" : [""], "partlistfile" : ["Source/particle_list.dat"], "sumlog" : "Events/miniboone.dat", "outlog" : "Events/miniboone_events.dat", "output_mode" :"summary", "samplesize" : 5000, "min_scatter_energy" : 0.035, "max_scatter_energy" : 1.0, "dm_energy_resolution" : 0.01, "efficiency" : 0.35, "beam_energy" : 8.9, "n_num_target" : 4, "p_num_target" : 4, "max_trials" : 80e6, "ptmax" : 0.2, "zmin" : 0.3, "zmax" : 0.7, "run" : -1, "POT" : 2e20, "pi0_per_POT" : 0.9, "p_cross" : 25*mb, "meson_per_pi0" : meson_per_pi0_miniboone, "min_scatter_angle" : 0.0, "max_scatter_angle" : 2.1*pi, "repeat" : 1, "timing" : 0.0, "burn_max" : -1,"inelastic_dist" : "data/DIS.dat", "coherent" : 'false', "model" : "Dark_Photon_DM", "gagg" : 0, "gagpg" : 0, "gagpgp" : 0,'weighted' : "false", "particle_list_position" : ["false"], "epsilon_W" : 0, "epsilon_l" : 1e-3, "epsilon_q" : 1e-3, "dark_scalar_mass" : 0.01}

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
    gagpgp = context["gagpgp"]; POT = context["POT"];
    particle_list_position = context["particle_list_position"]
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
            if len(particle_list_position) >= i:
                f.write("particle_list_position {}\n".format(particle_list_position))
            if prod_chan[i] in meson_per_pi0:
                f.write('meson_per_pi0 {}\n'.format(str(meson_per_pi0[prod_chan[i]])))
            if proddist[i]=='proton_brem' or proddist[i]=='proton_brem_baryonic' or prod_chan[i]=="proton_brem":
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
        f.write('weighted {}\n'.format(context['weighted']))
        f.write('epsilon {}\n'.format(str(eps)))
        f.write('n_num_target {}\n'.format(str(n_num_target)))
        f.write('p_num_target {}\n'.format(str(p_num_target)))
        f.write('beam_energy {}\n'.format(str(beam_energy)))
        write_model(f,context)
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
        if "kinetic_energy_cut" in context:
            f.write('kinetic_energy_cut {}'.format(context["kinetic_energy_cut"]))
        write_detector(f)
        f.close()

def write_model(f,context):
    model = str(context["model"])
    f.write('model {}\n'.format(model));
    if(model=="Dark_Scalar"):
        f.write("epsilon_l {}\n".format(str(context["epsilon_l"])))
        f.write("epsilon_q {}\n".format(str(context["epsilon_q"])))
        f.write("epsilon_W {}\n".format(str(context["epsilon_W"])))
        f.write("dark_scalar_mass {}\n".format(str(context["dark_scalar_mass"])))
    if(model=="Inelastic_Dark_Matter"):
        f.write('dark_matter_1_mass {}\n'.format(str(context["mdm1"])))
        f.write('dark_matter_2_mass {}\n'.format(str(context["mdm2"])))
    if(model=="Dark_Photon_DM" or model == "Dark_Photon" or model=="Inelastic_Dark_Matter" or model == "Axion_Dark_Photon"):
        f.write('dark_photon_mass {}\n'.format(str(context["mv"])))
        f.write('epsilon {}\n'.format(str(context["eps"])))
    if(model=="Dark_Photon_DM" or model=="Inelastic_Dark_Matter" or model == "Dark_Photon"):
        f.write('alpha_D {}\n'.format(context["alpha_D"]))
    if(model == "Dark_Photon_DM" or model == "Axion_Dark_Photon" or model == "Dark_Photon"):
        f.write('dark_matter_mass {}\n'.format(str(context["mdm"])))
    if(model=="Axion_Dark_Photon"):
        f.write('gagg {}\n'.format(str(context["gagg"])))
        f.write('gagpg {}\n'.format(str(context["gagpg"])))
        f.write('gagpgp {}\n'.format(str(context["gagpgp"])))
    if(model=="Axion_Dark_Photon"):
        f.write("axion_mass {}\n".format(str(context["axion_mass"])))

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

LANL_test_default = {"coherent" : "true", "proddist" : ["burmansmith"], "partlistfile" : ["data/particle_list_coherent.dat"], "sumlog" : "Events/coherent.dat", "outlog" : "Events/coherent_events.dat", "min_scatter_energy" : 0.018, "max_scatter_energy" : 0.05, "dm_energy_resolution" : 0.001, "efficiency" : 0.5, "beam_energy" : 0.8, "p_num_target" : 74, "n_num_target" : 109, "POT" : 1.15e21, "pi0_per_POT" : 0.0425, "p_cross" : 30*mb, "burn_max" : 100, "meson_per_pi0" :
        meson_per_pi0_coherent, "outlog" : "Events/coherent_events.dat"}

def write_lanl(d={}, det=lanl_detector):
    context = LANL_test_default.copy()
    context.update(d)
    write_experiment(det,context)

SHIP_Energy = 400

ship_default = {"proddist" : ["bmpt"], "partlistfile" : ["data/particle_list_ship.dat"], "sumlog" : "Events/ship.dat", "outlog" : "Events/ship_events.dat", "beam_energy" : SHIP_Energy, "n_num_target" : 54, "p_num_target" : 42, "ptmax" : 1, "zmin" : 0.1, "zmax" : 0.9, "signal_chan" : "NCE_electron", "min_scatter_energy" : 1, "max_scatter_energy" : 20, "min_scatter_angle" : 0.01, "max_scatter_angle" : 0.02, "samplesize" : 1000, "efficiency" : 0.5, "POT" : 2e20, "pi0_per_POT" : 1.8, "p_cross" : 11*mb, "meson_per_pi0" : meson_per_pi0_ship}

def write_ship(d={}, det=ship_detector):
    context = ship_default.copy()
    context.update(d)
    write_experiment(det,context)

lsnd_default = {"proddist" : ["burmansmith"], "partlistfile" : ["data/particle_list_lsnd.dat"], "signal_chan" : "NCE_electron", "sumlog" : "Events/lsnd.dat", "outlog" : "Events/lsnd_events.dat", "samplesize" : 5000, "min_scatter_energy" : 0.018, "max_scatter_energy" : 0.05, "efficiency" : 0.19*0.793, "beam_energy" : 0.8, "n_num_target" : 0, "p_num_target" : 1, "POT" : 1.8e23, "pi0_per_POT" : 0.06, "p_cross" : 30*mb}

def write_lsnd(d={}, det=lsnd_detector):
    context = lsnd_default.copy()
    context.update(d)
    write_experiment(det,context)

coherent_default = {"coherent" : "true", "proddist" : ["burmansmith"], "partlistfile" : ["data/particle_list_coherent.dat"], "sumlog" : "Events/coherent.dat", "outlog" : "Events/coherent_events.dat", "min_scatter_energy" : 15e-6, "max_scatter_energy" : 0.05, "dm_energy_resolution" : 0.001, "efficiency" : 0.5, "beam_energy" : 1.0, "p_num_target" : 80, "n_num_target" : 0, "POT" : 1e23, "pi0_per_POT" : 0.06, "p_cross" : 30*mb, "burn_max" : 100, "meson_per_pi0" :
        meson_per_pi0_coherent, "outlog" : "Events/coherent_events.dat"}

def write_coherent(d={}, det=coherent_detector_LAr):
    context = coherent_default.copy()
    context.update(d)
    write_experiment(det,context)

numi_energy = 120

numi_default = {"proddist" : ["bmpt"], "partlistfile" : ["data/particle_list_numi.dat"], "sumlog" : "Events/nova.dat", "outlog" : "Events/nova_events.dat", "output_mode" : "summary", "samplesize" : 5000, "min_scatter_energy" : 0.05, "max_scatter_energy" : 3.0, "dm_energy_resolution" : 0.01, "efficiency" : 0.35, "beam_energy" : 120, "n_num_target" : 8, "p_num_target" : 8, "ptmax" : 2, "zmin" : 0.1, "zmax" : 0.9, "POT" : 1e21,
        "pi0_per_POT" : 2.9, "p_cross" : 15*mb, "meson_per_pi0" : meson_per_pi0_numi}

def write_numi(d={}, det=NOvA_detector):
    context = numi_default.copy()
    context.update(d)
    write_experiment(det,context)

CHARM2_Energy=450

#No idea what the target is made of.
charm2_default = {"proddist" : ["bmpt"], "partlistfile" : ["data/particle_list_charm2.dat"], "sumlog" : "Events/charm2.dat", "outlog" : "Events/charm2_events.dat", "beam_energy" : CHARM2_Energy, "n_num_target" : 34, "p_num_target" : 29, "ptmax" : 1, "zmin" : 0.1, "zmax" : 0.9, "signal_chan" : "NCE_electron", "min_scatter_energy" : 3, "max_scatter_energy" : 24, "min_scatter_angle" : 0, "max_scatter_angle" : pi, "samplesize" : 1000, "efficiency" : 1, "POT" : 2.5e19, "pi0_per_POT" : 2.4, "p_cross" : 11*mb, "meson_per_pi0" : meson_per_pi0_ship}

def write_charm2(d={}, det=CHARM2_detector):
    context = charm2_default.copy()
    context.update(d)
    write_experiment(det,context)

CHARM_Energy=400

charm_default = {"proddist" : ["bmpt"], "partlistfile" : ["data/particle_list_charm.dat"], "sumlog" : "Events/charm.dat", "outlog" : "Events/charm_events.dat", "beam_energy" : CHARM_Energy, "n_num_target" : 34, "p_num_target" : 29, "ptmax" : 1, "zmin" : 0.1, "zmax" : 0.9, "signal_chan" : "Signal_Decay", "min_scatter_energy" : 1, "max_scatter_energy" : 400, "min_scatter_angle" : 0, "max_scatter_angle" : pi, "samplesize" : 1000, "efficiency" : 1, "POT" : 2.4e18, "pi0_per_POT" : 2.4, "p_cross" : 30*mb, "meson_per_pi0" : meson_per_pi0_ship}

def write_charm_decay(d={}, det=CHARM_decay_detector):
    context = charm_default.copy()
    context.update(d)
    write_experiment(det,context)


pip2_default = {"proddist" : ["burmansmith"], "partlistfile" : ["data/particle_list_pip3.dat"], "signal_chan" : "NCE_nucleon", "sumlog" : "Events/pip2.dat", "outlog" : "Events/pip2_events.dat", "samplesize" : 5000, "min_scatter_energy" : 20e-6, "max_scatter_energy" : 0.05, "efficiency" : 0.8, "beam_energy" : 0.8, "n_num_target" : 0, "p_num_target" : 1, "POT" : 2e22, "pi0_per_POT" : 0.12, "p_cross" : 30*mb, "kinetic_energy_cut" : "true"}

def write_pip2(d={}, det = pip2_coherent_det):
    context = pip2_default.copy()
    context.update(d)
    write_experiment(det,context)

nucal_energy=70
#How many pi0/POT? 38 mb pp cross from pdg, 78 mb pp->pi^0 X from https://arxiv.org/pdf/1104.2747.pdf
nucal_default = {"proddist" : ["bmpt"], "meson_per_pi0" : meson_per_pi0_numi, "partlistfile" : ["data/particle_list_nucal.dat"], "sumlog" : "Events/nucal.dat", "outlog" : "Events/nucal_events.dat", "output_mode" : "summary", "samplesize" : 100000, "min_scatter_energy" : 3, "max_scatter_energy" : 70, "dm_energy_resolution" : 0.01, "efficiency" : 1, "beam_energy" : 70, "n_num_target" : 30, "p_num_target" : 26, "ptmax" : 2, "zmin" : 0.1, "zmax" : 0.9, "POT" : 1.7e18, "pi0_per_POT" : 74/38, "p_cross" : 38*mb, 'max_scatter_angle' : 0.05}

def write_nucal(d={}, det=nucal):
    context = nucal_default.copy()
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

