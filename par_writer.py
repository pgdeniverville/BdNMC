import subprocess as subp
from contextlib import contextmanager
import os
import pandas as pd
import numpy as np
from multiprocessing import Pool

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
pi=3.14159

meson_per_pi0_miniboone = {'pi0_decay' : '1.0', 'eta_decay' : str(1.0/30.0), 'rho_decay' : str(1.0/20.0), 'omega_decay' : '0.046', 'phi_decay' : str(1.0/150.0), 'pi0_decay_baryonic' : '1.0', 'eta_decay_baryonic' : str(1.0/30.0), 'rho_decay_baryonic' : str(1.0/20.0), 'omega_decay_baryonic' : '0.046', 'phi_decay_baryonic' : str(1.0/150.0)}

meson_per_pi0_lsnd = {'pi0_decay' : '1.0'}

meson_per_pi0_ship = {'pi0_decay' : '1.0', 'eta_decay' : str(0.078), 'rho_decay' : str(0.11), 'omega_decay' : '0.11', 'phi_decay' : str(0.02)}

Hydrogen_string = "material Hydrogen\nnumber_density 7.26942e22\nproton_number 1\nneutron_number 0\nelectron_number 1\nmass 0.945778\n"

Carbon_string = "material Carbon\nnumber_density 3.63471e22\nproton_number 6\nneutron_number 6\nelectron_number 6\nmass 11.2593\n"

Argon_string = "material Argon\nnumber_density 2.11e22\nproton_number 18\nneutron_number 22\nelectron_number 18\nmass {0}\n".format(str(39.948*0.938))

ND280_string = "material nd280stuff\nnumber_density 3.7e23\nproton_number 1\nneutron_number 1\nelectron_number 1\nmass 0.945778\n"

def write_experiment(write_detector,eps=1e-3, mdm = 0.03, mv = 0.1, alpha_D = 0.1, prod_chan = ["pi0_decay"], signal_chan = "NCE_nucleon", outfile="parameter_run.dat", proddist=[""], partlistfile=["Source/particle_list.dat"],sumlog="Events/miniboone.dat",outlog="Events/miniboone_events.dat", output_mode="summary",samplesize=5000, min_scatter_energy=0.035, max_scatter_energy=1.0, dm_energy_resolution=0.01, efficiency=0.35,beam_energy=8.9, n_num_target=4,p_num_target=4,max_trials=80e6,ptmax=0.2,zmin=0.3,zmax=0.7,run=-1,POT=2e20,pi0_per_POT=0.9,p_cross=25*mb,meson_per_pi0=meson_per_pi0_miniboone,min_scatter_angle=0.0,max_scatter_angle=2.1*pi):
    with open(outfile,'w') as f:
        if run>=0:
            f.write('run {}\n'.format(run))
	for i in xrange(len(prod_chan)):
	    f.write('production_channel {}\n'.format(prod_chan[i]))
	    if(prod_chan[i]=="parton_production"):
            	prepare_parton(mA=mv,energy=beam_energy,file_path=v_parton_kinetic)
                f.write('parton_V_neutron_file {}\n'.format("Source/parton_V_n.dat"))
                f.write('parton_V_proton_file {}\n'.format("Source/parton_V_p.dat"))
            elif(prod_chan[i]=="parton_production_baryonic"):
            	prepare_parton(mA=mv,energy=beam_energy,file_path=v_parton_baryonic)
                f.write('parton_V_neutron_file {}\n'.format("Source/parton_V_n.dat"))
                f.write('parton_V_proton_file {}\n'.format("Source/parton_V_p.dat"))
            if(proddist[i]!=""):
            	f.write("production_distribution {}\n".format(proddist[i]))
	    if(partlistfile[i]!=""):
		f.write("particle_list_file {}\n".format(partlistfile[i]))
	    if(prod_chan[i]!='parton_production_baryonic' and prod_chan[i]!='parton_production' and prod_chan[i]!='V_decay'):
            	f.write('meson_per_pi0 {}\n'.format(str(meson_per_pi0[prod_chan[i]])))
	    if(proddist[i]=='proton_brem'):
		f.write('ptmax {}\n'.format(str(ptmax)))	
		f.write('zmax {}\n'.format(str(zmax)))	
	        f.write('zmin {}\n'.format(str(zmin)))
            f.write('\n')
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
        f.write('POT {}\n'.format(str(POT)))
        f.write('signal_channel {}\n'.format(signal_chan))
        f.write('output_file {}\n'.format(outlog))
        f.write('summary_file {}\n'.format(sumlog))
        f.write('output_mode {}\n'.format(output_mode))
        f.write('samplesize {}\n'.format(str(samplesize)))
        f.write('pi0_per_POT {}\n\n'.format(str(pi0_per_POT)))
        write_detector(f)
        f.close()

def miniboone_detector(f,xpos=0.0,ypos=-1.9,zpos=491.0,radius=5.0):
    f.write("\ndetector sphere\n")
    f.write("x-position {0}\ny-position {1}\nz-position {2}\nradius {3}\n".format(str(xpos),str(ypos),str(zpos),str(radius)))
    f.write('\n')
    f.write(Hydrogen_string)
    f.write('\n')
    f.write(Carbon_string)

#temp detector until I implement proper geometry handling
#This is actually the P0D, with only the proper number of neutrons and protons. NO ATOMS IMPLEMENTED
def t2k_ND280(f):
    xpos=11;ypos=0;zpos=240;detphi=0;radius=0.9413;dettheta=0.0436332;length=1.7;
    f.write("\ndetector cylinder\n")
    f.write("x-position {0}\ny-position {1}\nz-position {2}\nradius {3}\ndet-theta {4}\ndet-phi {5}\nlength {6}\n".format(str(xpos),str(ypos),str(zpos),str(radius),str(dettheta),str(detphi),str(length),str(length)))
    f.write('\n')
    f.write(ND280_string)

def write_miniboone(eps=1e-3, mdm = 0.03, mv = 0.1, alpha_D = 0.1, prod_chan = ["pi0_decay"], signal_chan = "NCE_nucleon", outfile="parameter_run.dat", proddist=[""], partlistfile=["Source/particle_list.dat"],sumlog="Events/miniboone.dat",outlog="Events/miniboone_events.dat", output_mode="summary",samplesize=5000, min_scatter_energy=0.035, max_scatter_energy=1.0, dm_energy_resolution=0.01, efficiency=0.35,beam_energy=8.9, n_num_target=4,p_num_target=4,max_trials=80e6,ptmax=0.2,zmin=0.3,zmax=0.7,run=-1,min_scatter_angle=0.0,max_scatter_angle=2.1*pi):
    POT=2e20
    pi0_per_POT=0.9
    p_cross=25*mb
    write_experiment(write_detector=miniboone_detector,eps=eps,mdm=mdm,mv=mv,alpha_D=alpha_D,prod_chan=prod_chan,signal_chan = signal_chan, outfile=outfile, proddist=proddist, partlistfile=partlistfile,sumlog=sumlog,outlog=outlog, output_mode=output_mode,samplesize=samplesize, min_scatter_energy=min_scatter_energy, max_scatter_energy=max_scatter_energy, dm_energy_resolution=dm_energy_resolution, efficiency=efficiency,beam_energy=beam_energy, n_num_target=n_num_target,p_num_target=p_num_target,max_trials=max_trials,ptmax=ptmax,zmin=zmin,zmax=zmax,run=run,POT=POT,pi0_per_POT=pi0_per_POT,p_cross=p_cross)

def write_nd280(eps=1e-3, mdm = 0.03, mv = 0.1, alpha_D = 0.1, prod_chan = ["pi0_decay"], signal_chan = "NCE_nucleon", outfile="parameter_run.dat", proddist=["bmpt"], partlistfile=["Source/particle_list.dat"],sumlog="Events/t2k.dat",outlog="Events/t2k_events.dat", output_mode="summary",samplesize=500, min_scatter_energy=0.035, max_scatter_energy=2.0, dm_energy_resolution=0.01, efficiency=0.35,beam_energy=30, n_num_target=6,p_num_target=6,max_trials=80e6,ptmax=1,zmin=0.2,zmax=0.8,run=-1):
    POT=5e21
    pi0_per_POT=1.0
    #Don't take this too seriously
    p_cross=15*mb
    write_experiment(write_detector=t2k_ND280,eps=eps,mdm=mdm,mv=mv,alpha_D=alpha_D,prod_chan=prod_chan,signal_chan = signal_chan, outfile=outfile, proddist=proddist, partlistfile=partlistfile,sumlog=sumlog,outlog=outlog, output_mode=output_mode,samplesize=samplesize, min_scatter_energy=min_scatter_energy, max_scatter_energy=max_scatter_energy, dm_energy_resolution=dm_energy_resolution, efficiency=efficiency,beam_energy=beam_energy, n_num_target=n_num_target,p_num_target=p_num_target,max_trials=max_trials,ptmax=ptmax,zmin=zmin,zmax=zmax,run=run,POT=POT,pi0_per_POT=pi0_per_POT,p_cross=p_cross)

def ship_detector(f,xpos=0.0,ypos=0,zpos=100.0,radius=0.52,length=2.1,theta=0,phi=0):
    f.write("\ndetector cylinder\n");
    f.write("x-position {0}\ny-position {1}\nz-position {2}\nradius {3}\nlength {4}\ndet-theta {5}\ndet-phi {6}\n".format(str(xpos),str(ypos),str(zpos),str(radius),str(length),str(theta),str(phi)))
    f.write('\n')
    f.write(Argon_string)

def write_ship(eps=1e-3, mdm = 0.03, mv = 0.1, alpha_D = 0.1, prod_chan = ["pi0_decay"], signal_chan = "NCE_electron", outfile="parameter_run.dat", proddist=["bmpt"], partlistfile=["Source/particle_list_ship.dat"],outlog="Events/ship_events.dat", output_mode="summary",samplesize=1000,beam_energy=400,n_num_target=54,p_num_target=42,min_scatter_energy=2,max_scatter_energy=20,min_scatter_angle=0.01,max_scatter_angle=0.02,efficiency=0.5,sumlog="Events/ship_run.dat",max_trials=80e6,ptmax=1,zmin=0.1,zmax=0.9,run=-1,dm_energy_resolution=0.01):
    POT=2e20
    pi0_per_POT=1.8
    p_cross=11*mb
    write_experiment(write_detector=ship_detector,eps=eps,mdm=mdm,mv=mv,alpha_D=alpha_D,prod_chan=prod_chan,signal_chan = signal_chan, outfile=outfile, proddist=proddist, partlistfile=partlistfile,sumlog=sumlog,outlog=outlog, output_mode=output_mode,samplesize=samplesize, min_scatter_energy=min_scatter_energy, max_scatter_energy=max_scatter_energy, dm_energy_resolution=dm_energy_resolution, efficiency=efficiency,beam_energy=beam_energy, n_num_target=n_num_target,p_num_target=p_num_target,max_trials=max_trials,ptmax=ptmax,zmin=zmin,zmax=zmax,run=run,POT=POT,pi0_per_POT=pi0_per_POT,p_cross=p_cross,meson_per_pi0=meson_per_pi0_ship,min_scatter_angle=min_scatter_angle,max_scatter_angle=max_scatter_angle)

def lsnd_detector(f,xpos=0.0,ypos=-4.65,zpos=29.8,radius=5.7/2.0-0.35,length=8.3,theta=0,phi=0):
    f.write("\ndetector cylinder\n");
    f.write("x-position {0}\ny-position {1}\nz-position {2}\nradius {3}\nlength {4}\ndet-theta {5}\ndet-phi {6}\n".format(str(xpos),str(ypos),str(zpos),str(radius),str(length),str(theta),str(phi)))
    f.write('\n')
    f.write(Carbon_string)
    f.write('\n')
    f.write(Hydrogen_string)

def write_lsnd(eps=1e-3, mdm = 0.03, mv = 0.1, alpha_D = 0.1, prod_chan = ["pi0_decay"], signal_chan = "NCE_electron", outfile="parameter_run.dat", proddist=["burmansmith"], partlistfile=["Source/particle_list_lsnd.dat"],sumlog="Events/lsnd.dat",outlog="Events/lsnd_events.dat", output_mode="summary",samplesize=5000, min_scatter_energy=0.018, max_scatter_energy=0.05, dm_energy_resolution=0.01, efficiency=0.19*0.793,beam_energy=0.8, n_num_target=0,p_num_target=1,max_trials=80e6,ptmax=0.2,zmin=0.3,zmax=0.7,run=-1,min_scatter_angle=0.0,max_scatter_angle=2.1*pi,POT=1.8e23): 
    pi0_per_POT=0.06
    #Just a random value. I should code the actual function in at some point.
    p_cross=30*mb
    write_experiment(write_detector=lsnd_detector,eps=eps,mdm=mdm,mv=mv,alpha_D=alpha_D,prod_chan=prod_chan,signal_chan = signal_chan, outfile=outfile, proddist=proddist, partlistfile=partlistfile,sumlog=sumlog,outlog=outlog, output_mode=output_mode,samplesize=samplesize, min_scatter_energy=min_scatter_energy, max_scatter_energy=max_scatter_energy, dm_energy_resolution=dm_energy_resolution, efficiency=efficiency,beam_energy=beam_energy, n_num_target=n_num_target,p_num_target=p_num_target,max_trials=max_trials,ptmax=ptmax,zmin=zmin,zmax=zmax,run=run,POT=POT,pi0_per_POT=pi0_per_POT,p_cross=p_cross, min_scatter_angle=min_scatter_angle,max_scatter_angle=max_scatter_angle)

v_parton_kinetic = "~/Code/DMcode/"
v_parton_baryonic = "~/Code/DMcodeBaryon/"

SHIP_Energy = 400
miniboone_Energy= 8.9

local_directory=os.getcwd()

#write_lsnd(output_mode="particle_list",proddist="burmansmith",samplesize=1e6)
#subp.call(["./Source/main", "parameter_run.dat"])
#write_lsnd(proddist="particle_list",mdm=0.005,p_num_target=70)
#subp.call(["./Source/main", "parameter_run.dat"])

def prepare_parton(mA=1,energy=miniboone_Energy,file_path=v_parton_kinetic):
    with cd(file_path):
        subp.call(["./scalarDM",str(mA),"1",str(energy)])
        subp.call(["mv","dsigdpV.out","{}/Source/parton_V_p.dat".format(local_directory)])
        subp.call(["./scalarDM",str(mA),"2",str(energy)])
        subp.call(["mv","dsigdpV.out","{}/Source/parton_V_n.dat".format(local_directory)])
    textarr=["Source/parton_V_p.dat","Source/parton_V_n.dat"]
    for sfile in textarr:
        df = pd.read_csv(sfile, sep=" ", engine='python',skipinitialspace=True,header=None)
        arr= df.values.copy()
        arr2 = [[i, np.interp(i,arr[:,0],arr[:,1])] for i in np.arange(arr[0,0],arr[-1,0],0.1)]
        np.savetxt(sfile,arr2,delimiter=" ")

#write_miniboone(mdm=0.005,mv=0.4,proddist=["","","proton_brem"],prod_chan=["pi0_decay","eta_decay","V_decay"],partlistfile=["","",""])
