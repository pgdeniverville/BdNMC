from par_writerbak import *
import numpy as np
import sys

import time

_pion_decay=1
_eta_decay=2
_rho_decay=3
_parton=4
_brem=5
_piminus_cap=6

def coherent_eval(mass_arr,det_switch="LAr",coherent="true",channels={_pion_decay,_piminus_cap},signal_channel="NCE_nucleon",alpha_D=0.5,sumlog="Events/coherent_y.dat"):
    MV=mass_arr[0]
    MX=mass_arr[1]
    t0=time.time()
    parfile="parameter_run_{0}_{1}.dat".format(str(MV),str(MX))

    proddist = []
    prodchan = []
    partlistfile = [    ]
    executing=False
    if signal_channel=="NCE_nucleon_baryonic":
        if MX/1000.0<mpi0/2.0 and MV<600.0 and _pion_decay in channels:
            proddist.append("particle_list")
            prodchan.append("pi0_decay_baryonic")
            partlistfile.append("data/particle_list.dat")
            executing = True

    else:
        if MX/1000.0<mpi0/2.0 and MV<600.0 and _pion_decay in channels:
            proddist.append("particle_list")
            prodchan.append("pi0_decay")
            partlistfile.append("data/particle_list.dat")
            executing = True
        if MX/1000.0<0.129/2.0 and MV<600 and MV>2*MX and _piminus_cap in channels:
            proddist.append("")
            prodchan.append("piminus_capture")
            partlistfile.append("")
            executing = True

    if not executing:
        return

    if det_switch == "LAr":
        det=coherent_detector_LAr
    elif det_switch == "NaI":
        det=coherent_detector_NaI

    if signal_channel=="NCE_nucleon":
        write_coherent(mdm=MX/1000.0,mv=MV/1000.0,proddist=proddist,\
                prod_chan=prodchan,\
                partlistfile=partlistfile,outfile=parfile,samplesize=100,\
                signal_chan=signal_channel,sumlog=sumlog,min_scatter_energy=15e-6,\
                max_scatter_energy=0.2,alpha_D=alpha_D,eps=1e-3,efficiency=0.5,\
                det=det,output_mode="summary",\
                outlog="Events/coherent_events.dat",coherent=coherent)
    elif signal_channel=="NCE_nucleon_baryonic":
        write_coherent(mdm=MX/1000.0,mv=MV/1000.0,proddist=proddist,\
                prod_chan=prodchan,\
                partlistfile=partlistfile,outfile=parfile,samplesize=200,\
                signal_chan=signal_channel,sumlog=sumlog,\
                min_scatter_energy=15e-6,\
                max_scatter_energy=0.2,alpha_D=1e-4,eps=0.0,efficiency=0.5,\
                det=det,output_mode="summary",\
                outlog="Events/coherent_events.dat",coherent=coherent)

    subp.call(["./build/main", parfile])
    t1 = time.time()
    print("\ntime={}\n".format(t1-t0))
    t0 = time.time()
    subp.call(["rm", parfile])

def execute_coherent(genlist=True):
    if genlist:
        write_coherent(prod_chan=["pi0_decay"],proddist=["burmansmith"],samplesize=1e6,output_mode="particle_list",partlistfile=["data/particle_list_coherent.dat"],p_num_target=80)
        subp.call(["./build/main","parameter_run.dat"])
    #vmassarr=[i for i in range(11,30,2)]+[i for i in range(30,130,10)]+[129,131,132,134,136,138]+[i for i in range(140,700,10)]
    #chimassarr=[5]
    #massarr=[[MV,MX] for MV in vmassarr for MX in chimassarr]
    vmassarr=[3,6,9]
    massarr=[[MV,MV/3.0] for MV in vmassarr]
    #massarr=massarr+massarr2
    for marr in massarr:
        coherent_eval(marr,signal_channel="NCE_nucleon",sumlog="Events/coherent_LAr_test1.dat",coherent="true",det_switch="LAr")
        #coherent_eval(marr,signal_channel="NCE_nucleon",sumlog="Events/coherent_NaI2.dat",coherent="true",det_switch="NaI")
        coherent_eval(marr,signal_channel="NCE_nucleon_baryonic",sumlog="Events/coherent_LAr_test1.dat",coherent="true",det_switch="LAr")
        #coherent_eval(marr,signal_channel="NCE_nucleon_baryonic",sumlog="Events/coherent_NaI2.dat",coherent="true",det_switch="NaI")


execute_coherent(genlist=False)

