#include "DMgenerator.h"
#include "branchingratios.h"
#include "decay.h"
#include "constants.h"
#include <cstdlib>
#include <iostream>
#include <fstream>

//std::ofstream datalog("piminuslog.dat", std::ios::out);

const double mass = 0.129;//From MacDonald et al. 1976. Need to implement a six percent spread.
//const double pi = 3.14159;
const int BURN_MAX = 1000;

void IsotropicParticleGen(Particle &part, double Energy){
    double theta = acos(Random::Flat(-1,1));
    double phi = Random::Flat(0,1)*2*pi;
    double mom;
    //datalog << "Isotropic Report\n";
    //part.report(datalog);
    if(Energy < part.m)
        mom=0;
    else
        mom = sqrt(pow(Energy,2)-pow(part.m,2));
                            
    part.ThreeMomentum(mom*sin(theta)*cos(phi),mom*sin(theta)*sin(phi),mom*cos(theta));
}

//Off-shell equations for this code are dodgy. I need some theory input on these branching ratios.

piminus_capture_gen::piminus_capture_gen(double MV, double MX, double kap, double alp){
    set_model_params(MV, MX, kap, alp);
    chan_name="Piminus_capture";
}

void piminus_capture_gen::burn_in(int runs){
    double s, theta;
    for(int i=0; i<runs; i++)
       sample_dist(s, theta);
}

void piminus_capture_gen::Evaluate_Branching_Ratio(){
    branchingratio = brmass_to_dm_dm(mass, mv, mx, kappa, alphaD);
    pmax = 0;//resetting pmax to zero, as new model parameters have been chosen.

    if(branchingratio>(brmasstoVgamma(mass,mv,mx,kappa,alphaD)*brV_to_dm_dm(mv, mx, kappa, alphaD)*1.3)){
        OFF_SHELL=true;
        if(mv*mv<=mass*mass && mv*mv >= 4*mx*mx)
            pmax =  d2brmass_to_dm_dm(mass, mv, mx, kappa, alphaD, mv*mv, pi/2);
        burn_in(BURN_MAX);
    }
    else
        OFF_SHELL=false;
    branchingratio*=0.5*0.5;//Photon production only occurs half the time according to Maxim.
}

void piminus_capture_gen::sample_dist(double& s, double& theta){
    while(true){
        s = Random::Flat(4*mx*mx,pow(mass,2));
        theta = Random::Flat(0,pi);
        if(d2brmass_to_dm_dm(mass, mv, mx, kappa, alphaD, s, theta)>pmax*Random::Flat(0,1)){
            if(d2brmass_to_dm_dm(mass, mv, mx, kappa, alphaD, s, theta)>pmax)
                pmax = d2brmass_to_dm_dm(mass, mv, mx, kappa, alphaD, s, theta);
            break;    
        }
    }
}

bool piminus_capture_gen::GenDM(std::list<Particle>& vec, std::function<double(Particle&)> det_int, Particle &darkphoton){
    double intersect1=0;
    double intersect2=0;

    //I don't think this is necessary, but it shouldn't change anything.
    darkphoton.Set_Mass(mv);
    darkphoton.name = "V";
    Particle darkmatter1(mx);
    darkmatter1.name = "DM";
    Particle darkmatter2(mx);
    darkmatter2.name = "DM";
    
    if(OFF_SHELL){
        double s, theta;//theta is the angle relative to the gamma in the V* rest frame
        sample_dist(s, theta);
        darkphoton.m = sqrt(s);
        IsotropicParticleGen(darkphoton,mass);
        Meson_Capture_Off_Shell(darkmatter1, darkmatter2, darkphoton, theta);
    }
    else{
        IsotropicParticleGen(darkphoton,mass);
        TwoBodyDecay(darkphoton, darkmatter1, darkmatter2);
    }
    //darkphoton.report(datalog);
    //darkmatter1.report(datalog);
    //darkmatter2.report(datalog);
    vec.push_back(darkphoton);
    if((intersect1=det_int(darkmatter1))>0 || (intersect2=det_int(darkmatter2))>0){
        if(intersect1>0){
            //datalog << "DM1 hit\n";
            vec.push_back(darkmatter1);
        }
        if(intersect2>0){
            //datalog << "DM2 hit\n";
            vec.push_back(darkmatter2);
        }
        //darkphoton.report(datalog);
        //darkmatter1.report(datalog);
        //darkmatter2.report(datalog);
        return true;
    }
    else
        return false;
}
