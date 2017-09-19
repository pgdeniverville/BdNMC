#include "DMgenerator.h"
#include "branchingratios.h"
#include "decay.h"
#include "constants.h"
#include <cstdlib>
#include <iostream>
#include <fstream>
#include <exception>
//std::ofstream datalogpion("pionlog.dat",std::ios::out);

/*
 *
 *PI0 DECAY GEN
 *
 */

using std::cout;
using std::cerr;
using std::endl;
const int BURN_MAX = 1000;

pion_decay_gen::pion_decay_gen(double MV, double MX, double kap, double alp){
    set_model_params(MV, MX, kap, alp);
    chan_name="pi0_decay";
}

void pion_decay_gen::burn_in(int runs){
    double s, theta;
    for(int i=0; i<runs; i++)
       sample_dist(s, theta);
}

void pion_decay_gen::Evaluate_Branching_Ratio(){
    if(2*mx>mpi0){
        branchingratio = 0;
        cerr << chan_name << " forbidden by energy conservation.";
		throw std::domain_error("mx>mpi0/2, out of domain"); 
    }

    branchingratio = brpi0_to_gamma_dm_dm(mv, mx, kappa, alphaD);
    pmax = 0;//resetting pmax to zero, as new model parameters have been chosen.
    
    if(branchingratio<brpi0toVgamma(mv,mx,kappa,alphaD)*brV_to_dm_dm(mv, mx, kappa, alphaD))
            std::cerr << "Waring: Off-Shell branching ratio calculated at less than on-shell branching-ratio.\nOff-Shell Branching Ratio: " << branchingratio  << " On-Shell: " << brpi0toVgamma(mv,mx,kappa,alphaD)*brV_to_dm_dm(mv, mx, kappa, alphaD) << std::endl;

    if(branchingratio>brpi0toVgamma(mv,mx,kappa,alphaD)*brV_to_dm_dm(mv, mx, kappa, alphaD)*1.30){
                OFF_SHELL=true;
        if(mv*mv<=mpi0*mpi0 && mv*mv >= 4*mx*mx)
            pmax =  d2brpi0_to_gamma_dm_dm(mv, mx, kappa, alphaD, mv*mv, pi/2);
		else
        	burn_in(BURN_MAX);//skip if on shell available?
    }
    else
        OFF_SHELL=false;
}

void pion_decay_gen::sample_dist(double& s, double& theta){
    while(true){
        s = Random::Flat(4*mx*mx,pow(mpi0,2));
        theta = Random::Flat(0,pi);
		double hold;
        if((hold=d2brpi0_to_gamma_dm_dm(mv, mx, kappa, alphaD, s, theta))>pmax*Random::Flat(0,1)){
            if(hold>pmax)
                pmax = hold;
            break;    
        }
    }
}

bool pion_decay_gen::GenDM(std::list<Particle>& vec, std::function<double(Particle&)> det_int, Particle& part){
    double intersect1=0;
    double intersect2=0;
    Particle meson = part;
    meson.name = "pion";

    Particle darkphoton(mv);//In off-shell case this is changed to sqrt(s)
    darkphoton.name = "V";
    Particle darkmatter1(mx);
    darkmatter1.name = "DM";
    Particle darkmatter2(mx);
    darkmatter2.name = "DM";

    if(OFF_SHELL&&!((mx*2<=mv)&&(mv<meson.m)&&(2*mx>=0.97*mv))){
		double s, theta;//theta is the angle relative to the gamma in the V* rest frame
        sample_dist(s, theta);
        darkphoton.Set_Mass(sqrt(s));
        DecayDM_Off_Shell(darkmatter1, darkmatter2, darkphoton, meson, theta);
    }
    else
        DecayDM(darkmatter1, darkmatter2, darkphoton, meson);
/*
    darkmatter1.report(cout);
    darkmatter2.report(cout);
    darkphoton.report(cout);
    meson.report(cout);
*/
    vec.push_back(meson);
 	intersect1=det_int(darkmatter1);
    intersect2=det_int(darkmatter2);

    if((intersect1)>0 || (intersect2)>0){
        vec.push_back(darkphoton);
        if(intersect1>0){
            vec.push_back(darkmatter1);
        }
        if(intersect2>0){
            vec.push_back(darkmatter2);
        }
        return true;
    }  
    else
        return false;
}
/*
 *
 *ETA DECAY GEN
 *
 */

eta_decay_gen::eta_decay_gen(double MV, double MX, double kap, double alp){
    set_model_params(MV, MX, kap, alp);
    chan_name="eta_decay";
}


void eta_decay_gen::burn_in(int runs){
    double s, theta;
    for(int i=0; i<runs; i++)
       sample_dist(s, theta);
}

void eta_decay_gen::Evaluate_Branching_Ratio(){
    if(2*mx>meta){
        branchingratio = 0;
        cerr << chan_name << " forbidden by energy conservation.";
        return;
           
    }
    
    branchingratio = breta_to_gamma_dm_dm(mv, mx, kappa, alphaD);
    //cout << branchingratio << endl;
	double onshellbr = bretatoVgamma(mv,mx,kappa,alphaD)*brV_to_dm_dm(mv, mx, kappa, alphaD);
    //cout << onshellbr << endl;
	pmax = 0;//resetting pmax to zero, as new model parameters have been chosen.
    
    if(branchingratio<onshellbr)
            std::cerr << "Off-Shell Branching Ratio: " << branchingratio  << " On-Shell: " << bretatoVgamma(mv,mx,kappa,alphaD)*brV_to_dm_dm(mv, mx, kappa, alphaD) << std::endl;

    if(branchingratio>(onshellbr*1.30)){
        OFF_SHELL=true;
        if(mv<=meta && mv > 2*mx)
            pmax =  d2breta_to_gamma_dm_dm(mv, mx, kappa, alphaD, mv*mv, pi/2);
        burn_in(BURN_MAX);//skip if on-shell available?
    }
    else
        OFF_SHELL=false;
}

void eta_decay_gen::sample_dist(double& s, double& theta){
    while(true){
        s = Random::Flat(4*mx*mx,pow(meta,2));
        theta = Random::Flat(0,pi);
        double hold;
		if((hold=d2breta_to_gamma_dm_dm(mv, mx, kappa, alphaD, s, theta))>pmax*Random::Flat(0,1)){
            if(hold>pmax)
                pmax = hold;
            break;    
        }
    }
}

bool eta_decay_gen::GenDM(std::list<Particle>& vec, std::function<double(Particle&)> det_int, Particle& part){
    double intersect1=0;
    double intersect2=0;
    
    Particle meson = part;
    meson.name = "eta";

    Particle darkphoton(mv);//In off-shell case this is changed to sqrt(s)
    darkphoton.name = "V";
    Particle darkmatter1(mx);
    darkmatter1.name = "DM";
    Particle darkmatter2(mx);
    darkmatter2.name = "DM";

    if(OFF_SHELL&&!((mx*2<=mv)&&(mv<meson.m)&&(2*mx>=0.97*mv))){
        double s, theta;//theta is the angle relative to the gamma in the V* rest frame
        sample_dist(s, theta);
        darkphoton.Set_Mass(sqrt(s));
        DecayDM_Off_Shell(darkmatter1, darkmatter2, darkphoton, meson, theta);
    }
    else{
        DecayDM(darkmatter1, darkmatter2, darkphoton, meson);
    }

    intersect1=det_int(darkmatter1);
    intersect2=det_int(darkmatter2);
/*
    meson.report(std::cout);
    darkphoton.report(std::cout);
    darkmatter1.report(std::cout);
    darkmatter2.report(std::cout);
*/

    vec.push_back(meson);
    if((intersect1)>0 || (intersect2)>0){
		vec.push_back(darkphoton);
        if(intersect1>0){
            vec.push_back(darkmatter1);
        }
        if(intersect2>0){
            vec.push_back(darkmatter2);
        }
        return true;
    }	
    else
        return false;
}

/*
 *
 *RHO DECAY GEN
 *
 */



rho_decay_gen::rho_decay_gen(double MV, double MX, double kap, double alp){
    if(MV<=2*MX)
        std::cerr << "invalid model parameters. M_A <= 2*M_DM\n";
    set_model_params(MV, MX, kap, alp);
    chan_name="rho_decay";
    OFF_SHELL = true;
}

void rho_decay_gen::Evaluate_Branching_Ratio(){
    branchingratio = brrho_to_V(mv, mx, kappa,alphaD)*brV_to_dm_dm(mv,mx,kappa,alphaD);
}

bool rho_decay_gen::GenDM(std::list<Particle>& vec, std::function<double(Particle&)> det_int, Particle& part ){
    double intersect1=0;
    double intersect2=0;
    
    Particle meson = part;
    meson.name = "rho";

    Particle darkphoton(mv);//In off-shell case this is changed to sqrt(s)
    darkphoton.name = "V";
    Particle darkmatter1(mx);
    darkmatter1.name = "DM";
    Particle darkmatter2(mx);
    darkmatter2.name = "DM";
    
    darkphoton.ThreeMomentum(meson.px, meson.py, meson.pz);//Omega oscillates into V

    TwoBodyDecay(darkphoton, darkmatter1, darkmatter2);
    intersect1=det_int(darkmatter1);
    intersect2=det_int(darkmatter2);
    vec.push_back(meson);
    if((intersect1)>0 || (intersect2)>0){
        vec.push_back(darkphoton);
        if(intersect1>0){
            vec.push_back(darkmatter1);
        }
        if(intersect2>0){
            vec.push_back(darkmatter2);
        }
        return true;
    }
    else
        return false;
}

/*
 *
 *OMEGA DECAY GEN
 *
 */


omega_decay_gen::omega_decay_gen(double MV, double MX, double kap, double alp){
    if(MV<=2*MX)
        std::cerr << "invalid model parameters. M_A <= 2*M_DM\n";
    set_model_params(MV, MX, kap, alp);
    chan_name="omega_decay";
    OFF_SHELL = true;
}

void omega_decay_gen::Evaluate_Branching_Ratio(){
    branchingratio = bromega_to_V(mv, mx, kappa,alphaD)*brV_to_dm_dm(mv,mx,kappa,alphaD);
}

bool omega_decay_gen::GenDM(std::list<Particle>& vec, std::function<double(Particle&)> det_int, Particle& part ){
    double intersect1=0;
    double intersect2=0;
    
    Particle meson = part;
    meson.name = "omega";

    Particle darkphoton(mv);//In off-shell case this is changed to sqrt(s)
    darkphoton.name = "V";
    Particle darkmatter1(mx);
    darkmatter1.name = "DM";
    Particle darkmatter2(mx);
    darkmatter2.name = "DM";
    
    darkphoton.ThreeMomentum(meson.px, meson.py, meson.pz);//Omega oscillates into V

    TwoBodyDecay(darkphoton, darkmatter1, darkmatter2);
    intersect1=det_int(darkmatter1);
    intersect2=det_int(darkmatter2);
    vec.push_back(meson);
    if((intersect1)>0 || (intersect2)>0){
        vec.push_back(darkphoton);
        if(intersect1>0){
            vec.push_back(darkmatter1);
        }
        if(intersect2>0){
            vec.push_back(darkmatter2);
        }
        return true;
    }
    else
        return false;
}

/*
 *
 *PHI DECAY GEN
 *
 */


phi_decay_gen::phi_decay_gen(double MV, double MX, double kap, double alp){
    if(MV<=2*MX)
        std::cerr << "invalid model parameters. M_A <= 2*M_DM\n";
    set_model_params(MV, MX, kap, alp);
    chan_name="phi_decay";
    OFF_SHELL = true;
}

void phi_decay_gen::Evaluate_Branching_Ratio(){
    branchingratio = brphi_to_V(mv, mx, kappa,alphaD)*brV_to_dm_dm(mv,mx,kappa,alphaD);
}

bool phi_decay_gen::GenDM(std::list<Particle>& vec, std::function<double(Particle&)> det_int, Particle& part){
    double intersect1=0;
    double intersect2=0;
    
    Particle meson = part;
    meson.name = "phi";

    Particle darkphoton(mv);//In off-shell case this is changed to sqrt(s)
    darkphoton.name = "V";
    Particle darkmatter1(mx);
    darkmatter1.name = "DM";
    Particle darkmatter2(mx);
    darkmatter2.name = "DM";
    
    darkphoton.ThreeMomentum(meson.px, meson.py, meson.pz);//Omega oscillates into V

    TwoBodyDecay(darkphoton, darkmatter1, darkmatter2);
    intersect1=det_int(darkmatter1);
    intersect2=det_int(darkmatter2);
    vec.push_back(meson);
    if((intersect1)>0 || (intersect2)>0){
        vec.push_back(darkphoton);
        if(intersect1>0){
            vec.push_back(darkmatter1);
        }
        if(intersect2>0){
            vec.push_back(darkmatter2);
        }
        return true;
    }
    else
        return false;
}
