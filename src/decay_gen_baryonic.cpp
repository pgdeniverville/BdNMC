#include "DMgenerator.h"
#include "constants.h"
#include "branchingratios.h"
#include "decay.h"
#include <cstdlib>
#include <iostream>
#include <fstream>

/*
 *
 *Pion
 *
 */

const int BURN_MAX = 1000;

pion_decay_gen_baryonic::pion_decay_gen_baryonic(double MV, double MX, double kap, double alp){
    if(MV<=2*MX)
        std::cerr << "invalid model parameters. M_A <= 2*M_DM\n";
    set_model_params(MV, MX, kap, alp);
    chan_name="Pion_decay";
}

void pion_decay_gen_baryonic::burn_in(int runs){
    double s, theta;
    for(int i=0; i<runs; i++)
       sample_dist(s, theta);
}

void pion_decay_gen_baryonic::Evaluate_Branching_Ratio(){
    branchingratio = brpi0_to_gamma_dm_dm_baryonic(mv, mx, kappa, alphaD);
    pmax = 0;//resetting pmax to zero, as new model parameters have been chosen.
    double brtemp; 
    if(branchingratio<(brtemp=brpi0toVBgamma(mv,mx,kappa,alphaD)*brVB_to_dm_dm(mv, mx, kappa, alphaD)))
            std::cerr << "Branching ratios might not be calculated correctly!\nOff-Shell Branching Ratio: " << branchingratio  << " On-Shell: " << brtemp << std::endl;

    if(branchingratio>brtemp*1.30){
                OFF_SHELL=true;
        if(mv*mv<=mpi0*mpi0 && mv*mv >= 4*mx*mx)
            pmax =  d2brpi0_to_gamma_dm_dm_baryonic(mv, mx, kappa, alphaD, mv*mv, pi/2);
        burn_in(BURN_MAX);//skip if on shell available?
    }
    else
        OFF_SHELL=false;
}

void pion_decay_gen_baryonic::sample_dist(double& s, double& theta){
    while(true){
        s = Random::Flat(4*mx*mx,pow(mpi0,2));
        theta = Random::Flat(0,pi);
        if(d2brpi0_to_gamma_dm_dm_baryonic(mv, mx, kappa, alphaD, s, theta)>pmax*Random::Flat(0,1)){
            if(d2brpi0_to_gamma_dm_dm_baryonic(mv, mx, kappa, alphaD, s, theta)>pmax)
                pmax = d2brpi0_to_gamma_dm_dm_baryonic(mv, mx, kappa, alphaD, s, theta);
            break;    
        }
    }
}

bool pion_decay_gen_baryonic::GenDM(std::list<Particle>& vec, std::function<double(Particle&)> det_int, Particle& part){
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

    if(OFF_SHELL){
        double s, theta;//theta is the angle relative to the gamma in the V* rest frame
        sample_dist(s, theta);
        darkphoton.Set_Mass(sqrt(s));
        DecayDM_Off_Shell(darkmatter1, darkmatter2, darkphoton, meson, theta);
    }
    else
        DecayDM(darkmatter1, darkmatter2, darkphoton, meson);
/*
    darkmatter1.report(datalogpion);
    darkmatter2.report(datalogpion);
    meson.report(datalogpion);
*/
    vec.push_back(meson);
    if((intersect1=det_int(darkmatter1))>0 || (intersect2=det_int(darkmatter2))>0){
        vec.push_back(darkphoton);
        if(intersect1>0){
            vec.push_back(darkmatter1);
        }
        if(intersect2>0){
            vec.push_back(darkmatter2);
        }
        //datalog << "Event_Start" << endl;
        
        return true;
    }
    else
        return false;
}

/*
 *
 *Eta Decay Gen
 *
 */


eta_decay_gen_baryonic::eta_decay_gen_baryonic(double MV, double MX, double kap, double alp){
    set_model_params(MV, MX, kap, alp);
    chan_name="eta_decay";
}

void eta_decay_gen_baryonic::burn_in(int runs){
    double s, theta;
    for(int i=0; i<runs; i++)
       sample_dist(s, theta);
}

void eta_decay_gen_baryonic::Evaluate_Branching_Ratio(){
    branchingratio = breta_to_gamma_dm_dm_baryonic(mv, mx, kappa, alphaD);
    double onshellbr = bretatoVBgamma(mv,mx,kappa,alphaD)*brVB_to_dm_dm(mv, mx, kappa, alphaD);
    pmax = 0;//resetting pmax to zero, as new model parameters have been chosen.
    
    if(branchingratio<onshellbr)
            std::cerr << "Off-Shell Branching Ratio: " << branchingratio  << " On-Shell: " << bretatoVBgamma(mv,mx,kappa,alphaD)*brVB_to_dm_dm(mv, mx, kappa, alphaD) << std::endl;

    if(branchingratio>onshellbr*1.30){
                OFF_SHELL=true;
        if(mv*mv<=meta*meta && mv*mv >= 4*mx*mx)
            pmax =  d2breta_to_gamma_dm_dm_baryonic(mv, mx, kappa, alphaD, mv*mv, pi/2);
        burn_in(BURN_MAX);//skip if on shell available?
    }
    else
        OFF_SHELL=false;
}

void eta_decay_gen_baryonic::sample_dist(double& s, double& theta){
    while(true){
        s = Random::Flat(4*mx*mx,pow(meta,2));
        theta = Random::Flat(0,pi);
        if(d2breta_to_gamma_dm_dm_baryonic(mv, mx, kappa, alphaD, s, theta)>pmax*Random::Flat(0,1)){
            if(d2breta_to_gamma_dm_dm_baryonic(mv, mx, kappa, alphaD, s, theta)>pmax)
                pmax = d2breta_to_gamma_dm_dm_baryonic(mv, mx, kappa, alphaD, s, theta);
            break;    
        }
    }
}

bool eta_decay_gen_baryonic::GenDM(std::list<Particle>& vec, std::function<double(Particle&)> det_int, Particle& part){
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

    if(OFF_SHELL){
        double s, theta;//theta is the angle relative to the gamma in the V* rest frame
        sample_dist(s, theta);
        darkphoton.Set_Mass(sqrt(s));
        DecayDM_Off_Shell(darkmatter1, darkmatter2, darkphoton, meson, theta);
    }
    else
        DecayDM(darkmatter1, darkmatter2, darkphoton, meson);

    vec.push_back(meson);
    if((intersect1=det_int(darkmatter1))>0 || (intersect2=det_int(darkmatter2))>0){
        vec.push_back(darkphoton);
        if(intersect1>0){
            vec.push_back(darkmatter1);
        }
        if(intersect2>0){
            vec.push_back(darkmatter2);
        }
        //datalog << "Event_Start" << endl;
        
        return true;
    }
    else
        return false;
}
/*
 *
 *Omega 
 *
 */

omega_decay_gen_baryonic::omega_decay_gen_baryonic(double MV, double MX, double kap, double alp){
    if(MV<=2*MX)
        std::cerr << "invalid model parameters. M_A <= 2*M_DM\n";
    set_model_params(MV, MX, kap, alp);
    chan_name="omega_decay";
    OFF_SHELL = true;
}

void omega_decay_gen_baryonic::Evaluate_Branching_Ratio(){
    branchingratio = bromega_to_Vb(mv, mx, kappa,alphaD);
}

bool omega_decay_gen_baryonic::GenDM(std::list<Particle>& vec, std::function<double(Particle&)> det_int, Particle& part ){
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

phi_decay_gen_baryonic::phi_decay_gen_baryonic(double MV, double MX, double kap, double alp){
    if(MV<=2*MX)
        std::cerr << "invalid model parameters. M_A <= 2*M_DM\n";
    set_model_params(MV, MX, kap, alp);
    chan_name="phi_decay";
    OFF_SHELL = true;
}

void phi_decay_gen_baryonic::Evaluate_Branching_Ratio(){
    branchingratio = brphi_to_Vb(mv, mx, kappa,alphaD)*brVB_to_dm_dm(mv,mx,kappa,alphaD);
}

bool phi_decay_gen_baryonic::GenDM(std::list<Particle>& vec, std::function<double(Particle&)> det_int, Particle& part){
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
