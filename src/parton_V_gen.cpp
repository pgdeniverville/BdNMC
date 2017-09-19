#include "DMgenerator.h"
#include "branchingratios.h"
#include "decay.h"
#include <cstdlib>
#include <iostream>
#include <fstream>
#include "constants.h"
#include "Integrator.h"

parton_V_gen::parton_V_gen(double MV, double MX, double kap, double alp, const std::string chan){
    set_model_params(MV, MX, kap, alp);
	chan_name=std::string(chan);
}

void parton_V_gen::Evaluate_Branching_Ratio(){
    if(chan_name=="parton_production_baryonic")
		branchingratio = brVB_to_dm_dm(mv, mx, kappa, alphaD);
	else
		branchingratio = brV_to_dm_dm(mv, mx, kappa, alphaD);
	OFF_SHELL = false;
}

bool parton_V_gen::GenDM(std::list<Particle>& vec, std::function<double(Particle&)> det_int, Particle &part){
    double intersect1=0;
    double intersect2=0;

    Particle darkphoton = part;
    darkphoton.name = "V";
    Particle darkmatter1(mx);
    darkmatter1.name = "DM";
    Particle darkmatter2(mx);
    darkmatter2.name = "DM";

	//darkphoton.report(std::cout);
	//std::cout << darkphoton.m << std::endl;
	//darkmatter1.report(std::cout);
	//darkmatter2.report(std::cout);


	//darkphoton.report(std::cout);
    double thetad;

    while(true){
        thetad=Random::Flat(0,pi/2.0);
        if((1-pow(cos(thetad),2))>Random::Flat(0,1))
            break;
    }
	//std::cout << thetad << std::endl;
    TwoBodyDecay(darkphoton, darkmatter1, darkmatter2, thetad);
  	//darkmatter1.report(std::cout);
	//darkmatter2.report(std::cout);

    vec.push_back(darkphoton);
    if((intersect1=det_int(darkmatter1))>0 || (intersect2=det_int(darkmatter2))>0){
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
