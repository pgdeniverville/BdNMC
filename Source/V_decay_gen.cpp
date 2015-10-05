#include "DMgenerator.h"
#include "branchingratios.h"
#include "decay.h"
#include <cstdlib>
#include <iostream>
#include <fstream>

const double pi = 3.14159;

proton_brem::proton_brem(double MV, double MX, double kap, double alp, const std::string chan){
    set_model_params(MV, MX, kap, alp);
	chan_name=std::string(chan);
}

void proton_brem::Evaluate_Branching_Ratio(){
	branchingratio = brV_to_dm_dm(mv, mx, kappa, alphaD);
	OFF_SHELL = false;
}

bool proton_brem::GenDM(std::list<Particle>& vec, std::function<double(Particle)> det_int, std::shared_ptr<Particle_Generator> V_prod){
    double intersect1=0;
    double intersect2=0;

    Particle darkphoton(mv);
    darkphoton.name = "V";
    Particle darkmatter1(mx);
    darkmatter1.name = "DM";
    Particle darkmatter2(mx);
    darkmatter2.name = "DM";

    V_prod->Generate_Particle(darkphoton);
	//darkphoton.report(std::cout);

    TwoBodyDecay(darkphoton, darkmatter1, darkmatter2);
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
