#include "DMgenerator.h"
#include "branchingratios.h"
#include "decay.h"
#include <cstdlib>
#include <iostream>
#include <fstream>
#include "constants.h"

V_decay_gen::V_decay_gen(double MV, double MX, double kap, double alp, const std::string chan){
    set_model_params(MV, MX, kap, alp);
	chan_name=std::string(chan);
}

void V_decay_gen::Evaluate_Branching_Ratio(){
	if(chan_name=="V_decay_baryonic"){
        branchingratio= brVB_to_dm_dm(mv,mx,kappa,alphaD);
    }
    else{
        branchingratio = brV_to_dm_dm(mv, mx, kappa, alphaD);
    }
    OFF_SHELL = false;
}

bool V_decay_gen::GenDM(std::list<Particle>& vec, std::function<double(Particle)> det_int, Particle& part){
    double intersect1=0;
    double intersect2=0;

    Particle darkphoton = part;
    darkphoton.name = "V";
    Particle darkmatter1(mx);
    darkmatter1.name = "DM";
    Particle darkmatter2(mx);
    darkmatter2.name = "DM";



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

//Do_Nothing_Gen does not decay the received particle! This whole 
//file should be renamed something more general.

Do_Nothing_Gen::Do_Nothing_Gen(const std::string chan){
    chan_name = std::string(chan);
    branchingratio=1;
}

void Do_Nothing_Gen::Evaluate_Branching_Ratio(){
    return;
}

bool Do_Nothing_Gen::GenDM(std::list<Particle>& vec, std::function<double(Particle)> det_int, Particle& part){
    vec.push_back(part);

    if(det_int(part)>0)
        return true;
}
