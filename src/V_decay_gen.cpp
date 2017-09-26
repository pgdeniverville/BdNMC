#include "DMgenerator.h"
#include "branchingratios.h"
#include "decay.h"
#include <cstdlib>
#include <iostream>
#include <fstream>
#include "constants.h"

using std::cout; using std::endl;

//std::ofstream logging("log.dat", std::ofstream::out);

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

bool V_decay_gen::GenDM(std::list<Particle>& vec, std::function<double(Particle&)> det_int, Particle& part){
    double intersect1=0;
    double intersect2=0;

    Particle darkphoton = part;
    darkphoton.name = "V";
    Particle darkmatter1(mx);
    darkmatter1.name = "DM";
    Particle darkmatter2(mx);
    darkmatter2.name = "DM";

	//darkphoton.report(logging);

    TwoBodyDecay(darkphoton, darkmatter1, darkmatter2);
  	//darkmatter1.report(logging);
	//darkmatter2.report(logging);
    


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

Do_Nothing_Gen::Do_Nothing_Gen(const std::string chan, const std::string part_name){
    Part_Name = std::string(part_name);
    chan_name = std::string(chan);
    branchingratio=1;
}

void Do_Nothing_Gen::Evaluate_Branching_Ratio(){
    branchingratio=1;
}


bool Do_Nothing_Gen::GenDM(std::list<Particle>& vec, std::function<double(Particle&)> det_int, Particle& part){
    part.name = Part_Name;
    
    if(det_int(part)>0){
        vec.push_back(part);
        return true;
    }

    return false;
}

//A general isotropic two body decay. Only cares about the kinematics!

Two_Body_Decay_Gen::Two_Body_Decay_Gen(double branching_ratio, double parent_mass, std::string part_name, Particle Daughter1, Particle Daughter2){
    if(parent_mass < Daughter1.m+Daughter2.m){
        std::cerr << "Parent_Mass is smaller than daughter masses, invalid decay!\n";
        throw -30;
    }
    Part_Name = std::string(part_name);
    daughter1 = Particle(Daughter1);
    daughter2 = Particle(Daughter2);
    branchingratio=branching_ratio;
}

bool Two_Body_Decay_Gen::GenDM(std::list<Particle>& vec, std::function<double(Particle&)> det_int, Particle& part){
    Particle parent = Particle(part);
    parent.name = Part_Name;
    vec.push_back(parent);

    TwoBodyDecay(parent, daughter1, daughter2);
    bool intersect=false;
    if(det_int(daughter1)>0){
        intersect=true;
        vec.push_back(daughter1);
    }
    if(det_int(daughter2)>0){
        intersect=true;
        vec.push_back(daughter2);
    }
    return intersect;
}
