#include <vector>
#include "Random.h"
#include "partonsample.h"
#include <iostream>

using std::shared_ptr;
//No error checking here! It's awful!
parton_sample::parton_sample(std::string &pfilename, std::string &nfilename, double proton_n, double neutron_n){
	pmax_n=0;
	pmax_p=0;
	parse_V_dist(pfilename, V_p_dist, cross_section_p,pmax_p,mom_min_p,mom_max_p);
	parse_V_dist(nfilename, V_n_dist, cross_section_n,pmax_n,mom_min_n,mom_max_n);
	neutron_number=neutron_n;
	proton_number=proton_n;
}

void parton_sample::parse_V_dist(std::string &filename, shared_ptr<Linear_Interpolation> &V_dist, double &cross_section, double &pmax, double &mom_min, double &mom_max){
    std::ifstream instream(filename);
	if(!instream.is_open()){
		std::cerr << "parton_sample cannot open " << filename << std::endl;
		throw -1;
	}

	std::vector<double> v_dist;
    double in;
    instream >> mom_min;
    instream >> in;
    v_dist.push_back(in);
    pmax=in;
    double total = in;
    while(instream >> mom_max){
        instream >> in;
        v_dist.push_back(in);
        total+=in;
        if(in>pmax)
            pmax=in;
    }
   
    cross_section = total*(mom_max-mom_min)/v_dist.size();

    V_dist = shared_ptr<Linear_Interpolation>( new Linear_Interpolation(v_dist, mom_min, mom_max));
}

void parton_sample::sample_particle(Particle &part){
	double mom, theta, phi;
	sample_momentum(mom, theta, phi);
	part.ThreeMomentumPolar(mom, theta, phi);
}

void parton_sample::sample_momentum(double &pmom, double &theta, double &phi){
    theta=0; 
	phi=0; 
	pmom=gen_V_mom();
}

double parton_sample::gen_V_mom(){
    double mom;
	if(Random::Flat(0,1)<proton_number/(proton_number+neutron_number))
    	while(true){
        	mom = Random::Flat(mom_min_p, mom_max_p);
        	if(V_p_dist->Interpolate(mom)>Random::Flat(0,1)*pmax_p)
            	break;
    	}
	else
		while(true){
        	mom = Random::Flat(mom_min_n, mom_max_n);
        	if(V_n_dist->Interpolate(mom)>Random::Flat(0,1)*pmax_n)
            	break;
    	}
    return mom;
}

