//Declared in DMgenerator.h

#include "DMgenerator.h"
#include "Integrator.h"
#include "Random.h"
#include "Kinematics.h"
#include <cmath>
#include "constants.h"

using std::endl;
using std::cout;

Generic_Bremsstrahlung::Generic_Bremsstrahlung(Particle& brem_product, double Beam_energy, double Ptmin, double Ptmax, double Zmin, double Zmax, std::string prodchoice, std::function<double(double, double)> dist){

	product = brem_product;
	beam_energy = Beam_energy;
	ZMIN=Zmin;
	ZMAX=Zmax;
	PTMIN=Ptmin;
	PTMAX=Ptmax;
	split_dist = dist;

	Set_Channel_Name(prodchoice);

	calc_prod_rate();
}

void Generic_Bremsstrahlung::Burn_In(int iterations){
	for(int i =0; i < iterations; i++){
		sample_particle(product);
	}
}

void Generic_Bremsstrahlung::Toggle_Decay(std::shared_ptr<DMGenerator> decay){
	unstable = true;
	decay_gen = decay;

	decay_gen->record_parent = false;
}

void Generic_Bremsstrahlung::calc_prod_rate(){
	if(ZMIN>ZMAX){
		prod_rate=0;
	}
	else{
		prod_rate = SimpsonCubature_adapt(split_dist,ZMIN,ZMAX,100,PTMIN,PTMAX,100,0.01);
		//prod_rate = RandomIntegrate2(split_dist,ZMIN,ZMAX,PTMIN,PTMAX,1e6);
	}
	pmax = split_dist(ZMIN,PTMIN);
}

bool Generic_Bremsstrahlung::GenDM(std::list<Particle>& vec, std::function<double(Particle&)> det_int, Particle& part){

	//Generate a bremsstrahlung product
	sample_particle(product);

	//Rotate it so its z-axis is parallel to the trajectory of its parent.
	product.Rotate(part);

	if(unstable){
		//This particle decays or does something else transformative, we're looking for its decay products in the detector.
		vec.push_back(product);
		return decay_gen->GenDM(vec, det_int, product);
	}
	else{
		//This particle does not decay, it needs to intersect the detector
		if(det_int(product)>0){
			vec.push_back(product);
			return true;
		}
	}
	return false;
}

void Generic_Bremsstrahlung::sample_particle(Particle &part){
	double z, pt, hold;
	while(true){
		z = Random::Flat(ZMIN,ZMAX);
		pt = Random::Flat(PTMIN,PTMAX);
		if((hold=split_dist(z,pt))>Random::Flat()*pmax){
			if(pmax<hold)
				pmax=hold;
			break;
		}
	}
	double phi = Random::Flat(0,2*pi);
	double pz = sqrt(pow(beam_energy,2)-pow(MASS_PROTON,2))*z;
	double theta=atan2(pt,pz);
	double pmom = sqrt(pz*pz+pt*pt);

	part.ThreeMomentumPolar(pmom, theta, phi);
}