#include "Random.h"
#include "Proton_Brem_Distribution.h"
#include <math.h>
#include <iostream>
#include <functional>
#include "Integrator.h"
#include "branchingratios.h"
#include "constants.h"


using std::bind;
using namespace std::placeholders;
const double mp = MASS_PROTON;


Proton_Brem_Distribution::Proton_Brem_Distribution(double Beam_E, double epsilon, double mA, double ptmax, double zmax, double zmin, double ptmin){
	Beam_Energy=Beam_E; kappa=epsilon; PTMIN=ptmin; PTMAX=ptmax; ZMAX = zmax; ZMIN = zmin; MA=mA;
	sppM = pow(2*mp+Mpp,2);
	calc_V_prod_rate();
}

//Total proton-proton scattering cross section
double Proton_Brem_Distribution::sigmapp(double s){
	return Hpp*pow(log(s/sppM),2)+Ppp+R1pp*pow(s/sppM,-eta1pp)-R2pp*pow(s/sppM,-eta2pp);
}

//Very basic version of form factor. Maybe should move form factors from DMNscattering.cpp into their own file and include them? Probably do this when I do the big constants reorganization.

double Proton_Brem_Distribution::F_1_proton(double q2){
	return pow(1+q2/mD2,-2);
}

double Proton_Brem_Distribution::d2N_proton_brem_to_V(double z, double pt2){
	return pow(F_1_proton(MA*MA),2)*sigmapp(2*mp*(Beam_Energy-sqrt(MA*MA+pt2+z*z*(pow(Beam_Energy,2)-mp*mp))))/sigmapp(2*mp*Beam_Energy)*wpp(z,pt2, MA, kappa);
}

void Proton_Brem_Distribution::set_fit_parameters(production_channel &par){
	par.query_dist_param("rD",rD);
	par.query_dist_param("mD2",mD2);
	par.query_dist_param("H",Hpp);
	par.query_dist_param("M",Mpp);
	par.query_dist_param("eta1",eta1pp);
	par.query_dist_param("eta2",eta2pp);
	par.query_dist_param("R1pp",R1pp);
	par.query_dist_param("R2pp",R2pp);
	par.query_dist_param("Ppp",Ppp);
	sppM = pow(2*mp+Mpp,2);
	max_prod = d2N_proton_brem_to_V(ZMIN,PTMIN);
}

void Proton_Brem_Distribution::calc_V_prod_rate(){
	std::function<double(double, double)> func = std::bind(&Proton_Brem_Distribution::d2N_proton_brem_to_V,this,_1,_2);
	vprodrate = SimpsonCubature(func,ZMIN,ZMAX,100,PTMIN,PTMAX,100);//Will need to tweak this. Hopefully come up with a more general algorithm at some point.
	max_prod=d2N_proton_brem_to_V(ZMIN,PTMIN);
};

void Proton_Brem_Distribution::sample_particle(Particle &part){
	double mom, theta, phi;
	sample_momentum(mom, theta, phi);
	part.ThreeMomentumPolar(mom, theta, phi);
}

void Proton_Brem_Distribution::sample_momentum(double &pmom, double &theta, double &phi){
	double z, pt, hold;
	while(true){
		z = Random::Flat(ZMIN,ZMAX);
		pt = Random::Flat(PTMIN,PTMAX);
		if((hold=d2N_proton_brem_to_V(z,pt))>Random::Flat()*max_prod){
			if(max_prod<hold)
				max_prod=hold;
			break;
		}
	}
	phi = Random::Flat(0,2*pi);
	double pz = sqrt(pow(Beam_Energy,2)-mp*mp)*z;
	theta=atan2(pt,pz);
	pmom = sqrt(pz*pz+pt*pt);
}
