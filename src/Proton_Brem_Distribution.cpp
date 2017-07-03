#include "Random.h"
#include "Proton_Brem_Distribution.h"
#include <math.h>
#include <iostream>
#include <functional>
#include "Integrator.h"
#include "branchingratios.h"
#include "constants.h"

#include <complex>

#include <iostream>

using std::complex;

using std::bind;
using std::cout;
using std::endl;
using namespace std::placeholders;
const double mp = MASS_PROTON;


complex<double> F_1_proton_baryonic(double q2, double kappa, double alD);
complex<double> F_1_proton(double);

double gu(double alD, double kappa){
    return sqrt(4*pi*alD)/3-2*kappa*sqrt(4*pi*alphaEM)/3;
}

double gd(double alD, double kappa){
    return sqrt(4*pi*alD)/3+kappa*sqrt(4*pi*alphaEM)/3;
}

Proton_Brem_Distribution::Proton_Brem_Distribution(double Beam_E, double epsilon, double mA, double ptmax, double zmax, double zmin, double alphaD,  std::string &mode, double ptmin){
	Beam_Energy=Beam_E; kappa=epsilon; PTMIN=ptmin; PTMAX=ptmax; ZMAX = zmax; ZMIN = zmin; MA=mA; model=mode; alpha_D = alphaD;
	sppM = pow(2*mp+Mpp,2);
	set_mass(mA);
	calc_V_prod_rate();
}

//Total proton-proton scattering cross section
double Proton_Brem_Distribution::sigmapp(double s){
	return Hpp*pow(log(s/sppM),2)+Ppp+R1pp*pow(s/sppM,-eta1pp)-R2pp*pow(s/sppM,-eta2pp);
}

//Very basic version of form factor. Maybe should move form factors from DMNscattering.cpp into their own file and include them? Probably do this when I do the big constants reorganization.
/*
double Proton_Brem_Distribution::F_1_proton(double q2){
	return pow(1+q2/mD2,-2);
}
*/

/*
	This is the timelike electromagnetic form factor established in 
	arXiv:0910.5589. 
*/

double mrh = 0.77, mome = 0.77, mrhop = 1.25, mrhopp = 1.45, momep = 1.25, momepp = 1.45, gammarho= 0.15, gammaome = 0.0085, gammarhop=0.3, gammaomep= 0.3, gammarhopp = 0.5, gammaomepp = 0.5;

double f1ra = 0.6165340033101271; double f1rb = 0.22320420111672623; double f1rc = -0.33973820442685326;
double f1wa = 1.0117544786579074; double f1wb = -0.8816565944110686; double f1wc = 0.3699021157531611;

complex<double> F1r(double q2){
 	return f1ra*mrh*mrh/(mrh*mrh-complex<double>(q2,mrh*gammarho))+f1rb*mrhop*mrhop/(mrhop*mrhop-complex<double>(q2,mrhop*gammarhop))\
 	+f1rc*mrhopp*mrhopp/(mrhopp*mrhopp-complex<double>(q2,mrhopp*gammarhopp));
 } 

complex<double> F1w(double q2){
 	return f1wa*mome*mome/(mome*mome-complex<double>(q2,mome*gammaome))+f1wb*momep*momep/(momep*momep-complex<double>(q2,momep*gammaomep))\
 	+f1wc*momepp*momepp/(momepp*momepp-complex<double>(q2,momepp*gammaomepp));
 }

complex<double> F_1_proton(double q2){
	return F1r(q2)+F1w(q2);
}

complex<double> F_1_proton_baryonic(double q2, double kappa, double alD){
	return 0.5*(gu(alD, kappa)-gd(alD, kappa))*F1r(q2) + 1.5*(gu(alD, kappa)+gd(alD, kappa))*F1w(q2);
}

double Proton_Brem_Distribution::d2N_proton_brem_to_V(double z, double pt2){
	if(model=="proton_brem_baryonic"){
		return pow(std::abs(F_1_proton_baryonic(MA*MA,kappa,alpha_D)),2)*sigmapp(2*mp*(Beam_Energy-sqrt(MA*MA+pt2+z*z*(pow(Beam_Energy,2)-mp*mp))))/sigmapp(2*mp*Beam_Energy)*wpp(z,pt2, MA)/4.0/pi;
	}
	else 
		return pow(std::abs(F_1_proton(MA*MA)),2)*sigmapp(2*mp*(Beam_Energy-sqrt(MA*MA+pt2+z*z*(pow(Beam_Energy,2)-mp*mp))))/sigmapp(2*mp*Beam_Energy)*pow(kappa,2)*alphaEM*wpp(z,pt2,MA);
}

void Proton_Brem_Distribution::set_fit_parameters(production_channel &par){
	par.query_dist_param("rD",rD);
	// par.query_dist_param("mD2",mD2);
	par.query_dist_param("H",Hpp);
	par.query_dist_param("M",Mpp);
	par.query_dist_param("eta1",eta1pp);
	par.query_dist_param("eta2",eta2pp);
	par.query_dist_param("R1pp",R1pp);
	par.query_dist_param("R2pp",R2pp);
	par.query_dist_param("Ppp",Ppp);
	sppM = pow(2*mp+Mpp,2);
	calc_V_prod_rate();
	max_prod = d2N_proton_brem_to_V(ZMIN,PTMIN);
}

void Proton_Brem_Distribution::calc_V_prod_rate(){
	std::function<double(double, double)> func = std::bind(&Proton_Brem_Distribution::d2N_proton_brem_to_V,this,_1,_2);
	if(ZMIN>ZMAX){
		vprodrate=0;
	}
	else
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
