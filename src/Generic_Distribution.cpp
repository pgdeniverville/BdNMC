#include "Generic_Distribution.h"

#include "constants.h"
#include "Random.h"

#include <cmath>

using std::string; using std::bind;


//From Saeid.
const double JPsi_norm = 0002421037;
const double JPsi_n = 0.54;

using namespace std::placeholders;

//Beam_Energy is given in GeV, assumed to be a proton impacting on proton.
//Could add a switch later to recalculate for an electron beam.
Generic_Distribution::Generic_Distribution(double _Beam_Energy, string distribution){
	Beam_Energy = _Beam_Energy;
	s = 2*MASS_PROTON*(Beam_Energy+MASS_PROTON);
	dist_name = distribution;
	if(distribution == "J/Psi_Distribution"){
		dist_func = bind(this,Generic_Distribution::JPsi,_1,_2);
		Part_Mass=MASS_JPSI;
		mom_max = 0.5*sqrt(s)*(1-pow(MASS_JPSI,2)/s);
	}
}


double Generic_Distribution::JPsi_dist_x(double x, double s, double gamma_CM){
    double c_JPsi = 2+(gamma_CM-5)/5;
    return (c_JPsi+1)/2*pow(1-abs(x),c_JPsi);
}

double Generic_Distribution::JPsi_dist_pT(double p_Trans,double s){
    double b = 2.66 + 3e-4*(s-1000);
    double Normalization = b**2.0/5.0;
    return (1.0/Normalization)*2.0*p_Trans*pow((1.0+(p_Trans/b)**2),-6);
}

//Saeid gave me this distribution in a python file.
double Generic_Distribution::JPsi(double p, double theta){
	double E = sqrt(p*p - pow(MASS_JPSI,2));
	double p_Trans = p*sin(theta);
	double gamma_CM = 0.5*sqrt(s)/MASS_PROTON;
	double beta_CM = sqrt(1-1.0/pow(gamma_CM,2));
	double Jacobian = -gamma_CM/mom_max*(beta_CM*p*cos(theta) - E);
	double x_F = -gamma_CM/mom_max*(beta_CM*E - p*cos(theta));
	//return sin(theta)*(1/Psi_norm)*cos(theta)*exp(-Psi_n*(p*sin(theta))**2);
	return Jacobian*JPsi_dist_x(x_F,s,gamma_CM)*JPsi_dist_pT(p_Trans,s);
}

void Generic_Distribution::sample_momentum(double &p, double &theta, double &phi){
    double prob;
    while(true){
        p = Random::Flat(0,1)*p_max;
        theta = Random::Flat(0,1)*theta_max;
        while(true){//I'm a little uncertain of this sampling scheme.
	        prob = dist_func()
	        double u = Random::Flat(0,1);
	        if(prob>(u*prob_max)){
	            if(prob>prob_max)
	                prob_max = prob;
	            phi = Random::Flat(0,1)*2*pi;
	            break;
	        }
    	}
	}
}

void Generic_Distribution::sample_particle(Particle &part){
	double mom, theta, phi;
	sample_momentum(mom, theta, phi);
	part.ThreeMomentumPolar(mom, theta, phi);
}
