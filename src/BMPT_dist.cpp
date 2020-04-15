#include "BMPT_dist.h"

#include "constants.h"
#include "Random.h"
#include <math.h>
#include <iostream>

const double mp = MASS_PROTON;
const double mn = MASS_NEUTRON;

const int pi0_choice=0;
const int k0_choice=1;

//BMPT fit parameters

const double A = 62.3;
const double B = 1.57;
const double alpha = 3.45;
const double beta = 0.517;
const double a = 6.10;
const double gammafit = 0.153;
const double deltafit = 0.478;
const double r0 = 1.05;
const double r1 = 2.65;
const double MXBAR = 1.88;//GeV

namespace Kaon{
    const double A = 7.74;
    const double B = 0;
    const double alpha = 2.45;
    const double beta = 0.444;
    const double a = 5.04;
    const double gammafit = 0.121;
    const double deltafit = 2*gammafit;
    const double r0 = 1.15;
    const double r1 = -3.17;
}

//Should be moved into constants.h?
const double Be_Mass_Number = 8; //This distribution is calibrated for Beryllium, but can be scaled to other materials.

using std::cout;
using std::endl;

BMPT::BMPT(double beamE, int mass_number, int choice){
    Beam_Energy = beamE;
    Mass_Number = mass_number;
    particle_choice=choice;

    if(choice==pi0_choice){
        Meson_Mass =  MASS_PION_CHARGED;//not really sure about this, but I think it's the safest way to do it.
    }
    else if(choice==k0_choice){
        Meson_Mass = MASS_KAON_CHARGED;
    }
    else{
        std::cerr << "BMPT supplied a non-existent particle_choice number " << choice << endl;
        throw -1;
    }
    sBMPT = 2*(mp*mp+mp*Beam_Energy);
    Beta_CM = sqrt(pow(Beam_Energy,2)-mp*mp)/(Beam_Energy+mp);
    prob_max = 0.0;
    p_max = Beam_Energy-2*mp-Meson_Mass;//Do not have a good reason for this, but the parameterization seems to work pretty well no matter the p value. Real problem is that thetamax is p dependent.
    theta_max = pi/2.0;//Kinematics often make this smaller, hopefully I am sampling correctly.
    double p, t, phi;
    //Burn-In
    for(int i=0; i<10000; i++){
        sample_momentum(p, t, phi);
    }
}

double BMPT::xR(double p, double theta){
    return (sBMPT - (2*mp*mp + pow(Meson_Mass,2) + 2*Beam_Energy*mp -\
                2*Beam_Energy*sqrt(p*p+pow(Meson_Mass,2))+2*sqrt(pow(Beam_Energy,2)-mp*mp)*p*cos(theta)\
                - 2*sqrt(p*p+pow(Meson_Mass,2))*mp)+pow(Meson_Mass,2))/(sBMPT-pow(MXBAR,2)+ pow(Meson_Mass,2));
}

double gammaboost(double beta_boost){
    return 1.0/sqrt(1-pow(beta_boost,2));
}

double BMPT::xF(double p, double theta){
    return gammaboost(Beta_CM)*(p*cos(theta)-Beta_CM*sqrt(p*p+pow(Meson_Mass,2)))/sqrt(sBMPT);
}

double aprime(double x){
    return a/pow(x,gammafit);
}

double bprime(double x){
    return a*a/(2*pow(x,deltafit));
}

double aprime_Kaon(double x){
    return Kaon::a/pow(x,Kaon::gammafit);
}

double bprime_Kaon(double x){
    return pow(Kaon::a,2)/(2*pow(x,Kaon::deltafit));
}

double alphafunc(double _xF, double pT){
    return (0.74-0.55*_xF+0.26*pow(_xF,2))*(0.98+0.21*pT);
}

double BMPT::Invariant_Cross_Section(double p, double theta){
    double xRhold = xR(p,theta);
    double pT = p*sin(theta);
    return A*pow(1.0-xRhold,alpha)*(1.0+B*xRhold)*pow(xRhold,-1.0*beta)*(1+aprime(xRhold)*pT+bprime(xRhold)*pT*pT)*exp(-aprime(xRhold)*pT)*pow(Mass_Number/Be_Mass_Number,alphafunc(xF(p, theta), pT));
}

double BMPT::Invariant_Cross_Section_Kaon(double p, double theta){
    double xRhold = xR(p,theta);
    double pT = p*sin(theta);
    return A*pow(1.0-xRhold,Kaon::alpha)*(1.0+Kaon::B*xRhold)*pow(xRhold,-1.0*Kaon::beta)*(1+aprime_Kaon(xRhold)*pT+bprime_Kaon(xRhold)*pT*pT)*exp(-aprime_Kaon(xRhold)*pT)*pow(Mass_Number/Be_Mass_Number,alphafunc(xF(p, theta), pT));
}

double BMPT::pratio(double p, double theta){
    return r0*pow(1+xR(p,theta),r1); 
}

double BMPT::pratio_Kaon(double p, double theta){
    return Kaon::r0*pow(1+xR(p,theta),Kaon::r1); 
}

double BMPT::Invariant_Cross_Section_pi_minus(double p, double theta){
    return Invariant_Cross_Section(p, theta)/r0/pow(1+xR(p,theta),r1);
}

double BMPT::Invariant_Cross_Section_K_minus(double p, double theta){
    return Invariant_Cross_Section_Kaon(p, theta)/Kaon::r0/pow(1+xR(p,theta),Kaon::r1);
}

void BMPT::sample_particle(Particle &part){
	double mom, theta, phi;
	sample_momentum(mom, theta, phi);
	part.ThreeMomentumPolar(mom, theta, phi);
}

void BMPT::sample_momentum(double &p, double &theta, double &phi){
    double prob;
    while(true){
        p = Random::Flat(0,1)*p_max;
        theta = Random::Flat(0,1)*theta_max;
        while(true){//I'm a little uncertain of this sampling scheme.
            if(xR(p,theta)<=1)
                break;
            theta = (theta*Random::Flat(0,1));
        }
        if(particle_choice==pi0_choice){
            prob = 2*pi*sin(theta)*pow(p,2)/sqrt(p*p+pow(Meson_Mass,2))*Invariant_Cross_Section(p, theta)*(1.0+1.0/pratio(p, theta));
        }
        else if(particle_choice==k0_choice){
            prob = 2*pi*sin(theta)*pow(p,2)/sqrt(p*p+pow(Meson_Mass,2))*Invariant_Cross_Section_Kaon(p, theta)*(1.0+3.0/pratio(p, theta))/2;
        }
        else{
            std::cerr << "sample_momentum supplied an invalid particle_choice " << particle_choice << endl;
            throw -1;
        }
        double u = Random::Flat(0,1);
        if(prob>(u*prob_max)){
            if(prob>prob_max)
                prob_max = prob;

            phi = Random::Flat(0,1)*2*pi;
            break;
        }
    }
}
