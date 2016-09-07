#include "sanfordwang.h"
#include "Random.h"
#include <math.h>
#include <iostream>
#include "constants.h"

using std::cout;
using std::endl;

//const double MPI = mpi0;
//const double META = meta;
const double thetamin=0, thetamax=pi/2.0;
//const double sigswpip=158.93469618126042;
//const double sigswpim=144.63;
//piplus values

const int burn_max=20000;

void sanfordwang::set_fit_parameters(production_channel &san){
	if(!(san.query_dist_param()))
		return;
	bool cval = false;
	bool dval = false;
	bool eval = false;
	if(san.query_dist_param(pB_key,pB)){
		cval=true;
		dval=true;
		eval=true;
	}
	if(san.query_dist_param(mommax_key,ppimax)){
		cval=true;
		dval=true;
		eval=true;
	}
	for(int i=1;i<10;i++){
		if(san.query_dist_param(ckey[i],c[i]))
			cval=true;
		if(san.query_dist_param(dkey[i],d[i]))
			dval=true;
		if(san.query_dist_param(ekey[i],e[i]))
			eval=true;
	}
	if(cval||dval){
		fpi0max=0;
		for(int i=0; i<burn_max;i++){
			pi0burnin();	
		}
	}
	if(eval){
		sigKmax=0;
		for(int i=0; i<burn_max; i++){
			K0burnin();
		}
	}
}

void sanfordwang::report(){
	for(int i=1; i<10; i++){
		cout << ckey[i] << " = " << c[i] << endl;
	}
	for(int i=1; i<10; i++){
		cout << dkey[i] << " = " << d[i] << endl;
	}
	for(int i=1; i<10; i++){
		cout << ekey[i] << " = " << e[i] << endl;
	}
	cout << "fpi0max = " << fpi0max << endl;
	cout << "sigKmax = " << sigKmax << endl;
	cout << "pB = " << pB << endl;
	cout << "ppimax = " << ppimax << endl;
}

// Differential cross section for Pi+ production 
// d^2 sigma / dOmega dp
double sanfordwang::swpip (const double p, const double theta) {
	return(c[1]*pow(p,c[2])*(1-p/(pB-c[9]))*exp(-c[3]*pow(p,c[4])/pow(pB,c[5])-c[6]*theta*(p-c[7]*pB*pow(cos(theta),c[8]))));
}

// Differential cross section for Pi- production 
// d^2 sigma / dOmega dp
double sanfordwang::swpim (const double p, const double theta) {
	return(d[1]*pow(p,d[2])*(1-p/(pB-d[9]))*exp(-d[3]*pow(p,d[4])/pow(pB,d[5])-d[6]*theta*(p-d[7]*pB*pow(cos(theta),d[8]))));
}

// Differential cross section for K0 production 
// d^2 sigma / dOmega dp
double sanfordwang::swK (const double p, const double theta) {
	return(e[1]*pow(p,e[2])*(1-p/(pB-e[9]))*exp(-e[3]*pow(p,e[4])/pow(pB,e[5])-e[6]*theta*(p-e[7]*pB*pow(cos(theta),e[8]))));
}

void sanfordwang::pi0burnin(){
	double hold, theta, pmom;
	while(true){
		theta = Random::Flat(thetamin,thetamax);
		pmom = Random::Flat(0, ppimax);
		if((hold=sin(theta)*(swpim(pmom, theta)+swpip(pmom, theta)))>(fpi0max*Random::Flat(0,1))){
			if(hold > fpi0max){
				fpi0max = hold;
			}
			break;
		}
	}
	
}

void sanfordwang::pi0sample(double &pmom, double &theta, double &pphi){
	double hold;
	while(true){
		theta = Random::Flat(thetamin,thetamax);
		pmom = Random::Flat(0, ppimax);
		if((hold=sin(theta)*(swpim(pmom, theta)+swpip(pmom, theta)))>(fpi0max*Random::Flat(0,1))){
			pphi=Random::Flat(0,2*pi);
			break;
		}
	}
}

void sanfordwang::K0burnin(){
	double hold, theta, pmom;
	while(true){
		theta = Random::Flat(thetamin,thetamax);
		pmom = Random::Flat(0, ppimax);
		if((hold=sin(theta)*(swK(pmom, theta)))>(sigKmax*Random::Flat(0,1))){
			if(hold > sigKmax)
				sigKmax = hold;
			break;
		}
	}
}

void sanfordwang::K0sample(double &pmom, double &theta, double &pphi){
	double hold;
	while(true){
		theta = Random::Flat(thetamin,thetamax);
		pmom = Random::Flat(0, ppimax);
		if((hold=sin(theta)*(swK(pmom, theta)))>(sigKmax*Random::Flat(0,1))){
			pphi=Random::Flat(0,2*pi);
			break;
		}
	}
}

void sanfordwang::sample_momentum(double &pmom, double &theta, double &phi){
	if(production_choice=="pi0_sanfordwang")
		pi0sample(pmom, theta, phi);
	else if(production_choice=="k0_sanfordwang")
		K0sample(pmom, theta, phi);
	else
		std::cerr << "Unknown production distribution" << endl;
}

void sanfordwang::sample_particle(Particle &part){
	double mom, theta, phi;
	sample_momentum(mom, theta, phi);
	part.ThreeMomentumPolar(mom,theta,phi);
}

/*
void sanfordwang::etaGen(Particle &eta){
    eta.m = META;
    double theta, peta;
    while(true){
        theta = Random::Flat(thetamin,thetamax);
        peta = Random::Flat(0, ppimax);
        if((sin(theta)*swK(peta, theta))>(sigKmax*Random::Flat(0,1))){
            double pphi=Random::Flat(0,2*pi);
            eta.ThreeMomentum(peta*cos(pphi)*sin(theta), peta*sin(pphi)*sin(theta), peta*cos(theta));
            break;
        }
    }
}*/
