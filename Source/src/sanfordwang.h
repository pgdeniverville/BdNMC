#ifndef GUARD_sanford_h
#define GUARD_sanford_h

// class sanfordwang
// Taken from 
// http://arxiv.org/abs/0806.1449

#include "Parameter.h"//This is annoying. I should get rid of it if I refactor.
#include "Distribution.h"
#include <string>

class sanfordwang: public Distribution{
public:
    sanfordwang(const std::string p_choice){production_choice=p_choice;}
    //Many of these could be made private
	void sample_momentum(double &, double &, double &);
	void sample_particle(Particle &);
	void pi0sample(double &, double &, double &);
    void pi0burnin();
	void K0sample (double &, double &, double &);
	void K0burnin();	
	void set_fit_parameters(production_channel &);
	double swpip (const double p, const double theta);
	double swpim (const double p, const double theta);
	double swK(const double p, const double theta);
	void report();
private:
	//These aren't exact, they're a little higher than you'd think	
	double fpi0max=54.25;
	double sigKmax=0.6506;
	std::string production_choice;
	double pB = 8.89;
	std::string pB_key="Beam_Energy";
	double ppimax = 7.8900;
	std::string mommax_key="maximum_pi_momentum";
	//These are the default values for the sanford-wang fits. These may be changed after initialization.
	//pi+ fit
	double c [10] = {0,220.7000,1.080000,1.000000,1.978000,1.320000,5.572000,0.086800,9.686000,1.000000};
	const std::string ckey[10] = {"","c1","c2","c3","c4","c5","c6","c7","c8","c9"};
	//pi- fit
	double d [10] = {0,213.7000,0.937900,5.454000,1.210000,1.284000,4.781000,0.073380,8.329000,1.000000};
	const std::string dkey[10] = {"","d1","d2","d3","d4","d5","d6","d7","d8","d9"};
	//K0 fit
	double e [10] = {0,15.130,1.975,4.084,0.928,0.731,4.362,0.048,13.300,1.278};

	const std::string ekey[10] = {"","e1","e2","e3","e4","e5","e6","e7","e8","e9"};

};


#endif
