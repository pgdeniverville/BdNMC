#ifndef GUARD_BMPT_dist_h
#define GUARD_BMPT_dist_h

#include "Distribution.h"
//Distribution from arxiv.org/abs/hep-ph/0101163v3

class BMPT : public Distribution{

    public:
        //void pionGen(Particle &);
        //void pionGen(double &, double &, double &, double &);
        BMPT(double beam_Energy, int Mass_Number);
        void sample_momentum(double&,double&,double&);
		void sample_particle(Particle &);
		double Invariant_Cross_Section(double p, double theta);//E \times d^3\sigma/dp^3
        double Invariant_Cross_Section_pi_minus(double p, double theta);//E \times d^3\sigma/dp^3
    private:
        double Beam_Energy, Mass_Number, Meson_Mass;
        double sBMPT, Beta_CM;
        double theta_max, p_max;
        double xR(double p, double theta);
        double xF(double p, double theta);
        double pratio(double p, double theta);
        double prob_max;
};

#endif
