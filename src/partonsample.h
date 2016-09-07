#ifndef PARTONSAMPLE_H
#define PARTONSAMPLE_H

#include <iostream>
#include <fstream>
#include <string>
#include <memory>
#include "Integrator.h"
#include "Distribution.h"

class parton_sample : public Distribution{

    public:
        parton_sample(std::string &pfilenmame, std::string &nfilename, double, double);
		void sample_momentum(double&, double&, double&);
		void sample_particle(Particle &);
        double production_neutron_cross_section(){return cross_section_n;}
		double production_proton_cross_section(){return cross_section_p;}
		double production_cross_section(){return cross_section_p*proton_number+cross_section_n*neutron_number;}
        double gen_V_mom();
    private:
		void parse_V_dist(std::string &pfilenmame, std::shared_ptr<Linear_Interpolation> &V_dist, double &cross_section, double &pmax, double &minmom, double &maxmom);
		std::shared_ptr<Linear_Interpolation> V_n_dist;
		std::shared_ptr<Linear_Interpolation> V_p_dist;
        double pmax_n, pmax_p, mom_min_n, mom_min_p, mom_max_n, mom_max_p, cross_section_n, cross_section_p;
		double neutron_number, proton_number;
};

#endif
