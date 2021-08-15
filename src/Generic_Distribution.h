#ifndef GUARD_generic_dist_h
#define GUARD_generic_dist_h

#include "Distribution.h"

#include <functional>

class Generic_Distribution : public Distribution{

	public:
		Generic_Distribuiton(double _Beam_Energy, std::string distribution);
		void sample_momentum(Particle &);
		void sample_momentum(double&, double&, double&);
	private:
		double Beam_Energy, Part_Mass;
		double s, p_max, mom_max;
		std::string dist_name;
		std::function<double E, double theta> dist_func;
		double JPsi(double p, double theta);
		double JPsi_dist_x(double x, double s, double gamma_CM);
		double JPsi_dist_pT(double p_Trans,double s);
};

#endif