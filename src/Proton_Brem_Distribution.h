#ifndef GUARD_Proton_Brem_Dist_h
#define GUARD_Proton_Brem_Dist_h

#include "Parameter.h"
#include "Distribution.h"
#include <string>

class Proton_Brem_Distribution : public Distribution{
		
	public:
		Proton_Brem_Distribution(double Beam_E, double epsilon, double mA, double ptmax, double zmax, double zmin, double alphaD, std::string &mode, double ptmin=0);
		double V_prod_rate(){return vprodrate;}
		double d2N_proton_brem_to_V(double z, double pt2);
		void sample_particle(Particle &);
		void sample_momentum(double &, double &, double &);
		void set_fit_parameters(production_channel &); 
	private:
		void calc_V_prod_rate();
		double sigmapp(double s);
		std::string model;//proton_brem_baryonic triggers leptophobic behavior
		// double F_1_proton(double q2);
			
		double Beam_Energy, kappa, alpha_D, PTMIN, PTMAX, ZMAX, ZMIN, MA;
		double vprodrate, max_prod;
			
		//Proton Form Factor parameters
		double rD = 0.8/0.197; 
		//double mD2 = 12.0/rD/rD;
		//Proton Cross Section fit parameters
		double Hpp = 0.2704; 
		double Mpp=2.2127; 
		double sppM; 
		double eta1pp=0.451; 
		double eta2pp=0.549; 
		double R1pp=12.98; 
		double R2pp=7.38; 
		double Ppp = 34.49;
};

#endif
