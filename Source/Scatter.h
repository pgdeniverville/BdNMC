#ifndef GUARD_scatter_h
#define GUARD_scatter_h
#include <iostream>
#include "detector.h"
#include "Particle.h"
#include "Integrator.h"
#include "Random.h"
#include <memory>
#include <vector>
#include <string>
/*
 *
 * The scatter class is designed to handle all of the scattering for a single detector in a given
 * scattering channel. If, in the future, I want to handle multiple detectors at a time, I will
 * need to change this to a vector of pMaxes, each with some way of associating with a detector
 * of interest.
 */
class Scatter{
	public:
		Scatter(){pMax = 0;}
		virtual ~Scatter(){};
		virtual bool probscatter(std::shared_ptr<detector>& det, Particle &DM, Particle &target) = 0;
		 virtual bool probscatter(std::shared_ptr<detector>& det, Particle &DM) = 0;
		void set_Model_Parameters(double MDM, double MV, double alphaprime, double kappa){mdm=MDM; MDP=MV; alD=alphaprime; kap=kappa; set_pMax(0);}
		double get_pMax(){return pMax;}
		void set_pMax(double pm){pMax=pm;}
		void Generate_Position(std::shared_ptr<detector>& det, Particle &DM, Particle &scat){
			DM.Generate_Position(Random::Flat(det->cross_point[0],det->cross_point[1]));
            scat.Set_Position(DM.end_coords[0],DM.end_coords[1],DM.end_coords[2]);
            scat.Set_Time(DM.end_coords[3]);
            scat.EVENT_SET = true;
		}
		//void set_scattering_energy(double emin, double emax){Escatmin=emin; Escatemax=emax;}
		//double get_MDP(){return kap;}
	protected:	
		double pMax, MDP, mdm, alD, kap;
		double Escatmax, Escatmin;
};

/*
 * Nucleon_scatter also handles the generation of interpolation functions for the nucleon scattering
 * cross section. 
 * Edmres is the energy resolution of the linear interpolation in incoming dark matter energy. 
 * Must be tuned!
 */ 

class Nucleon_Scatter: public Scatter{
	public:
		Nucleon_Scatter(double Emin, double Emax, double Eres, double MDM, double MV, double alphaprime, double kappa, double NEmin, double NEmax);
		~Nucleon_Scatter(){}
		bool probscatter(std::shared_ptr<detector>& det, Particle &DM, Particle &nucleon);
		bool probscatter(std::shared_ptr<detector>& det, Particle &DM);
		void set_Model_Parameters(double MDM, double MV, double alphaprime, double kappa){
			mdm=MDM; MDP=MV; alD=alphaprime; kap=kappa; set_pMax(0);
			generate_cross_sections();
		}
	private:
		double scatmax(double);
		double scatmin(double DM_Energy, double DM_Mass, double Nucleon_Mass);
		double Edmmin, Edmmax, Edmres;
		std::unique_ptr<Linear_Interpolation> proton_cross_maxima;
		std::unique_ptr<Linear_Interpolation> neutron_cross_maxima;
		std::unique_ptr<Linear_Interpolation> proton_cross;
		std::unique_ptr<Linear_Interpolation> neutron_cross;
		void scatterevent(Particle &DM, Particle &nucleon, std::function<double(double)>, Linear_Interpolation&);
		void generate_cross_sections();
};

class Nucleon_Scatter_Baryonic: public Scatter{
	public:
		Nucleon_Scatter_Baryonic(double Emin, double Emax, double Eres, double MDM, double MV, double alphaprime, double kappa, double NEmin, double NEmax);
		~Nucleon_Scatter_Baryonic(){}
		bool probscatter(std::shared_ptr<detector>& det, Particle &DM, Particle &nucleon);
		bool probscatter(std::shared_ptr<detector>& det, Particle &DM);
		void set_Model_Parameters(double MDM, double MV, double alphaprime, double kappa){
			mdm=MDM; MDP=MV; alD=alphaprime; kap=kappa; set_pMax(0);
			generate_cross_sections();
		}
	private:
		double scatmax(double);
		double scatmin(double DM_Energy, double DM_Mass, double Nucleon_Mass);
		double Edmmin, Edmmax, Edmres;
		std::unique_ptr<Linear_Interpolation> proton_cross;
		std::unique_ptr<Linear_Interpolation> neutron_cross;
		std::unique_ptr<Linear_Interpolation> proton_cross_maxima;
		std::unique_ptr<Linear_Interpolation> neutron_cross_maxima;
		void scatterevent(Particle &DM, Particle &nucleon, std::function<double(double)>, Linear_Interpolation&);
		void generate_cross_sections();
};

class Electron_Scatter: public Scatter{
	public:
		Electron_Scatter(double MDM, double MV, double alphaprime, double kappa, double eEmin, double eEmax);
		~Electron_Scatter(){}
		bool probscatter(std::shared_ptr<detector>& det, Particle &DM, Particle &electron);
		bool probscatter(std::shared_ptr<detector>& det, Particle &DM);
		void get_pMax(){}
	private:
		void scatterevent(Particle &DM, Particle &electron);
		double scatmax(double, double);
		double scatmin(double, double);
};

#endif
