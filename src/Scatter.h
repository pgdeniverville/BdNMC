#ifndef GUARD_scatter_h
#define GUARD_scatter_h
#include <iostream>
#include "detector.h"
#include "Particle.h"
#include "Integrator.h"
#include "Random.h"
#include <memory>
#include <vector>
#include <list>
#include <string>
/*
 * The scatter class is designed to handle all of the scattering for a single detector in a given
 * scattering channel. If, in the future, I want to handle multiple detectors at a time, I will
 * need to change this to a vector of pMaxes, each with some way of associating with a detector
 * of interest.
 */
class Scatter{
	public:
		Scatter(){pMax = 0;}
		virtual ~Scatter(){};
		virtual bool probscatter(std::shared_ptr<detector>& det, std::list<Particle>& partlist, std::list<Particle>::iterator&) = 0;
		virtual bool probscatter(std::shared_ptr<detector>& det, Particle &DM) = 0;
		void set_Model_Parameters(double MDM, double MV, double alphaprime, double kappa){
			mdm=MDM; MDP=MV; alD=alphaprime; kap=kappa; set_pMax(0);}
		double get_pMax(){return pMax;}
		void set_pMax(double pm){pMax=pm;}
		void Generate_Position(std::shared_ptr<detector>& det, Particle &DM, Particle &scat){
			DM.Generate_Position(Random::Flat(det->cross_point[0],det->cross_point[1]));
            scat.Set_Creation_Time(DM.end_coords[3]);
            scat.Set_Origin(DM.end_coords[0],DM.end_coords[1],DM.end_coords[2]);
			scat.EVENT_SET = true;
		}
		void set_angle_limits(double min, double max){min_angle=min;max_angle=max;}
		//void set_scattering_energy(double emin, double emax){Escatmin=emin; Escatemax=emax;}
		//double get_MDP(){return kap;}
		//This tells the code to rotate the end state particle with 
		//bool set_end_state(bool t){_End_state_with_DM_parallel_to_z=t;}
	protected:	
		double pMax, MDP, mdm, alD, kap;
		double Escatmax, Escatmin,min_angle,max_angle;
		//bool _End_state_with_DM_parallel_to_z=false;
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
		bool probscatter(std::shared_ptr<detector>& det, std::list<Particle> &partlist, std::list<Particle>::iterator& it);
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

//This is good for inelastic
class Inelastic_Nucleon_Scatter: public Scatter{
	public:
		Inelastic_Nucleon_Scatter(double MDM, double mv, double alphaprime, double kappa, const std::string &chan_name, const std::string &filename);
		~Inelastic_Nucleon_Scatter(){}
		bool probscatter(std::shared_ptr<detector>& det, std::list<Particle> &partlist, std::list<Particle>::iterator& it);
		bool probscatter(std::shared_ptr<detector>& det, Particle &DM);
		void set_Model_Parameters(double MDM, double MV, double alphaprime, double kappa){
			mdm=MDM; MDP=MV; alD=alphaprime; kap=kappa; set_pMax(0);
		}
	private:
		void scatterevent(Particle &DM, Particle &nucleon, std::function<double(double)>, Linear_Interpolation&);
		void load_cross_sections(const std::string &filename);
		std::unique_ptr<Linear_Interpolation> scatter_dist;//Proton
		std::unique_ptr<Linear_Interpolation> scatter_dist_n;//Neutron
		std::string channel_name;
		double E_min, E_max;
};

class Nucleon_Scatter_Baryonic: public Scatter{
	public:
		Nucleon_Scatter_Baryonic(double Emin, double Emax, double Eres, double MDM, double MV, double alphaprime, double kappa, double NEmin, double NEmax);
		~Nucleon_Scatter_Baryonic(){}
		bool probscatter(std::shared_ptr<detector>& det, std::list<Particle> &partlist, std::list<Particle>::iterator& it);
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
		bool probscatter(std::shared_ptr<detector>& det, std::list<Particle> &partlist, std::list<Particle>::iterator& it);
		bool probscatter(std::shared_ptr<detector>& det, Particle &DM, Particle &nucleon);
		bool probscatter(std::shared_ptr<detector>& det, Particle &DM);
		void get_pMax(){}
	private:
		void scatterevent(Particle &DM, Particle &electron);
		double scatmax(double, double);
		double scatmin(double, double);
};

class Pion_Inelastic: public Scatter{
	public:
		Pion_Inelastic(double Emin, double Emax, double Eres, double MDM, double MV, double alphaprime, double kappa, double NEmin, double NEmax);
		~Pion_Inelastic(){}
		bool probscatter(std::shared_ptr<detector>& det, std::list<Particle> &partlist, std::list<Particle>::iterator& it);
		bool probscatter(std::shared_ptr<detector>& det, Particle &DM, Particle &nucleon);
		bool probscatter(std::shared_ptr<detector>& det, Particle &DM);
		void set_Model_Parameters(double MDM, double MV, double alphaprime, double kappa){
			mdm=MDM; MDP=MV; alD=alphaprime; kap=kappa; set_pMax(0);
			generate_cross_sections();
		}
	private:
		double Ermax(double DM_Energy, double DM_Mass, double Nucleon_Mass);
		double Ermin(double DM_Energy, double DM_Mass, double Nucleon_Mass);
		void prep_ab(double &, double &, double, double, double);
		double Er_to_theta(double DM_Energy, double Delta_Energy, double DM_Mass, double Nucleon_Mass);
		
		double dsigma_dER_N(double En, double ER, double mx, double mA, double alphaprime, double kappa, double mN);
		double GM(double q2);
		double Edmmin, Edmmax, Edmres;
		std::unique_ptr<Linear_Interpolation> proton_cross_maxima;
		std::unique_ptr<Linear_Interpolation> neutron_cross_maxima;
		std::unique_ptr<Linear_Interpolation> proton_cross;
		std::unique_ptr<Linear_Interpolation> neutron_cross;
		std::shared_ptr<Linear_Interpolation> form_factor;
		void scatterevent(Particle &DM, Particle &nucleon, std::function<double(double)>, Linear_Interpolation&, double nuclmass);
		void generate_cross_sections();
		void load_form_factor(const std::string &filename, std::shared_ptr<Linear_Interpolation>& form);
		double MAX_Q2;
		bool MAX_Q2_WARNING;
};
#endif
