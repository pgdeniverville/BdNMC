#ifndef GUARD_scatter_h
#define GUARD_scatter_h
#include <iostream>
#include "detector.h"
#include "Particle.h"
#include "Integrator.h"
#include "Random.h"
#include "DMgenerator.h"
#include <memory>
#include <vector>
#include <list>
#include <string>

void Prepare_Cross_Section(std::function<double(double,double)> f_init, std::function<double(double)> ER_Min, std::function<double(double)> ER_Max, std::shared_ptr<Linear_Interpolation>& cross, std::shared_ptr<Linear_Interpolation>& maxima, double E_in_min, double E_in_max, double E_in_res);

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
        void add_decay(const std::string parent_name, std::shared_ptr<DMGenerator> decay_gen){
            parent_name_vec.push_back(parent_name);
            decay_gen_vec.push_back(decay_gen);
            decay=true;
            }
		void set_Model_Parameters(double MDM, double MV, double alphaprime, double kappa){
			mdm=MDM; MDP=MV; alD=alphaprime; kap=kappa; set_pMax(0);}
		virtual double get_pMax(){return pMax;}
		void set_pMax(double pm){pMax=pm;}
		void set_angle_limits(double min, double max){min_angle=min;max_angle=max;}
		void set_energy_limits(double min, double max){Escatmin=min; Escatmax=max;}
        //void report(std::ostreami& out) const{out << pMax << std::endl;}
        //void set_scattering_energy(double emin, double emax){Escatmin=emin; Escatemax=emax;}
		//double get_MDP(){return kap;}
		//This tells the code to rotate the end state particle with 
		//bool set_end_state(bool t){_End_state_with_DM_parallel_to_z=t;}
	protected:	
        //If true, then one of the recoil particles can be expected to decay.
        bool decay;
		double pMax, MDP, mdm, alD, kap;
		double Escatmax, Escatmin,min_angle,max_angle;
		//bool _End_state_with_DM_parallel_to_z=false;
        std::vector<std::string> parent_name_vec;
        std::vector<std::shared_ptr<DMGenerator> > decay_gen_vec;
};

//This should replace all the elastic scattering channels.
//I would also like to generalize all inelastic scattering, but
//that seems much more complicated.
class Elastic_Scatter : public Scatter{
    public:
        Elastic_Scatter();
        ~Elastic_Scatter(){};
        void add_channel(Particle& end_state, std::function<double(double)>& cross_total, std::function<double(double)>& cross_max, std::function<double(double, double)>& dsig,double num_density);
        bool probscatter(std::shared_ptr<detector>& det, std::list<Particle>& partlist, std::list<Particle>::iterator&);
        bool probscatter(std::shared_ptr<detector>& det, Particle &part, Particle &recoil);
        bool probscatter(std::shared_ptr<detector>& det, Particle &part);
        void set_Model_Parameters(){};
        double get_pMax(){return pMax;}
    private:
        double scatmax(double E1, double m1, double m2);
        double scatmin(double E1, double m1, double m2);
        void scatterevent(Particle &particle, Particle &recoil, std::function<double(double)> cross_section, std::function<double(double)> cross_maxima);
        
        //total value of the cross section for a given signal channel
        std::vector<std::function<double(double)> > cross_tot; 
        //maximum value of the cross section for a given energy.
        std::vector<std::function<double(double)> > cross_maxima;
        //differential cross section dsigma(E_DM_i,E_recoil)/dE_recoil
        std::vector<std::function<double(double, double)> > dsigma;
        //Recoil particle
        std::vector<Particle> recoil_particle;
        //number density associated with a given scattering channel. 
        std::vector<double> number_density;
        int chan_number;
};

//Here's the effort to generalize all inelastic scattering.
class Two_to_Two_Scatter : public Scatter{
    public:
        Two_to_Two_Scatter();
        ~Two_to_Two_Scatter(){};
        void add_channel(Particle& out_state, Particle& end_state, double targ_mass, std::function<double(double)>& cross_total, std::function<double(double)>& cross_max, std::function<double(double, double)>& dsig, double num_density, const std::string in_state);
        void Build_Channel(Particle& out_state, Particle& end_state, double in_mass, double targ_mass, std::function<double(double,double)> &dsig, double num_density, const std::string in_state, double Max_Energy, double EDM_RES);
        bool probscatter(std::shared_ptr<detector>& det, std::list<Particle>& partlist, std::list<Particle>::iterator& partit);
        bool probscatter(std::shared_ptr<detector>& det, Particle& part);
        bool probscatter(std::shared_ptr<detector>& det, Particle &part, Particle &recoil, Particle &recoil2);
        double scatmax(double E1, double m1, double m2, double m3, double m4);
        double scatmin(double E1, double m1, double m2, double m3, double m4);
        double get_pMax(){return pMax;}
    private:
        int chan_number;

        void scatterevent (Particle &part, double targ_mass, Particle &recoil, Particle &recoil2, std::function<double(double)> Xsec,std::function<double(double)> Xmax);
        //total value of the cross section for a given signal channel.
        std::vector<std::function<double(double)> > cross_tot; 
        //maximum value of the cross section for a given incoming energy.
        std::vector<std::function<double(double)> > cross_maxima;
        //differential cross section dsigma(E_1,E_4)/dE_4.
        std::vector<std::function<double(double, double)> > dsigma;
        //Recoil Target particle. This is particle 4.
        std::vector<Particle> recoil_particle;
        //Recoiling other particle. If one is invisible, it's this one! This is particle 3
        std::vector<Particle> out_particle;
        //number density associated with a given scattering channel. 
        std::vector<double> number_density;
        //name of incoming particle.
        std::vector<std::string> in_name;
        //mass of target particle for a given channel. This is particle 2
        std::vector<double> target_mass;
};

/*
 * Nucleon_scatter also handles the generation of interpolation functions for the nucleon scattering
 * cross section. 
 * Edmres is the energy resolution of the linear interpolation in incoming dark matter energy. 
 * Must be tuned!
 */ 

class Nucleon_Scatter: public Scatter{
	public:
		//Nucleon_Scatter(double Emin, double Emax, double Eres, double MDM, double MV, double alphaprime, double kappa, double NEmin, double NEmax, std::string model);
		Nucleon_Scatter(double Emin, double Emax, double Eres, double MDM, double MV, double alphaprime, double kappa, double NEmin, double NEmax, std::string model, const bool cohere, std::shared_ptr<detector> det);
		~Nucleon_Scatter(){}
		bool probscatter(std::shared_ptr<detector>& det, std::list<Particle> &partlist, std::list<Particle>::iterator& it);
		bool probscatter(std::shared_ptr<detector>& det, Particle &DM, Particle &nucleon);
		bool probscatter(std::shared_ptr<detector>& det, Particle &DM);
		void set_Model_Parameters(double MDM, double MV, double alphaprime, double kappa){
			mdm=MDM; MDP=MV; alD=alphaprime; kap=kappa; set_pMax(0);
			generate_cross_sections();
		}
		void set_Model_Parameters(double MDM, double MV, double alphaprime, double kappa, std::shared_ptr<detector>& det){
			mdm=MDM; MDP=MV; alD=alphaprime; kap=kappa; set_pMax(0);
            generate_coherent_cross_sections(det);
        }
        double get_pMax(){return pMax;}
	private:
        bool coherent;
		double scatmax(double);
		double scatmin(double DM_Energy, double DM_Mass, double Nucleon_Mass);
		double Edmmin, Edmmax, Edmres;
		Linear_Interpolation proton_cross_maxima;
		Linear_Interpolation neutron_cross_maxima;
		Linear_Interpolation proton_cross;
		Linear_Interpolation neutron_cross;
       
        //Function pointers to cross section functions!
        double (*dsig_prot)(double, double, double, double, double, double);
        double (*dsig_neut)(double, double, double, double, double, double);
        double (*dsig_cohp)(double, double, double, double, double, double, double, double);

        //Coherent Variables
        std::vector<Linear_Interpolation> atom_maxima;
        std::vector<Linear_Interpolation> atom_cross;
        
        void scatterevent(Particle &DM, Particle &nucleon, std::function<double(double)>, Linear_Interpolation&);
		void generate_cross_sections();
		void generate_coherent_cross_sections(std::shared_ptr<detector>& det);
         
        void cross_gen_handler(std::function<double(double)> fp, std::function<double(double)> fplim,
                std::vector<double> &cross_vec, std::vector<double> &cross_vec_maxima,
                double iter, double m1, const double m2); 
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
        double get_pMax(){return pMax;}
	private:
		void scatterevent(Particle &DM, Particle &nucleon, std::function<double(double)>, Linear_Interpolation&);
		void load_cross_sections(const std::string &filename);
		std::unique_ptr<Linear_Interpolation> scatter_dist;//Proton
		std::unique_ptr<Linear_Interpolation> scatter_dist_n;//Neutron
		std::string channel_name;
		double E_min, E_max;
};
/*Deprecated!
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
*/
class Electron_Scatter: public Scatter{
	public:
		Electron_Scatter(double MDM, double MV, double alphaprime, double kappa, double eEmin, double eEmax);
		~Electron_Scatter(){}
		bool probscatter(std::shared_ptr<detector>& det, std::list<Particle> &partlist, std::list<Particle>::iterator& it);
		bool probscatter(std::shared_ptr<detector>& det, Particle &DM, Particle &nucleon);
		bool probscatter(std::shared_ptr<detector>& det, Particle &DM);
        double get_pMax(){return pMax;}
	private:
		void scatterevent(Particle &DM, Particle &electron);
		double scatmax(double, double);
		double scatmin(double, double);
};


//final_state=0 -> Inelastic pi0 production
//final_state=1 -> Inelastic photon production
class Pion_Inelastic: public Scatter{
	public:
		Pion_Inelastic(double Emin, double Emax, double Eres, double MDM, double MV, double alphaprime, double kappa, double NEmin, double NEmax,int final_state=0);
		~Pion_Inelastic(){}
		bool probscatter(std::shared_ptr<detector>& det, std::list<Particle> &partlist, std::list<Particle>::iterator& it);
		bool probscatter(std::shared_ptr<detector>& det, Particle &DM,Particle &pion, Particle &Delta,Particle &nucleon);
		bool probscatter(std::shared_ptr<detector>& det, Particle &DM);
		void set_Model_Parameters(double MDM, double MV, double alphaprime, double kappa){
			mdm=MDM; MDP=MV; alD=alphaprime; kap=kappa; set_pMax(0);
			generate_cross_sections();
		}
        double get_pMax(){return pMax;}
	private:
		double Ermax(double DM_Energy, double DM_Mass, double Nucleon_Mass);
		double Ermin(double DM_Energy, double DM_Mass, double Nucleon_Mass);
		void prep_ab(double &, double &, double, double, double);
		double Er_to_theta(double DM_Energy, double Delta_Energy, double DM_Mass, double Nucleon_Mass);
        double final_mass,final_branch;
        std::string final_name;		
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
