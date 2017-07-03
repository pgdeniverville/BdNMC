#ifndef Parameter_h
#define Parameter_h

#include <ios>
#include <iostream>
#include <fstream>
#include <string>
#include <map>
#include <list>
#include "detector.h"
#include <memory>

class production_distribution{
	public:
		production_distribution();
		std::string dist_mod;
		double offset_coords[4] = {0,0,0,0};
		double get_offset(int i);
		std::string name(){return dist_mod;}
		std::string DIST_MOD(){return dist_mod;}
};

class production_channel{
	public:
		std::string prod_chan, prod_dist, particle_list_file;
		std::string sanfordwang_file;
		std::string parton_V_n_file, parton_V_p_file;
		//This map holds sanfordwang parameters.
		std::map<std::string, std::string> dist_param_map;
		double meson_per_pi0, PTMAX, PTMIN, ZMIN, ZMAX;
		bool particle_list_position;
		production_channel();
		bool query_dist_param(){return (dist_param_map.size()!=0);}
		bool query_dist_param(const std::string &, double&);
		bool par_list_pos(){return particle_list_position;} 
		double Meson_Per_Pi0(){return meson_per_pi0;}
        double ptmin(){return PTMIN;}
		double ptmax(){return PTMAX;}
		double zmin(){return ZMIN;}
		double zmax(){return ZMAX;}
		std::string Prod_Dist(){return prod_dist;}
		std::string Part_List_File(){return particle_list_file;}
		std::string Parton_V_Neutron_File(){return parton_V_n_file;}
		std::string Parton_V_Proton_File(){return parton_V_p_file;}
        std::string Production_Channel(){return prod_chan;}
		std::shared_ptr<std::list<production_distribution> > Get_Dist_Mods_List(){return dist_mods;}
		std::shared_ptr<std::list<production_distribution> > dist_mods;//Should probably be private
};

class Parameter{
    public:
        Parameter(std::ifstream&);
        ~Parameter(){}
       
        //Allows the user to query the map directly.
        //This will be used by future model implementations.
        bool Query_Map(const std::string, double&);

        int Integrity(){return integrity;}
        double MassDM(){return mass_dm;}
        double MassDP(){return mass_V;}
        double kap(){return epsilon;}
        double Epsilon(){return epsilon;}
        double alD(){return alpha_D;}
        std::string Model_Name(){return model_name;}


		double Beam_Energy(){return beam_energy;}
		double Max_DM_Energy(){return max_dm_energy;}
		double EDM_RES(){return edmres;}
		double Min_Scatter_Energy(){return min_scatter_energy;}
		double Max_Scatter_Energy(){return max_scatter_energy;}
		double Max_Angle(){return angle_upper_limit;}
		double Min_Angle(){return angle_lower_limit;}
		double Timing_Cut(){return timing_cut;}
		std::string Scatter_Dist_Filename(){return scatter_dist_filename;}

		double Target_E_Num(){return target_e_num;}
		double Target_N_Num(){return target_n_num;}
		double Target_P_Num(){return target_p_num;}
		double Target_N_Dens(){return target_n_density;}
		double Target_Length(){return target_length;}
		double Target_Mass(){return target_mass;}
	
		double P_Cross(){return p_cross;}

		double Efficiency(){return efficiency;}
		int Sample_Size(){return samplesize;}
		int Burn_In(){return burn_in;}
		int Burn_Timeout(){return burn_timeout;}
		int Max_Trials(){return max_trials;}
		int Repeat(){return repeat;}
        int Seed(){return seed;}
		double Protons_on_Target(){return POT;}
		double Pi0_per_POT(){return pi0ratio;}
		//Indicates if any custom sanfordwang parameters have been declared.
		std::string Signal_Channel(){return sig_channel;}
        bool Coherent(){return coherent;}
        std::string Output_File(){return output_file;}
		std::string Output_Mode(){return output_mode;}
		std::string Summary_File(){return summary_file;}
        std::shared_ptr<detector> Get_Detector(){return det;}
		std::shared_ptr<std::list<production_channel> > Get_Production_List(){return prodlist;}
		std::string Run_Name(){return run_name;}
    private:
        std::map<std::string, std::string> keymap;
            
        std::shared_ptr<detector> det;
        //std::string prod_channel;
		//std::string prod_dist;
		std::string sig_channel;
		//std::string particle_list_file;
		//std::string parton_V_neutron_file;
		//std::string parton_V_proton_file;
		std::string run_name;	
		std::shared_ptr<std::list<production_channel> >  prodlist; 
	

		double angle_upper_limit;
		double angle_lower_limit;
		double edmres;
		int burn_in;
		int burn_timeout;
		int repeat;
		int seed;
		int sample_gen;
		int max_trials;
		double meson_per_pi0;
        std::string output_file;
        std::string summary_file;
		std::string output_mode;
        bool coherent;

		double beam_energy;
		double max_scatter_energy;
		double min_scatter_energy;
		double max_dm_energy;
		double timing_cut;
		std::string scatter_dist_filename="";

		double target_length;
		double target_n_density;
		double target_mass;
		double target_e_num;
		double target_n_num;
		double target_p_num;

		double p_cross;	
	
        std::string model_name = "Hidden_Sector_DM";

		double efficiency;
		double mass_dm;
        double mass_V;
        double epsilon;
        double alpha_D;
        double POT;
        double pi0ratio;
		int samplesize;
        int integrity;
		//void parse_parameter_file(const std::string &, std::map<std::string, std::string> &parammap);
        void Set_Double(const std::string &, double &, std::map<std::string, std::string> &keymap);
        void Set_Integer(const std::string &, int &, std::map<std::string, std::string> &keymap);
        void Set_Bool(const std::string &, bool &, std::map<std::string, std::string> &keymap, const bool &def);
        void Set_Double(const std::string &, double &, std::map<std::string, std::string> &keymap, const double &def);
        void Set_Integer(const std::string &, int &, std::map<std::string, std::string> &keymap, const int &def);
		void Set_String(const std::string &, std::string &, std::map<std::string, std::string> &keymap, const std::string &def);
		void Set_Model_Parameters(std::map<std::string, std::string> &keymap);
        void Build_Detector(std::map<std::string, std::string> &keymap);
		void Build_Distribution(std::map<std::string,std::string> &keymap);
		void Set_Target_Parameters(std::map<std::string,std::string> &keymap);
};

#endif
