#ifndef MODEL_H_GUARD
#define MODEL_H_GUARD

#include <vector>
#include <string>
#include <iostream>
#include <fstream>
#include <memory>

#include "Parameter.h"
#include "Particle.h"
#include "Scatter.h"
#include "DMgenerator.h"
#include "Distribution.h"

//The decays of each particle in the scenario are captured by
//a Decay_Channels object. These can be used by Signal_Decay
//later on.
struct Decay_Channels{
    Decay_Channels();
    //A negative lifetime will indicate stability.
    double lifetime = -1;
    std::string name = "";
    std::vector<double> partial_widths;
    //Each vector of particles provides one possible decay channel
    std::vector<std::vector<Particle> > end_states;
};

double load_form_factor(const std::string& filename, std::shared_ptr<Linear_Interpolation>& ff);

class Model{
    public:
        Model(Parameter& par){Model_Name=par.Model_Name(); Vnumtot=0;}
        ~Model(){};
        //This is the only point at which the Model class gets access to a
        //Parameter object. It should use this opportunity to prepare DMGen_list, PartDist_list and Signal_list.
        //Should probably move this implementation into a Model.cpp
        void Prepare_Model(Parameter& par);
        void Report_Model(std::ostream);
        virtual bool Set_Model_Parameters(Parameter& par) = 0;
        virtual bool Prepare_Signal_Channel(Parameter& par) = 0;
        virtual bool Prepare_Production_Channel(std::string prodchoice, std::string proddist, production_channel& prodchan, std::shared_ptr<DMGenerator>& DMGen, std::shared_ptr<Distribution>&, double& Vnum, Parameter& par) = 0;
        //Model can handle this itself, since Distributions use SM physics.
        bool Prepare_Production_Distribution(std::string prodchoice, std::string proddist, production_channel& prodchan, std::shared_ptr<Distribution>& Dist, Parameter& par);
        virtual void Report(std::ostream& out) = 0;
        virtual void Report_Model() = 0;
        virtual void Branching_Ratios() = 0;
        //These functions supply the prepared DMGenerator, Distribution and Scatter lists.
        void get_DMGen(std::vector<std::shared_ptr<DMGenerator> >& DMGen_list){DMGen_list = Gen_list;}
        std::vector<std::shared_ptr<DMGenerator> > get_DMGen(){return Gen_list;}
        void get_Distribution(std::vector<std::shared_ptr<Distribution> >& PartDist_list){PartDist_list = Dist_list;}
        std::vector<std::shared_ptr<Distribution> > get_Distribution(){return Dist_list;}
        void get_SigGen(std::vector<std::shared_ptr<Scatter> >& Signal_list){Signal_list = Sig_list;}
        void get_first_SigGen(std::shared_ptr<Scatter>& Signal){if(Sig_list.size()==0){std::cerr << "Sig_list is empty, get_first_SigGen will fail\n";}
            Signal = Sig_list.front();}
        std::shared_ptr<Scatter> get_SigGen(int i){return Sig_list[i];}
        void get_Vnum(std::vector<double>& vnum){vnum = Vnum_list;}
        std::vector<double> get_Vnum(){return Vnum_list;}
        double get_Vnumtot(){return Vnumtot;}
        std::string model_name(){return Model_Name;}
        std::vector<std::string> get_sig_part_vec(){return sig_part_vec;}
    protected:
        //std::string sig_part_name="DM";
        std::vector<std::string> sig_part_vec;
        std::string Model_Name = "MODEL_NAME";
        std::vector<Decay_Channels> decay_channels;
        std::vector<std::shared_ptr<Scatter> > Sig_list;
        std::vector<std::shared_ptr<Distribution> > Dist_list;
        std::vector<std::shared_ptr<DMGenerator> > Gen_list;
        std::vector<double> Vnum_list;
        //Minimum and maximum scattering energies;
        double scat_max, scat_min,Vnumtot;
};

class Inelastic_Dark_Matter : public Model{
    public:
        Inelastic_Dark_Matter(Parameter& par) : Model(par){};
        bool Prepare_Signal_Channel(Parameter& par);
        bool Prepare_Production_Channel(std::string prodchoice, std::string proddist, production_channel& prodchan, std::shared_ptr<DMGenerator>& DMGen, std::shared_ptr<Distribution>&, double& Vnum, Parameter& par);
        bool Set_Model_Parameters(Parameter& par);
        void Report_Model();
        void Report(std::ostream& out);
        void Branching_Ratios(){};
        double meson_decay_amplitude2(double m12s, double m23s, double m0, double m1, double m2, double m3);
        double Amplitude2_dm2_to_lepton_lepton_dm1(double m12s, double m23s, double m0, double m1, double m2, double m3);
        //m1=m2=MASS_MUON
        double Amplitude2_dm2_to_hadrons_dm1(double m12s, double m23s, double m0, double m3);
        double A_width();
        double dm2_width();
    private:
        double dsigma_dm_e_to_dm_e(double E1lab, double E4, double mass_dm_in,double mass_dm_out, double mR);
        double dm_e_to_dm_e_amp(double s, double t, double mass_dm_in, double mass_dm_out, double mR);
        double dsigma_dm_e_to_dm_e_2(double E1lab, double t, double mass_dm_in,double mass_dm_out, double mR);
        double dsigma_dm_N_to_dm_Delta(double E1lab, double E4, double mass_dm_in, double mass_dm_out, double mN, double delta);
        double amp_dm_N_to_dm_Delta(double s, double t, double m1, double m2, double mN, double mD);
        void Evaluate_Widths();
        double GM(double q2);
        //mass_dm2 is assumed to be heavier.
        std::shared_ptr<Linear_Interpolation> GM_form_factor;
        bool MAX_Q2_WARNING = false;
        double MAX_Q2=0;
        double Awidth, dm2width;
        double width_dm2_to_elec_elec_dm1;
        double width_dm2_to_muon_muon_dm1;
        double width_dm2_to_dm1_hadrons;
        double mass_dp, mass_dm1, mass_dm2, epsilon, alpha_D;
};

class Axion_Dark_Photon : public Model{
    public:
        //Might not be enough, we'll see.
        Axion_Dark_Photon(Parameter& par) : Model(par){};
        bool Prepare_Signal_Channel(Parameter& par);
        bool Prepare_Production_Channel(std::string prodchoice, std::string proddist, production_channel& prodchan, std::shared_ptr<DMGenerator>& DMGen, std::shared_ptr<Distribution>&, double& Vnum, Parameter& par);
        bool Set_Model_Parameters(Parameter& par);
        void Report_Model();
        void Report(std::ostream& out);
        void Branching_Ratios(){};
    private:
        double mass_axion, mass_dp, Gagg, Gagpg, Gagpgp, epsilon;
        double pi0_decay_amplitude2(double m12s, double s, double mgamma, double ma, double mA);
        double pi0_decay_amplitude2_b(double m12, double CosTheta, double m0, double m1, double m2, double m3);
        double pi0_decay_amplitude2_c(double m12s, double s, double m0, double m1, double m2, double m3);
        double pi0_decay_amplitude2(double m12, double CosTheta, double m0, double m1, double m2, double m3);
        double pi0_decay_width(double m12s, double s, double mgamma, double ma, double mA);
        double pi0_decay_width_2(double m12, double ct1, double mpion, double mgamma, double maxion, double mdarkphoton);

        double eta_decay_amplitude2_c(double m12s, double s, double m0, double m1, double m2, double m3);
        double meson_decay_amplitude2(double m12s, double s, double mmeson, double mgamma, double ma, double mA);
        double eta_decay_width_2(double m12, double CosTheta, double m0, double m1, double m2, double m3);

        double A_to_a_gamma_width(double ma, double mA);
        double A_to_a_elec_pos_amplitude(double m12s, double m23s, double mA, double me, double me2, double ma);
        double A_width(double ma, double mA);


        double Axion_DP_electron_Amp(double s, double t, double mA, double mR);
        double dsigma_a_to_DP(double E1lab, double E4, double mA, double mR);
        double dsigma_DP_to_a(double E1lab, double E4, double mA, double mR);
};

class Pseudoscalar : public Model{
    public:
        //Calling Model constructor, it does everything we need.
        Pseudoscalar(Parameter& par) : Model(par){};
        //~Pseudoscalar(){};
        //bool Set_Model_Parame ters(Parameter& par);
        bool Prepare_Signal_Channel(Parameter& par);
        bool Prepare_Production_Channel(std::string prodchoice, std::string proddist, production_channel& prodchan, std::shared_ptr<DMGenerator>& DMGen, std::shared_ptr<Distribution>&, double& Vnum, Parameter& par);
        bool Set_Model_Parameters(Parameter& par);
        void Report_Model();
        void Report(std::ostream& out);
        void Branching_Ratios(){};
        //void get_DMGen(std::vector<std::shared_ptr<DMGenerator> > DMGen_list);
        //void get_Distribution(std::vector<std::shared_ptr<Distribution> > PartDist_list);
        //void get_SigGen(std::vector<std::shared_ptr<Scatter> > Signal_list);

        double dsigma_dEf_electron(double Ei, double Ef);
        double sigma_Ef_electron(double Ei, double Ef);
        double dsigma_dEk_qq_to_chichi(double Ek, double EA, double x,double y, double gf, double MASS);
        double dsigma_hat_dt_qq_to_chi_chi(double EA, double t, double x, double y, double qf, double MASS);
        double sigma_hat_tot_qq_to_chi_chi(double EA, double x, double y, double qf, double MASS);
        //p_1 = x P_A, p_2 = y P_B, E_k is outgoing energy of p_3 (one of the 
        //chi particles).
        double proton_form_factor(double q2);

    private :
        //gchi is the dark matter charge, mchi is the dark matter mass, 
        //ma is the mass of the pseudoscalar mediator
        //gq is the A-quark coupling, qe is the A-electron coupling.
        double gchi, gq, gae, ma, mchi;
        double gnu=1;
        double dsig_max(double Ei);
        double sigma_tot_electron(double Ei);

};

class Kinetic_Mixing : public Model{
    public:
        Kinetic_Mixing(Parameter& par) : Model(par){};
        bool Prepare_Signal_Channel(Parameter& par);
        bool Prepare_Production_Channel(std::string prodchoice, std::string proddist, production_channel& prodchan, std::shared_ptr<DMGenerator>& DMGen, std::shared_ptr<Distribution>&, double& Vnum, Parameter& par);
        bool Set_Model_Parameters(Parameter& par);
        void Report_Model();
        void Report(std::ostream& out);
        void Branching_Ratios(){};
    private:
        //model parameters
        double alpha_D, mv, mx, eps;
};

#endif
