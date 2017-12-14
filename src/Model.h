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

class Model{
    public:
        Model(Parameter& par){Model_Name=par.Model_Name(); Vnumtot=0; Prepare_Model(par);}
        ~Model(){};
        //This is the only point at which the Model class gets access to a
        //Parameter object. It should use this opportunity to prepare DMGen_list, PartDist_list and Signal_list.
        //Should probably move this implementation into a Model.cpp
        void Prepare_Model(Parameter& par);
        virtual bool Set_Model_Parameters(Parameter& par) = 0;
        virtual bool Prepare_Signal_Channel(Parameter& par) = 0;
        virtual bool Prepare_Production_Channel(std::string prodchoice, std::string proddist, production_channel& prodchan, std::shared_ptr<DMGenerator> DMGen, std::shared_ptr<Distribution>, double& Vnum, Parameter& par) = 0;
        //Model can handle this itself, since Distributions use SM physics.
        bool Prepare_Production_Distribution(std::string prodchoice, std::string proddist, production_channel& prodchan, std::shared_ptr<Distribution> Dist, Parameter& par);
        virtual void Report(std::ostream& out, double tot) = 0;
        virtual void Branching_Ratios() = 0;
        //These functions supply the prepared DMGenerator, Distribution and Scatter lists.
        void get_DMGen(std::vector<std::shared_ptr<DMGenerator> > DMGen_list){DMGen_list = Gen_list;}
        void get_Distribution(std::vector<std::shared_ptr<Distribution> > PartDist_list){PartDist_list = Dist_list;}
        void get_SigGen(std::vector<std::shared_ptr<Scatter> > Signal_list){Signal_list = Sig_list;}
        void get_first_SigGen(std::shared_ptr<Scatter> Signal){Signal = Sig_list.front();}
        void get_Vnum(std::vector<double> vnum){vnum = Vnum_list;}
        double get_Vnumtot(){return Vnumtot;}
        std::string model_name(){return Model_Name;}
    protected:
        std::string Model_Name = "MODEL_NAME";
        std::vector<Decay_Channels> decay_channels;
        std::vector<std::shared_ptr<Scatter> > Sig_list;
        std::vector<std::shared_ptr<Distribution> > Dist_list;
        std::vector<std::shared_ptr<DMGenerator> > Gen_list;
        std::vector<double> Vnum_list;
        //Minimum and maximum scattering energies;
        double scat_max, scat_min,Vnumtot;
};

class Pseudoscalar : public Model{
    public:
        //Calling Model constructor, it does everything we need.
        Pseudoscalar(Parameter& par) : Model(par){};
        //~Pseudoscalar(){};
        //bool Set_Model_Parame ters(Parameter& par);
        bool Prepare_Signal_Channel(Parameter& par);
        bool Prepare_Production_Channel(std::string prodchoice, std::string proddist, production_channel& prodchan, std::shared_ptr<DMGenerator> DMGen, std::shared_ptr<Distribution>, double& Vnum, Parameter& par);
        bool Set_Model_Parameters(Parameter& par);
        void Report(std::ostream& out, double tot){};
        void Branching_Ratios(){};
        //void get_DMGen(std::vector<std::shared_ptr<DMGenerator> > DMGen_list);
        //void get_Distribution(std::vector<std::shared_ptr<Distribution> > PartDist_list);
        //void get_SigGen(std::vector<std::shared_ptr<Scatter> > Signal_list);

        double dsigma_dEf_electron(double Ei, double Ef);
        double sigma_Ef_electron(double Ei, double Ef);
        double dsigma_dEk_qq_to_chichi(double Ek, double EA, double x,double y, double gf, double MASS);
        //p_1 = x P_A, p_2 = y P_B, E_k is outgoing energy of p_3 (one of the 
        //chi particles).
       
    private :
        //gchi is the dark matter charge, mchi is the dark matter mass, 
        //ma is the mass of the pseudoscalar mediator
        double gchi, ma, mchi;
        double dsig_max(double Ei);
        double sigma_tot_electron(double Ei);

};


#endif
