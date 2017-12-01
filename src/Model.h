#ifndef MODEL_H_GUARD
#define MODEL_H_GUARD

#include <vector>
#include <string>
#include <iostream>
#include <fstream>
#include <memory>

#include "Parameter.h"
#include "Particle.h""
#include "Scatter.h"
#include "DMGenerator.h"
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
    std::vector<vector<Particle> > end_states;
};

class Model{
    public:
        Model(Parameter& par){Model_Name=par->Model_Name(); Prepare_Model(par);}
        ~Model(){};
        //This is the only point at which the Model class gets access to a
        //Parameter object. It should use this opportunity to prepare DMGen_list, PartDist_list and Signal_list.
        void Prepare_Model(Parameter& par){
            if(!Set_Model_Parameters()){
                std::cerr << "Model " << Model_Name << " is missing required model parameters\n";
                throw -1;
            }
            std::cout << "Preparing signal channel " << par->Signal_Channel() << std::endl;
            scat_max = par->Max_Scatter_Energy();
            scat_min = par->Min_Scatter_Energy();
            if(!Prepare_Signal_Channel(par)){
                std::cerr << "Something is wrong with signal channel declaration\n";
                throw -1;
            }
        };
        virtual bool Set_Model_Parameters(Parameter& par) = 0;
        virtual bool Prepare_Signal_Channel(Parameter& par) = 0;
        virtual void Report(std::ostream& out, double tot) = 0;
        virtual void Branching_Ratios();
        //These functions supply the prepared DMGenerator, Distribution and Scatter lists.
        void get_DMGen(std::vector<std::shared_ptr<DMGenerator> > DMGen_list){DMGen_list = Gen_list;}
        void get_Distribution(std::vector<std::shared_ptr<Distribution> > PartDist_list){PartDist_list = Dist_list;}
        void get_SigGen(std::vector<std::shared_ptr<Scatter> > Signal_list){Signal_list = Sig_list;} 
        std::string model_name(){return Model_Name;}
    private:
        std::string Model_Name = "MODEL_NAME";
        std::vector<Decay_Channels> decay_channels;
        std::vector<std::shared_ptr<Scatter> > Sig_list;
        std::vector<std::shared_ptr<Distribution> > Dist_list;
        std::vector<std::shared_ptr<DMGenerator> > Gen_list;
        
        //Minimum and maximum scattering energies;
        double scat_max, scat_min;
};

class Pseudoscalar : public Model{
    public:
        Pseudoscalar(Parameter& par){Model_Name=par->Model_Name(); Prepare_Model(par)}
        ~Pseudoscalar(){};
        bool Set_Model_Parameters(Parameter& par);
        bool Prepare_Signal_Channel(Parameter& par);
        void Report(std::ostream& out, double tot){};
        void Branching_Ratios(){};
        //void get_DMGen(std::vector<std::shared_ptr<DMGenerator> > DMGen_list);
        //void get_Distribution(std::vector<std::shared_ptr<Distribution> > PartDist_list);
        //void get_SigGen(std::vector<std::shared_ptr<Scatter> > Signal_list);

        void dsigma_dEf_electron(double Ei, double Ef);
        void sigma_Ef_electron(double Ei, double Ef);
    private :
        //gchi is the dark matter charge, mchi is the dark matter mass, 
        //ma is the mass of the pseudoscalar mediator
        double gchi, ma, mchi;
        

}


#endif
