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
        Model();
        ~Model(){};
        virtual void Set_Model_Parameters(Parameter& par) = 0;
        virtual void Report(std::ostream& out, double tot) = 0;
        virtual void Branching_Ratios();
        //These functions supply the prepared DMGenerator, Distribution and Scatter lists.
        virtual void DMGen(std::vector<std::shared_ptr<DMGenerator> > DMGen_list);
        virtual void Distribution(std::vector<std::shared_ptr<Distribution> > PartDist_list);
        virtual void SigGen(std::vector<std::shared_ptr<Scatter> > Signal_list);        
    private:
        std::vector<Decay_Channels> decay_channels;
        std::vector<std::shared_ptr<Scatter> > Sig_list;
        std::vector<std::shared_ptr<Distribution> > Dist_list;
        std::vector<std::shared_ptr<DMGenerator> > Gen_list;
};

class Pseudoscalar : public Model{
    public:
        Pseudoscalar();
        ~Pseudoscalar(){};
        void Set_Model_Parameters(Parameter& par);
        void Report(std::ostream& out, double tot){};
        void Branching_Ratios(){};
        void DMGen(std::vector<std::shared_ptr<DMGenerator> > DMGen_list);
        void Distribution(std::vector<std::shared_ptr<Distribution> > PartDist_list);
        void SigGen(std::vector<std::shared_ptr<Scatter> > Signal_list);        
    private :
            
}


#endif
