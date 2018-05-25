#ifndef GUARD_Proton_Brem_h
#define GUARD_Proton_Brem_h

#include "Parameter.h"
#include "DMgenerator.h"
#include <functional>
#include <string>

class Proton_Brem : public DMGenerator{
        
    public:
        Proton_Brem(double Beam_E, std::function<double(double, double)> splitting_function, Particle &mediator, double ptmax, double zmax, double zmin, std::string &mode, std::shared_ptr<DMGenerator> V_decay, double ptmin=0);
        bool GenDM(std::list<Particle>& vec, std::function<double(Particle&)> det_int, Particle& part);
        double V_prod_rate(){return vprodrate*V_decay->BranchingRatio();}
        void sample_particle(Particle &);
        void sample_momentum(double &, double &, double &);
        void set_fit_parameters(production_channel &);
        double d2N_proton_brem_to_out(double z, double pt2);
    private:
        std::shared_ptr<DMGenerator> V_decay;
        void calc_V_prod_rate();
        double sigmapp(double s);
        std::string model;//proton_brem_baryonic triggers leptophobic behavior
        // double F_1_proton(double q2);
            
        double Beam_Energy, PTMIN, PTMAX, ZMAX, ZMIN;
        double vprodrate, max_prod;
        double MA;
        std::string med_name;
            
        //d^2sigma/dpt^2/dz
        std::function<double(double, double)> dsig;

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