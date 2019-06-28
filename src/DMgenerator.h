#ifndef GUARD_DMgenerator_h
#define GUARD_DMgenerator_h

#include "Particle.h"
#include <string>
#include <vector>
#include <list>
#include <functional>
#include <memory>
//I need DMGenerator to provide methods that can be called by its subclasses to implement commonly used code. I should make a DMGenerator.cpp file to hold it.
class DMGenerator{
    public:
        DMGenerator(){}
        virtual bool GenDM(std::list<Particle>& vec, std::function<double(Particle&)> det_int, Particle& part) { return false;}
        double BranchingRatio(){return branchingratio;}
        void set_BranchingRatio(double br){branchingratio=br;}
        void Set_Channel_Name(std::string ch_na){chan_name=ch_na;}
        std::string Channel_Name(){return chan_name;}
        void set_model_params(double MV, double MX, double kap, double alp){mv=MV; mx=MX; kappa=kap; alphaD=alp; Evaluate_Branching_Ratio();}
        bool query_off_shell(){return OFF_SHELL;}
        void set_Off_Shell(bool u){OFF_SHELL=u;}
        virtual ~DMGenerator(){};
        bool record_parent = true;
    protected:
        virtual void Evaluate_Branching_Ratio() {};
        double mv, mx, kappa, alphaD;
        bool OFF_SHELL = false;
        std::string chan_name;
        double branchingratio;
        double pmax;
};

class pion_decay_gen: public DMGenerator{
    public:
        pion_decay_gen(double MV, double MX, double kap, double alp);
        ~pion_decay_gen(){};
        bool GenDM(std::list<Particle>& vec, std::function<double(Particle&)> det_int, Particle& part);
        void sample_dist(double& s, double& theta);
        void burn_in(int runs);
    private:
        void Evaluate_Branching_Ratio();
};

class eta_decay_gen: public DMGenerator{
    public:
        eta_decay_gen(double MV, double MX, double kap, double alp);
        ~eta_decay_gen(){};
        bool GenDM(std::list<Particle>& vec, std::function<double(Particle&)> det_int, Particle& part);
        void sample_dist(double& s, double& theta);
        void burn_in(int runs);
    private:
        void Evaluate_Branching_Ratio();
};

class pion_decay_gen_baryonic: public DMGenerator{
    public:
        pion_decay_gen_baryonic(double MV, double MX, double kap, double alp);
        bool GenDM(std::list<Particle>& vec, std::function<double(Particle&)> det_int, Particle& part);
        void sample_dist(double& s, double& theta);
        void burn_in(int runs);
    private:
        void Evaluate_Branching_Ratio();
};

class eta_decay_gen_baryonic: public DMGenerator{
    public:
        eta_decay_gen_baryonic(double MV, double MX, double kap, double alp);
        bool GenDM(std::list<Particle>& vec, std::function<double(Particle&)> det_int, Particle& part);
        void sample_dist(double& s, double& theta);
        void burn_in(int runs);
    private:
        void Evaluate_Branching_Ratio();
};


class piminus_capture_gen: public DMGenerator{
    public:
        piminus_capture_gen(double MV, double MX, double kap, double alp);
        bool GenDM(std::list<Particle>& vec, std::function<double(Particle&)> det_int, Particle& part);
        void sample_dist(double& s, double& theta);
        void burn_in(int runs);
    private:
        void Evaluate_Branching_Ratio();
};

class rho_decay_gen: public DMGenerator{
    public:
        rho_decay_gen(double MV, double MX, double kap, double alp);
        bool GenDM(std::list<Particle>& vec, std::function<double(Particle&)> det_int, Particle& part);
    private:
        void Evaluate_Branching_Ratio();
};

class phi_decay_gen_baryonic: public DMGenerator{
    public:
        phi_decay_gen_baryonic(double MV, double MX, double kap, double alp);
        bool GenDM(std::list<Particle>& vec, std::function<double(Particle&)> det_int, Particle& part);
    private:
        void Evaluate_Branching_Ratio();
};

class phi_decay_gen: public DMGenerator{
    public:
        phi_decay_gen(double MV, double MX, double kap, double alp);
        bool GenDM(std::list<Particle>& vec, std::function<double(Particle&)> det_int, Particle& part);
    private:
        void Evaluate_Branching_Ratio();
};


class omega_decay_gen: public DMGenerator{
    public:
        omega_decay_gen(double MV, double MX, double kap, double alp);
        bool GenDM(std::list<Particle>& vec, std::function<double(Particle&)> det_int, Particle& part);
    private:
        void Evaluate_Branching_Ratio();
};

class omega_decay_gen_baryonic: public DMGenerator{
    public:
        omega_decay_gen_baryonic(double MV, double MX, double kap, double alp);
        bool GenDM(std::list<Particle>& vec, std::function<double(Particle&)> det_int, Particle& part);
    private:
        void Evaluate_Branching_Ratio();
};

class parton_V_gen: public DMGenerator{
    public:
        parton_V_gen(double MV, double MX, double kap, double alp, const std::string chan);
        bool GenDM(std::list<Particle>& vec, std::function<double(Particle&)> det_int, Particle& part);
    private:
        void Evaluate_Branching_Ratio();
};

class V_decay_gen: public DMGenerator{
	public:
		V_decay_gen(double MV, double MX, double kap, double alp, const std::string chan=""); 
        bool GenDM(std::list<Particle>& vec, std::function<double(Particle&)> det_int, Particle& part);
	private:
		void Evaluate_Branching_Ratio();
};

class Do_Nothing_Gen: public DMGenerator{
    public:
        Do_Nothing_Gen(const std::string chan="", const std::string part_name="");
        bool GenDM(std::list<Particle>& vec, std::function<double(Particle&)> det_int, Particle& part);
        void Set_Part_Name(const std::string part_name){Part_Name = std::string(part_name);}
    private:
        void Evaluate_Branching_Ratio();
        std::string Part_Name="";
};

//Defined in V_decay_gen.cpp
//Need to add unstable daughter particles
class Two_Body_Decay_Gen: public DMGenerator{
    public:
        Two_Body_Decay_Gen(double branching_ratio, double parent_mass, const std::string parent_name, Particle daughter1, Particle daughter2, double lifetime);
        Two_Body_Decay_Gen(double branching_ratio, double parent_mass, const std::string parent_name, Particle daughter1, Particle daughter2);
        bool GenDM(std::list<Particle>& vec, std::function<double(Particle&)> det_int, Particle& part);
        void Toggle_Daughter_Decay(int index, std::shared_ptr<DMGenerator> daughter_decay_gen);
        //If true, then this daughter is checked for intersection with the detector.
        bool d1 = true;
        bool d2 = true;
    private:
        void Evaluate_Branching(){return;}
        bool d1_unstable = false;
        bool d2_unstable = false;
        std::shared_ptr<DMGenerator> d1_decay;
        std::shared_ptr<DMGenerator> d2_decay;
        Particle daughter1, daughter2;
        std::string Part_Name;
        double tau;
};

//This is nearly complete!
//If I go above 3-body final states, it would make sense to turn everything into lists.
class Three_Body_Decay_Gen: public DMGenerator{
    public:
        //std::shared_ptr<DMGenerator> decaygen1, std::shared_ptr<DMGenerator> decaygen2, std::shared_ptr<DMGenerator> decaygen3, will be added in when I add decays of component particles.
        Three_Body_Decay_Gen(Particle& parent, Particle& daughter1, Particle& daughter2, Particle& daughter3, std::string prodchoice,  double lifetime, std::function<double(double, double, double, double, double, double)> &amp);
        Three_Body_Decay_Gen(Particle& parent, Particle& daughter1, Particle& daughter2, Particle& daughter3, std::string prodchoice,  double lifetime, double partial_width, std::function<double(double, double, double, double, double, double)> &amp, int burn=1000);
        bool GenDM(std::list<Particle>& vec, std::function<double(Particle&)> det_int, Particle& part);
        //Turns on decays for daughter 1, 2 or 3.
        void Toggle_Daughter_Decay(int index, std::shared_ptr<DMGenerator> daughter_decay_gen);
        //Are these particles checked for intersection with the detector?
        void Burn_In(int i);
        bool d1 = true;
        bool d2 = true;
        bool d3 = true;
        std::string prodchoice;
    private:
        void Evaluate_Branching(){return;}
        Particle mother, daughter1, daughter2, daughter3;
        //Lifetime of particle, not yet implemented.
        double tau;
        std::function<double(double, double, double, double, double, double)> amp;
        bool d1_unstable = false;
        bool d2_unstable = false;
        bool d3_unstable = false;
        std::shared_ptr<DMGenerator> d1_decay;
        std::shared_ptr<DMGenerator> d2_decay;
        std::shared_ptr<DMGenerator> d3_decay;
};

#endif
