#ifndef GUARD_DMgenerator_h
#define GUARD_DMgenerator_h

#include "Particle.h"
#include <string>
#include <vector>
#include <list>
#include <functional>
#include <memory>
//I need DMGenerator to provide methods that can be called by its child functions to implement commonly used code. I should make a DMGenerator.cpp file to hold it.
class DMGenerator{
    public:
        DMGenerator(){}
        virtual bool GenDM(std::list<Particle>& vec, std::function<double(Particle)> det_int, Particle& part) {return false;}
        double BranchingRatio(){return branchingratio;}
        std::string Channel_Name(){return chan_name;}
        void set_model_params(double MV, double MX, double kap, double alp){mv=MV; mx=MX; kappa=kap; alphaD=alp; Evaluate_Branching_Ratio();}
        bool query_off_shell(){return OFF_SHELL;}
        virtual ~DMGenerator(){};
    protected:
        virtual void Evaluate_Branching_Ratio() {};
        double mv, mx, kappa, alphaD;
        bool OFF_SHELL;
        std::string chan_name;
        double branchingratio;
        double pmax;
};

class pion_decay_gen: public DMGenerator{
    public:
        pion_decay_gen(double MV, double MX, double kap, double alp);
        ~pion_decay_gen(){};
        bool GenDM(std::list<Particle>& vec, std::function<double(Particle)> det_int, Particle& part);
        void sample_dist(double& s, double& theta);
        void burn_in(int runs);
    private:
        void Evaluate_Branching_Ratio();
};

class eta_decay_gen: public DMGenerator{
    public:
        eta_decay_gen(double MV, double MX, double kap, double alp);
        ~eta_decay_gen(){};
        bool GenDM(std::list<Particle>& vec, std::function<double(Particle)> det_int, Particle& part);
        void sample_dist(double& s, double& theta);
        void burn_in(int runs);
    private:
        void Evaluate_Branching_Ratio();
};

class pion_decay_gen_baryonic: public DMGenerator{
    public:
        pion_decay_gen_baryonic(double MV, double MX, double kap, double alp);
        bool GenDM(std::list<Particle>& vec, std::function<double(Particle)> det_int, Particle& part);
        void sample_dist(double& s, double& theta);
        void burn_in(int runs);
    private:
        void Evaluate_Branching_Ratio();
};

class eta_decay_gen_baryonic: public DMGenerator{
    public:
        eta_decay_gen_baryonic(double MV, double MX, double kap, double alp);
        bool GenDM(std::list<Particle>& vec, std::function<double(Particle)> det_int, Particle& part);
        void sample_dist(double& s, double& theta);
        void burn_in(int runs);
    private:
        void Evaluate_Branching_Ratio();
};


class piminus_capture_gen: public DMGenerator{
    public:
        piminus_capture_gen(double MV, double MX, double kap, double alp);
        bool GenDM(std::list<Particle>& vec, std::function<double(Particle)> det_int, Particle& part);
        void sample_dist(double& s, double& theta);
        void burn_in(int runs);
    private:
        void Evaluate_Branching_Ratio();
};

class rho_decay_gen: public DMGenerator{
    public:
        rho_decay_gen(double MV, double MX, double kap, double alp);
        bool GenDM(std::list<Particle>& vec, std::function<double(Particle)> det_int, Particle& part);
    private:
        void Evaluate_Branching_Ratio();
};

class phi_decay_gen_baryonic: public DMGenerator{
    public:
        phi_decay_gen_baryonic(double MV, double MX, double kap, double alp);
        bool GenDM(std::list<Particle>& vec, std::function<double(Particle)> det_int, Particle& part);
    private:
        void Evaluate_Branching_Ratio();
};

class phi_decay_gen: public DMGenerator{
    public:
        phi_decay_gen(double MV, double MX, double kap, double alp);
        bool GenDM(std::list<Particle>& vec, std::function<double(Particle)> det_int, Particle& part);
    private:
        void Evaluate_Branching_Ratio();
};


class omega_decay_gen: public DMGenerator{
    public:
        omega_decay_gen(double MV, double MX, double kap, double alp);
        bool GenDM(std::list<Particle>& vec, std::function<double(Particle)> det_int, Particle& part);
    private:
        void Evaluate_Branching_Ratio();
};

class omega_decay_gen_baryonic: public DMGenerator{
    public:
        omega_decay_gen_baryonic(double MV, double MX, double kap, double alp);
        bool GenDM(std::list<Particle>& vec, std::function<double(Particle)> det_int, Particle& part);
    private:
        void Evaluate_Branching_Ratio();
};

class parton_V_gen: public DMGenerator{
    public:
        parton_V_gen(double MV, double MX, double kap, double alp, const std::string chan);
        bool GenDM(std::list<Particle>& vec, std::function<double(Particle)> det_int, Particle& part);
    private:
        void Evaluate_Branching_Ratio();
};

class V_decay_gen: public DMGenerator{
	public:
		V_decay_gen(double MV, double MX, double kap, double alp, const std::string chan=""); 
        bool GenDM(std::list<Particle>& vec, std::function<double(Particle)> det_int, Particle& part);
	private:
		void Evaluate_Branching_Ratio();
};

#endif