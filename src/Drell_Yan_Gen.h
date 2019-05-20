#ifndef GUARD_DYG_H
#define GUARD_DYG_H

#include "DMgenerator.h"
#include "cteq_pdf_reader/ct11pdf.h"
#include "Integrator.h"
#include "Parameter.h"

//This class describes pp -> \chi \chi + X.
//Currently assuming a mono-energetic beam.

class Drell_Yan_Gen: public DMGenerator{
    public:
        Drell_Yan_Gen(double mchi, double BeamEnergy, std::function<double(double, double, double, double, double)> dsigma_hat, double target_p, double target_n, production_channel& prodchan); 
        //const std::string pdf_file_path= "", const std::string pdf_neutron_file_path = "", const std::string chan = "", double Tolerance = 5e-2, double energy_bins = 100, double);
        bool GenDM(std::list<Particle>& vec, std::function<double(Particle&)> det_int, Particle& part);
        void set_product_name(const std::string pn){product_name = pn;}
        //double Cross_Section(){return cross_section;}
        //This returns (sig_p*p_num+sig_n*n_num)
        double Interaction_Cross_Section();
        double Sig_P(){return sig_p;}
    private:
        double sig(double x, double y, int mode);

        double dsig(double Ek, double x, double y, double Q, int mode);
        
        double Monte_Carlo(int mode);
        double Monte_Carlo(double Ek, int mode, double& max);
        void Evaluate_Branching_Ratio() {Evaluate_Cross_Sections();}
        void Evaluate_Cross_Sections();
        double Q2_eval(double E1, double m_target, double x, double y);
        bool Verify(double Ek, double m_target, double x, double y, double &Q);
        double Y_Min_from_QMIN(double E1, double m1, double m2, double x);
        std::function<double(double, double, double, double, double)> dsig_hat;
        std::function<double(double, double, double, double)> sig_hat;
        Linear_Interpolation2 dsigma_dE_proton;
        Linear_Interpolation2 dsigma_dE_neutron;
        Linear_Interpolation dsigma_dE_proton_max;
        Linear_Interpolation dsigma_dE_neutron_max;
        std::string product_name="DM";
        double BeamEnergy, Efmax, Efmin, Efres, mx, sig_p, sig_n, QMIN;
        double pMax_n, pMax_p;
        cteqpdf proton_pdf;
        //This is not being used yet.
        cteqpdf neutron_pdf;
        double TOLERANCE; 
        double target_n, target_p;

        double sig_test(double x, double y, int mode);
        double Monte_Carlo_test(int mode);
        double Monte_Carlo_test_2(int mode);
        double dsigma_dt(double x, double y, double mass_target, double shat);
        double sig_tot(double x, double y, double mass_target);
        void Test();
};

#endif
