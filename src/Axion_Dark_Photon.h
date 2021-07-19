#ifndef AXION_DARK_PHOTON_GUARD
#define AXION_DARK_PHOTON_GUARD

#include "Particle.h"
#include "Parameter.h"

#include <vector>
#include <string>
#include <iostream>
#include <fstream>

/*
 *Axion_Dark_Photon is a class that prepares Axion Dark Photon physics for
 *the simulation 
 */
class Axion_Dark_Photon {
    public:
        Axion_Dark_Photon(Parameter* par);
        double Lifetime(){return lifetime;}
        void Branching_Ratios(std::vector<double> &br){br=branching_ratios;}
        void Final_States(std::vector<std::vector<Particle> >& fl){fl=final_states;}
        void Set_Model_Parameters(Parameter* par);
        //Axion-Photon-Photon
        double GAGG(){return Gagg;}
        //Axion-Photon-Dark Photon
        double GAGPG(){return Gagpg;}
        //Axion-Dark Photon-Dark Photon
        double GAGPGP(){return Gagpgp;}
        //Axion-Electron-Electron
        double GAEE(){return Gaee}
        //This will be expanded later to accept string arguments and
        //provide specific reports based on channel.
        void Report(std::ostream&, double tot=0);
    private:
        double mass_axion;
        double mass_dp;
        double epsilon;
        double eprime;
        double Gagg, Gagpg, Gagpgp;
        double Gaee=0;

        double lifetime;
        std::vector<double> branching_ratios;
        std::vector<std::vector<Particle> > final_states;
};

#endif
