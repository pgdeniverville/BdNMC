#include "Axion_Dark_Photon.h"
#include "constants.h"
#include "branchingratios.h"
#include <cmath>

using std::string;
using std::vector;
using std::cout;
using std::endl;

Axion_Dark_Photon::Axion_Dark_Photon(Parameter* par){
    if(par->Model_Name()=="Axion_Dark_Photon"){
        Set_Model_Parameters(par);
    }
}

void Axion_Dark_Photon::Report(std::ostream& out, double tot){

    std::cout << "Lifetime of Dark Photon: " << lifetime << std::endl;
    out << "Lifetime " << lifetime << "\n";
    for(unsigned i=0; i<branching_ratios.size(); i++){
        out << final_states[i][0].name;
        for(unsigned j=1; j<final_states[i].size(); j++){
            out << "_" << final_states[i][j].name;
        }
        out << " " << branching_ratios[i]<< " " << mass_dp << " " << mass_axion << " " << epsilon << " " << eprime << " " << Gagg << " " << Gagpg << " " << Gagpgp;
        if(tot>=0){
            out << " " << tot*branching_ratios[i];
        }
        out << "\n";
    }
}

void Axion_Dark_Photon::Set_Model_Parameters(Parameter* par){
    using namespace Ax_DP;
    mass_axion = par->MassDM();
    mass_dp = par->MassDP();
    epsilon = par->Epsilon();
    eprime = sqrt(4*pi*par->alD());
    par->Query_Map("gagg", Gagg);
    par->Query_Map("gagpg", Gagpg);
    par->Query_Map("gagpgp", Gagpgp);
    //Currently assume that axion is stable, so dark photon
    //is the only interesting lifetime.
    double Gamma_tot = Gamma_dp(mass_dp, mass_axion, Gagpg, epsilon, eprime); 
    lifetime = hbar/Gamma_tot;
    
    Particle gamma(0);
    gamma.name = "Photon";

    Particle axion(mass_axion);
    axion.name = "Axion";

    Particle electron(MASS_ELECTRON);
    electron.name = "Electron";

    Particle muon(MASS_MUON);
    muon.name = "Muon";

    Particle hadronic(0);
    hadronic.name="Hadronic Stuff";

    branching_ratios.clear();
    final_states.clear();
    double Br;
    if((Br=Br_dp_to_a_gamma(mass_dp, mass_axion, Gagpg, epsilon, eprime))!=0){
        branching_ratios.push_back(Br);
        vector<Particle> vec;
        vec.push_back(gamma);
        vec.push_back(axion);
        final_states.push_back(vec);
    }
    if((Br=Br_dp_to_lepton(mass_dp, mass_axion, MASS_ELECTRON, Gagpg, epsilon, eprime))!=0){
        branching_ratios.push_back(Br);
        vector<Particle> vec;
        vec.push_back(electron);
        vec.push_back(electron);
        final_states.push_back(vec);
    }
    if((Br=Br_dp_to_lepton(mass_dp, mass_axion, MASS_MUON, Gagpg, epsilon, eprime))!=0){
        branching_ratios.push_back(Br);
        vector<Particle> vec;
        vec.push_back(muon);
        vec.push_back(muon);
        final_states.push_back(vec);
    }
    if((Br=Br_dp_to_3gamma(mass_dp, mass_axion, Gagpg, epsilon, eprime))!=0){
        branching_ratios.push_back(Br);
        vector<Particle> vec;
        vec.push_back(gamma);
        vec.push_back(gamma);
        vec.push_back(gamma);
        final_states.push_back(vec);
    }
    if((Br=Br_dp_to_hadrons(mass_dp, mass_axion, Gagpg, epsilon, eprime))!=0){
        branching_ratios.push_back(Br);
        vector<Particle> vec;
        vec.push_back(hadronic);
        final_states.push_back(vec);
    }
}


