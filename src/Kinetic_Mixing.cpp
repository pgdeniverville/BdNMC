/* 
*This file will hold all the standard Dark Photon with scalar DM stuff.
*I will be migrating everything from main to here eventually.
*/

#include "Model.h"
#include "Distribution.h"
#include "Proton_Brem_Distribution.h"
#include "Position_Distributions.h"
#include "branchingratios.h"
#include "Scatter.h"
#include <cmath>
#include "decay.h"
#include "DMgenerator.h"
#include "Model.h"
#include "Kinematics.h"
#include "SignalDecay.h"
#include "constants.h"

using std::string;
using std::vector;
using std::cout; using std::cerr;
using std::endl;
using std::bind;
using std::function;

using namespace std::placeholders;

//The ratio at which off-shell mode is activated
const double off_shell_ratio=1.3;

const string GM_form_factor_filename="data/delta_production_form_factor.dat";

bool Kinetic_Mixing::Set_Model_Parameters(Parameter& par){
	if(!(par.Query_Map("dark_photon_mass",mass_dp))) {
		std::cerr << par.Model_Name() << " requires " << "dark_photon_mass" << " to be defined.\n";
        return false;
	}
	if(!(par.Query_Map("dark_matter_mass",mass_dm))) {
		std::cerr << par.Model_Name() << " requires " << "dark_matter_mass" << " to be defined.\n";
        return false;
	}
	if(!(par.Query_Map("alpha_D",alpha_D))) {
		std::cerr << par.Model_Name() << " requires " << "alpha_D" << " to be defined.\n";
        return false;
	}
	if(!(par.Query_Map("epsilon",epsilon)) and !(par.Query_Map("kappa",epsilon))) {
		std::cerr << par.Model_Name() << " requires " << "epsilon or kappa" << " to be defined.\n";
        return false;
	}
	return true;
}

void Kinetic_Mixing::Report(std::ostream& out){
	out << mass_dp << " " << mass_dm << " " << epsilon << " " << alpha_D << " ";
}

void Kinetic_Mixing::Report_Model(){
    cout << "Mass Dark Photon = " << mass_dp << " GeV" << endl;
    cout << "Mass Dark Matter = " << mass_dm << " GeV" << endl;
    cout << "epsilon = " << epsilon << endl;
    cout << "alpha_D = " << alpha_D << endl;
}

bool Kinetic_Mixing::Prepare_Production_Channel(std::string prodchoice, std::string proddist, 
	production_channel& prodchan, std::shared_ptr<DMGenerator>& DMGen, 
	std::shared_ptr<Distribution>& Dist, double& Vnum, Parameter& par){

	Particle photon(0);
    photon.name="Photon";

    Particle dark_photon(mass_dp);
    dark_photon.name="Dark_Photon";

    Particle dm1(mass_dm);
    dm1.name = "DM";

    Particle dm2(mass_dm);
    dm1.name = "DM";

    string sigchoice = par.Signal_Channel();

    double POT = par.Protons_on_Target();
    double num_pi0 = par.Pi0_per_POT()*POT;

    if(sigchoice=="Signal_Decay"){
        if(proddist=="proton_brem"){
			string signal_string = "Dark_Photon";//This thing is getting killed as soon as possible.
            DMGen = std::shared_ptr<DMGenerator>(new Do_Nothing_Gen("Dark_Photon_Bremsstrahlung", signal_string));
            cout << DMGen->Channel_Name() << endl;
            std::shared_ptr<Proton_Brem_Distribution> pbd(new Proton_Brem_Distribution(par.Beam_Energy(),epsilon,mass_dp,prodchan.ptmax(),prodchan.zmax(),prodchan.zmin(),alpha_D,proddist,prodchan.ptmin()));
            Vnum = pbd->V_prod_rate()*par.Protons_on_Target();//V_decay should return a branching ratio for this. If it doesn't, I will have to change it in the future. prodchan.Num_per_pot()
            Dist->set_mass(mass_dp);
            Dist->Add_Dist(pbd);
            return true;
        }
        else if(prodchoice=="pi0_decay"||prodchoice=="eta_decay"){
            //Code repeat here. Need to reorganize to eliminate this.
            if(prodchoice=="pi0_decay"){
                if(mass_dp>MASS_PION_NEUTRAL){
                    cerr << "Attempting to produce on-shell mediator with mass " << mass_dp << " larger than pion mass " << MASS_PION_NEUTRAL << ".";
                    return -1;
                }
                DMGen = std::shared_ptr<DMGenerator>(new 
                        Two_Body_Decay_Gen(brpi0toVgamma(mass_dp,mass_dm,epsilon,alpha_D),
                        MASS_PION_NEUTRAL,"Pion",dark_photon,photon,0.0));
                Dist->set_mass(MASS_PION_NEUTRAL);
            }
            else if (prodchoice=="eta_decay"){
                if(mass_dp>MASS_ETA){
                    cerr << "Attempting to produce on-shell mediator with mass " << mass_dp << " larger than eta " << MASS_ETA << ".";
                    return -1;
                }
                DMGen = std::shared_ptr<DMGenerator>(new 
                        Two_Body_Decay_Gen(bretatoVgamma(mass_dp,mass_dm,epsilon,alpha_D),
                        MASS_ETA,"Eta",dark_photon,photon,0.0));
                Dist->set_mass(MASS_ETA);
            }
            if(prodchan.Meson_Per_Pi0()<=0){
			    //default for eta
                Vnum = DMGen->BranchingRatio()*num_pi0/30.0;
            }
		    else{
			    Vnum = DMGen->BranchingRatio()*num_pi0*
                    (prodchan.Meson_Per_Pi0());
            }
            DMGen->Set_Channel_Name(prodchoice);
            return true;
        }
    }
    else if(prodchoice=="pi0_decay"||prodchoice=="pi0_decay_baryonic"){
    	Dist->set_mass(MASS_PION_NEUTRAL);
    	if(prodchoice=="pi0_decay"){
			DMGen = std::shared_ptr<DMGenerator>(new pion_decay_gen(mass_dp, mass_dm, epsilon, alpha_D));
		}
		else{
			DMGen = std::shared_ptr<DMGenerator>(new pion_decay_gen_baryonic(mass_dp, mass_dm, epsilon, alpha_D));
		}
		Vnum=DMGen->BranchingRatio()*par.Protons_on_Target()*prodchan.Meson_Per_Pi0()*par.Pi0_per_POT();
		return true;
    }
    else if(prodchoice=="eta_decay"||prodchoice=="eta_decay_baryonic"){
		Dist->set_mass(MASS_ETA);
		
		if(prodchoice=="eta_decay"){ 
			DMGen = std::shared_ptr<DMGenerator>(new eta_decay_gen(mass_dp, mass_dm, epsilon, alpha_D));
		}
		else{
			DMGen = std::shared_ptr<DMGenerator>(new eta_decay_gen_baryonic(mass_dp, mass_dm, epsilon, alpha_D));
		}
		Vnum=DMGen->BranchingRatio()*par.Protons_on_Target()*prodchan.Meson_Per_Pi0()*par.Pi0_per_POT();
		return true;
	}
	else if(prodchoice=="omega_decay"||prodchoice=="rho_decay"||
            prodchoice=="phi_decay"||prodchoice=="omega_decay_baryonic"||prodchoice=="phi_decay_baryonic"){
		if(prodchoice=="omega_decay"){
			Dist->set_mass(momega);
			DMGen = std::shared_ptr<DMGenerator>(new omega_decay_gen(mass_dp, mass_dm, epsilon, alpha_D));
		}
		else if(prodchoice=="rho_decay"){
			Dist->set_mass(mrho);
			DMGen = std::shared_ptr<DMGenerator>(new rho_decay_gen(mass_dp, mass_dm, epsilon, alpha_D));
		}
		else if(prodchoice=="phi_decay"){
			Dist->set_mass(mphi);
			DMGen = std::shared_ptr<DMGenerator>(new phi_decay_gen(mass_dp, mass_dm, epsilon, alpha_D));
		}
		else if(prodchoice=="omega_decay_baryonic"){
			Dist->set_mass(momega);
			DMGen = std::shared_ptr<DMGenerator>(new omega_decay_gen_baryonic(mass_dp, mass_dm, epsilon, alpha_D));
		}
		else if(prodchoice=="phi_decay_baryonic"){
			Dist->set_mass(mphi);
			DMGen = std::shared_ptr<DMGenerator>(new phi_decay_gen_baryonic(mass_dp, mass_dm, epsilon, alpha_D));
		}
		Vnum = DMGen->BranchingRatio()*par.Protons_on_Target()*prodchan.Meson_Per_Pi0()*par.Pi0_per_POT();
		return true;
	}
	else if(prodchoice=="parton_production_baryonic"||prodchoice=="parton_production"){
		Dist->set_mass(mass_dp);
		DMGen = std::shared_ptr<DMGenerator>(new parton_V_gen(mass_dp, mass_dm, epsilon, alpha_D, prodchoice));
		return true;
	}
	else if(prodchoice=="V_decay"||prodchoice=="Brem_V"||
            prodchoice=="V_decay_baryonic"){
		Dist->set_mass(mass_dp);
		DMGen = std::shared_ptr<DMGenerator>(new V_decay_gen(mass_dp,mass_dm,epsilon,alpha_D,proddist));
		 if(proddist=="proton_brem"){
            if(prodchan.ptmax()<prodchan.ptmin() || prodchan.zmax() < 0 || prodchan.zmax()<prodchan.zmin() || prodchan.zmin() < 0){
                    std::cerr << "Invalid properties for production_distribution proton_brem." << endl;
                    return -1;
            }   
            std::shared_ptr<Proton_Brem_Distribution> pbd(new Proton_Brem_Distribution(par.Beam_Energy(),epsilon,mass_dp,prodchan.ptmax(),prodchan.zmax(),prodchan.zmin(),alpha_D,proddist,prodchan.ptmin()));
            Vnum = pbd->V_prod_rate()*par.Protons_on_Target();//V_decay should return a branching ratio for this. If it doesn't, I will have to change it in the future. prodchan.Num_per_pot()
            DMGen->Set_Channel_Name("Dark_Bremsstrahlung");
            Dist->Add_Dist(pbd);
        }
        return true;
	}
    else if(prodchoice=="piminus_capture"){
        Dist = std::shared_ptr<Distribution>(new DoNothingDist());
        if(proddist!="default"){
            cerr << "No production distribution supported for piminus_capture";
        }
		DMGen = std::shared_ptr<DMGenerator>(new piminus_capture_gen(mass_dp,mass_dm,epsilon,alpha_D));
        Vnum = DMGen->BranchingRatio()*par.Pi0_per_POT()*par.Protons_on_Target();
        return true;
    }
	return false;
}

bool Kinetic_Mixing::Prepare_Signal_Channel(Parameter& par){

	string sigchoice = par.Signal_Channel();

	double mdm = mass_dm;
	double mv = mass_dp;
	double alD = alpha_D;
	double kappa = epsilon;
	double max_scatter_energy = par.Max_Scatter_Energy();
	double min_scatter_energy = par.Min_Scatter_Energy();
	double EDMRES = par.EDM_RES();
	double max_dm_energy = par.Max_DM_Energy();

	std::shared_ptr<Scatter> SigGen;
	string sig_part_name;

	if(sigchoice=="NCE_electron"){
		SigGen = std::shared_ptr<Scatter>(new Electron_Scatter(mdm, mv, alD, kappa,max_scatter_energy,min_scatter_energy));
		sig_part_name = "DM";
	}
	else if(sigchoice=="NCE_nucleon"){
        SigGen = std::shared_ptr<Scatter>(new Nucleon_Scatter(mdm+EDMRES/100.0,max_dm_energy,EDMRES,mdm,mv,alD,kappa,max_scatter_energy,min_scatter_energy,"Kinetic_V",par.Coherent(),par.Get_Detector()));
        sig_part_name = "DM";
	}
    //Eventually these won't have to be separate channels, they'll just check the model!
	else if(sigchoice=="NCE_nucleon_baryonic"){
		SigGen = std::shared_ptr<Scatter>(new Nucleon_Scatter(mdm+EDMRES/100.0,max_dm_energy,EDMRES,mdm,mv,alD,kappa,max_scatter_energy,min_scatter_energy,"Baryonic_V",par.Coherent(),par.Get_Detector()));
		sig_part_name = "DM";
	}
	else if(sigchoice=="Pion_Inelastic"||sigchoice=="Pion_Inelastic_Charged"||sigchoice=="Inelastic_Delta_to_Gamma"){
		//I might need some checking for allowed energies.
        if(sigchoice=="Pion_Inelastic"){
		    SigGen = std::shared_ptr<Scatter>(new Pion_Inelastic(mdm+EDMRES/100.0,max_dm_energy,EDMRES,mdm,mv,alD,kappa,max_scatter_energy,min_scatter_energy));
        }
        else if(sigchoice=="Inelastic_Delta_to_Gamma"){
		    SigGen = std::shared_ptr<Scatter>(new Pion_Inelastic(mdm+EDMRES/100.0,max_dm_energy,EDMRES,mdm,mv,alD,kappa,max_scatter_energy,min_scatter_energy,1));
        }
        //This covers both charged and neutral pion production
        else if(sigchoice=="Pion_Inelastic_Charged"){
		    SigGen = std::shared_ptr<Scatter>(new Pion_Inelastic(mdm+EDMRES/100.0,max_dm_energy,EDMRES,mdm,mv,alD,kappa,max_scatter_energy,min_scatter_energy,2));
        }
		sig_part_name = "DM";
	}
	else if(sigchoice=="Inelastic_Nucleon_Scattering_Baryonic" || sigchoice=="Inelastic_Nucleon_Scattering"){
		if(par.Scatter_Dist_Filename()==""){
			cerr << "No scatter_dist_filename provided for " << sigchoice << endl;
			return -1;
		}
		if(par.Output_Mode()=="comprehensive"){
			cerr << sigchoice << " does not support comprehensive output. Use summary mode instead.\n";
			return -1;
		}
		SigGen = std::shared_ptr<Scatter>(new Inelastic_Nucleon_Scatter(mdm,mv,alD,kappa,sigchoice,par.Scatter_Dist_Filename()));
		min_scatter_energy=0;
		max_scatter_energy=1e9;
		sig_part_name = "DM";
	}
	else if(sigchoice=="Signal_Decay"){
		sig_part_name="Dark_Photon";
        double lifetime;
        vector<double> Branching_Ratios;
        vector<vector<Particle> > Final_States;
        cout << "Setting up signal decay!"  << endl;
        Particle electron(MASS_ELECTRON);
        electron.name = "Electron";
       
        Particle positron(MASS_ELECTRON);
        positron.name = "Positron";

        Particle DM(mdm);
        DM.name="DM";

        Particle muon(MASS_MUON);
        muon.name = "Muon";

        Particle anti_muon(MASS_MUON);
        anti_muon.name = "Anti-muon";

        Particle hadronic(0);
        hadronic.name = "Hadronic_Stuff";
        
        double GV = Gamma_V(mv,mdm,kappa,alD);
        lifetime=hbar/GV;
        cout << "Width: " << GV << " Lifetime " << lifetime << endl;
        double br = 0;
        if((br=Gamma_V_to_leptons(mv,kappa,MASS_ELECTRON)/GV)>0){
            Branching_Ratios.push_back(br);
            cout << "BR(V->e e) = " << br << endl;
            vector<Particle> vec;
            vec.push_back(electron);
            vec.push_back(positron);
            Final_States.push_back(vec);
        }
        
        if((br=Gamma_V_to_leptons(mv,kappa,MASS_MUON)/GV)>0){
            Branching_Ratios.push_back(br);
            cout << "BR(V->mu mu) = " << br << endl;
            vector<Particle> vec;
            vec.push_back(muon);
            vec.push_back(anti_muon);
            Final_States.push_back(vec);
        }

        if((br=Gamma_V_to_hadrons(mv,kappa)/GV)>0){
            Branching_Ratios.push_back(br);
            cout << "BR(V->hadronic) = " << br << endl;
            vector<Particle> vec;
            vec.push_back(hadronic);
            Final_States.push_back(vec);
        }
        
        if((br=GammaV_to_dm_dm(mv,mdm,kappa,alD)/GV)>0){
            Branching_Ratios.push_back(br);
            cout << "BR(V->DM DM) = " << br << endl;
            vector<Particle> vec;
            vec.push_back(DM);
            vec.push_back(DM);
            Final_States.push_back(vec);
        }
        SigGen = std::shared_ptr<Scatter>(new SignalDecay(lifetime, Branching_Ratios, Final_States));
	}
	else{
		return false;
	}
	Sig_list.push_back(SigGen);
	sig_part_vec.push_back(sig_part_name);
	return true;
}


/*bool Kinetic_Mixing::Set_Model_Parameters(Parameter& par){
    if(par.Query_Map(g_chi_key,gchi)&&(par.Query_Map(dark_matter_mass_key,mchi))&&par.Query_Map(dark_mediator_mass_key, ma)&&par.Query_Map(g_quark_key,gq)&&par.Query_Map(g_electron_key,gae)){
        return true;
    }
    return false;
}*/