#include "Model.h"
#include "SignalDecay.h"
#include "constants.h"
#include "Kinematics.h"
#include "decay.h" //We'll see if this is needed.
#include "branchingratios.h"
#include "Distribution.h"

#include <cmath>
#include <iomanip>

using std::string;
using std::vector;
using std::cout;
using std::endl;
using std::bind;
using std::function;

const double Z = 35.45; const double B = 0.308;
const double C1 = 42.53; const double C2 = 33.34;
const double s0 = pow(5.38,2); const double s1 = 1;
const double eta1 = 0.458; const double eta2 = 0.545;

const double g_SNN = 1.2e-3; //Scalar Effective Coupling to Nucleon

using namespace std::placeholders;

double sigma_pp(double s){
	return Z + B*pow(log(s/s0),2) + C1*pow(s1/s,eta1) - C2*pow(s1/s,eta2);
}

//This is actually just the splitting probability, which is all we need if we assume that all protons interact with the target.

double Scalar_Mediator::differential_splitting_probability(double z, double pT){
	double num = z*(pow(MASS_PROTON*(2-z),2)+pT*pT);
	double den = pow(pow(MASS_PROTON*z,2) + m_S*m_S*(1-z)+pT*pT,2);
	return pow(g_SNN*epsilon_q,2)/(8*pi*pi)*(2*pT)*num/den;
}

double Scalar_Mediator::sigma_brem_differential(double z, double pT){
	double E_s = sqrt(pow(pT,2) + pow(P_beam*z,2) + pow(m_S,2));
	double s_prime = 2*MASS_PROTON*(E_beam - E_s + MASS_PROTON);
	//cout << differential_splitting_probability(z,pT) << endl; 
	return differential_splitting_probability(z,pT)*sigma_pp(s_prime);
}

double Scalar_Mediator::S_width(){
	return Dark_Scalar::S_to_2leptons(epsilon_l,m_S,MASS_ELECTRON)+Dark_Scalar::S_to_2leptons(epsilon_l,m_S,MASS_MUON)+Dark_Scalar::S_to_2leptons(epsilon_l,m_S,MASS_TAU)+Dark_Scalar::S_to_2photon(epsilon_l,epsilon_q,epsilon_W,m_S);
}

double Scalar_Mediator::S_to_2lepton_width(double m_lepton){
	return Dark_Scalar::S_to_2leptons(epsilon_l,m_S,m_lepton);
}

bool Scalar_Mediator::Set_Model_Parameters(Parameter& par){
	if(!par.Query_Map("dark_scalar_mass",m_S)){
		std::cerr << par.Model_Name() << " requires dark_scalar_mass to be defined.\n";
        return false;
	}
	//Need to know what to call the coupling.
	if(!par.Query_Map("epsilon_l",epsilon_l)){
		std::cerr << par.Model_Name() << " requires epsilon_l to be defined.\n";
        return false;
	}
	if(!par.Query_Map("epsilon_q",epsilon_q)){
		std::cerr << par.Model_Name() << " requires epsilon_q to be defined.\n";
        return false;
	}	
	if(!par.Query_Map("epsilon_W",epsilon_W)){
		epsilon_W=epsilon_q;
		std::cout << "Setting epsilon_W=epsilon_q.\n";
	}	
	E_beam = par.Beam_Energy();
	P_beam = sqrt(E_beam*E_beam-pow(MASS_PROTON,2));
	return true;
}

void Scalar_Mediator::Report(std::ostream& out){
    out << m_S << " " << epsilon_q << " " << epsilon_l << " " << epsilon_W << " " << S_width() << " ";
}

void Scalar_Mediator::Report_Model(){
    cout << "Mass Dark Scalar = " << m_S << " GeV" << endl;
    cout << "epsilon_q = " << epsilon_q << endl;
    cout << "epsilon_l = " << epsilon_l << endl;
    cout << "epsilon_w = " << epsilon_W << endl;
}

bool Scalar_Mediator::Prepare_Production_Channel(std::string prodchoice, std::string proddist, production_channel& prodchan, std::shared_ptr<DMGenerator>& DMGen, std::shared_ptr<Distribution>& Dist, double& Vnum, Parameter& par){

	Particle Scalar(m_S);
	Scalar.name="Dark_Scalar";

    string sig_choice = par.Signal_Channel();

	if(prodchoice=="proton_brem"){
		if(prodchan.ptmax()<prodchan.ptmin() || prodchan.zmax() < 0 || prodchan.zmax()<prodchan.zmin()){
                    std::cerr << "Invalid properties for production_channel proton_brem." << endl;
                    throw -1;
            }

        //std::shared_ptr<BeamDistribution> BD = std::shared_ptr<BeamDistribution>(new BeamDistribution(par.Beam_Energy(), MASS_PROTON));

        //Dist = BD;

        //May update to account for target length.
        double production_fraction = 1;
        //Beam_Energy is assumed to be kinetic.
        double sig_tot = sigma_pp(2*MASS_PROTON*(par.Beam_Energy()+2*MASS_PROTON));

        std::function<double(double, double)> sigma_brem = std::bind(&Scalar_Mediator::sigma_brem_differential,this,_1,_2);

        //Write a Generic distribution for use here. Production distribution is a function pointer to some pT,z distribution, production_fraction is the likelihood that a proton interacts with the target. Maybe default ot 1 if target details are missing. 
        std::shared_ptr<Generic_Bremsstrahlung> pbd(new Generic_Bremsstrahlung(Scalar, par.Beam_Energy(), prodchan.ptmin(), prodchan.ptmax(),prodchan.zmin(),prodchan.zmax(), "Dark_Bremsstrahlung", sigma_brem));

        Vnum = pbd->Prod_Rate()/sig_tot*par.Protons_on_Target()*production_fraction;

        DMGen = pbd;

    	return true; 
	}

	return false;
}

bool Scalar_Mediator::Prepare_Signal_Channel(Parameter& par){
	string sig_choice = par.Signal_Channel();

	if(sig_choice=="Scalar_Signal_Decay"){

		sig_part_vec.push_back(string("Dark_Scalar"));

		Particle Scalar(m_S);
		Scalar.name = "Dark_Scalar";

		sig_part_vec.push_back(Scalar.name);

        vector<std::shared_ptr<DMGenerator> > dec_vec;

        //Need to add in the Branching Ratio
        if(m_S > 2*MASS_ELECTRON){
        	Particle electron(MASS_ELECTRON);
			electron.name = "Decay_Electron";
			Particle positron(MASS_ELECTRON);
			positron.name = "Decay_Positron";
        	std::shared_ptr<Two_Body_Decay_Gen> s_2decay_gen = std::shared_ptr<Two_Body_Decay_Gen>(new Two_Body_Decay_Gen(S_to_2lepton_width(MASS_ELECTRON)/S_width(), m_S,Scalar.name, electron, positron));
        	s_2decay_gen->record_parent = false;
        	dec_vec.push_back(s_2decay_gen);
        }
        if(m_S > 2*MASS_MUON){
        	Particle muon(MASS_MUON);
			muon.name = "Decay_Muon";
			Particle antimuon(MASS_MUON);
			antimuon.name = "Decay_Antimuon";
        	std::shared_ptr<Two_Body_Decay_Gen> s_2decay_gen = std::shared_ptr<Two_Body_Decay_Gen>(new Two_Body_Decay_Gen(S_to_2lepton_width(MASS_ELECTRON)/S_width(), m_S,Scalar.name, muon, antimuon));
        	s_2decay_gen->record_parent = false;
        	dec_vec.push_back(s_2decay_gen);	
        }

        double lifetime = hbar/S_width();

        std::shared_ptr<SignalDecay_2> sig_dec(new SignalDecay_2(lifetime, dec_vec));
        Sig_list.push_back(sig_dec);

        return true;
	}

	return false;
}