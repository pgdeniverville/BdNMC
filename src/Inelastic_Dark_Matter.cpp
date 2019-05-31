#include "Model.h"
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
using std::cout;
using std::endl;
using std::bind;
using std::function;

using namespace std::placeholders;

//The ratio at which off-shell mode is activated.
const double off_shell_ratio=1.3;

bool Inelastic_Dark_Matter::Set_Model_Parameters(Parameter& par){
    if(!par.Query_Map("dark_photon_mass",mass_dp)){
        std::cerr << par.Model_Name() << " requires dark_photon_mass to be defined.\n";
        return false;
    }
    if(!par.Query_Map("epsilon", epsilon)){
        std::cerr << par.Model_Name() << " requires epsilon to be defined.\n";
        return false;
    }
    if(!par.Query_Map("dark_matter_mass",mass_dm1)){
        if(!par.Query_Map("dark_matter_1_mass", mass_dm1)||!par.Query_Map("dark_matter_2_mass", mass_dm2)){
            std::cerr << par.Model_Name() << " requires dark_matter_mass or dark_matter_1_mass and dark_matter_2_mass to be defined.\n";
            return false;
        }
    }
    else{
        mass_dm2 = mass_dm1;
    }
    if(mass_dm2<mass_dm1){
        cout << "dark_matter_2_mass is assumed to be heavier than dark_matter_1_mass, swapping masses.\n";
        double tmp = mass_dm2;
        mass_dm2 = mass_dm1;
        mass_dm1 = tmp;
    }
    if(!par.Query_Map("alpha_D", alpha_D)){
        std::cerr << par.Model_Name() << " requires alpha_D to be defined.\n";
        return false;
    }
    Evaluate_Widths();

    return true;
}

void Inelastic_Dark_Matter::Report(std::ostream& out){
    out << " " << mass_dp << " " << mass_dm1 << " " << mass_dm2 << " " << epsilon << " " << alpha_D << " ";
}

void Inelastic_Dark_Matter::Report_Model(){
    cout << "Mass Dark Photon = " << mass_dp << " GeV" << endl;
    cout << "Mass Dark Matter 1 = " << mass_dm1 << " GeV" << endl;
    cout << "Mass Dark Matter 2 = " << mass_dm2 << " GeV" << endl;
    cout << "epsilon = " << epsilon << endl;
    cout << "alpha_D = " << alpha_D << endl;
}

void Inelastic_Dark_Matter::Evaluate_Widths(){
    Awidth=Inelastic_DM::Gamma_A_to_dm1_dm2(mass_dp,mass_dm1,mass_dm2,alpha_D)+Gamma_V_to_visible(mass_dp,epsilon);

    cout << "Gamma(A->DM1+DM2) = " << Inelastic_DM::Gamma_A_to_dm1_dm2(mass_dp,mass_dm1,mass_dm2,alpha_D) << endl;

    dm2width=0;

    std::function<double(double, double)> d2widthe = std::bind(&Inelastic_Dark_Matter::Amplitude2_dm2_to_lepton_lepton_dm1,this,_1,_2,mass_dm2,MASS_ELECTRON,MASS_ELECTRON,mass_dm1);

    width_dm2_to_elec_elec_dm1=Three_Body_Decay_Space::integrate_decay_width(d2widthe,mass_dm2,MASS_ELECTRON,MASS_ELECTRON,mass_dm1);

    cout << "Gamma(dm2->e^+ e^- dm1) = " << width_dm2_to_elec_elec_dm1 << endl;

    std::function<double(double, double)> d2widthmu = std::bind(&Inelastic_Dark_Matter::Amplitude2_dm2_to_lepton_lepton_dm1,this,_1,_2,mass_dm2,MASS_MUON,MASS_MUON,mass_dm1);

    width_dm2_to_muon_muon_dm1=Three_Body_Decay_Space::integrate_decay_width(d2widthmu,mass_dm2,MASS_MUON,MASS_MUON,mass_dm1);

    cout << "Gamma(dm2->mu^+ mu^- dm1) = " << width_dm2_to_muon_muon_dm1 << endl;

    // std::function<double(double, double)> d2width_hadrons = std::bind(&Inelastic_Dark_Matter::Amplitude2_dm2_to_hadrons_dm1,this,_1,_2,mass_dm2,mass_dm1);
 
    // width_dm2_to_dm1_hadrons=Three_Body_Decay_Space::integrate_decay_width(d2width_hadrons,mass_dm2,MASS_MUON,MASS_MUON,mass_dm1);;
    
    dm2width = width_dm2_to_muon_muon_dm1+width_dm2_to_elec_elec_dm1;

}

//This should just be a function from branchingratios.cpp
double Inelastic_Dark_Matter::A_width(){
    return Awidth;
}

double Inelastic_Dark_Matter::dm2_width(){
    return dm2width;
}

//p1 -> dm1, p2 -> dm2, p3 -> photon
double Inelastic_Dark_Matter::meson_decay_amplitude2(double m12s, double m23s, double m0, double m1, double m2, double m3){
    return (4*(m12s*(pow(m12s,2) - m12s*(pow(m2,2) - 2*m23s) + 2*m23s*(-pow(m2,2) + m23s)) -2*(m12s + pow(m2,2))*(m12s - pow(m2,2) + m23s)*pow(m0,2) + 
       (m12s + pow(m2,2))*pow(m0,4) + 2*m1*m2*pow(m12s - pow(m0,2),2) - 
       pow(m1,2)*(m12s - pow(m0,2))*(m12s - 2*pow(m2,2) + 2*m23s - pow(m0,2)))*alpha_D*
     pow(alphaEM,2)*pow(epsilon,2))/
   (pow(PION_DECAY_CONSTANT,2)*pi*(pow(m12s - pow(mass_dp,2),2) + pow(mass_dp,2)*pow(A_width(),2)));
}

    //Includes spin averaging over initial states. 
    //p1, p2 -> lepton pair, p3 -> dm1
    //Assumes m1=m2
double Inelastic_Dark_Matter::Amplitude2_dm2_to_lepton_lepton_dm1(double m12s, double m23s, double m0, double m1, double m2, double m3){
    return (-64*(pow(m12s,2) + m12s*(2*m23s - pow(m3,2)) + pow(m0,2)*(-m12s - 2*m23s + 2*pow(m3,2)) + 
       2*m0*m3*(m12s + 2*pow(m1,2)) + 2*
        (pow(m23s,2) + pow(m1,4) - m23s*(pow(m3,2) + 2*pow(m1,2))))*pow(pi,2)*alpha_D*alphaEM*
     pow(epsilon,2))/(pow(m12s - pow(mass_dp,2),2) + pow(mass_dp,2)*pow(A_width(),2));
}

double Inelastic_Dark_Matter::Amplitude2_dm2_to_hadrons_dm1(double m12s, double m23s, double m0, double m3){
    return Amplitude2_dm2_to_lepton_lepton_dm1(m12s, m23s, m0, MASS_MUON, MASS_MUON, m3)*pow(RRATIO(mass_dp),2);
}

bool Inelastic_Dark_Matter::Prepare_Production_Channel(std::string prodchoice, std::string proddist, production_channel& prodchan, std::shared_ptr<DMGenerator>& DMGen, std::shared_ptr<Distribution>& Dist, double& Vnum, Parameter& par){

    Particle photon(0);
    photon.name="Photon";

    Particle dark_photon(mass_dp);
    dark_photon.name="Dark_Photon";

    Particle dm1(mass_dm1);
    dm1.name = "Dark_Matter_1";

    Particle dm2(mass_dm2);
    dm2.name = "Dark_Matter_2";

    string sig_choice = par.Signal_Channel();

    if(prodchoice=="pi0_decay"||prodchoice=="eta_decay"){
        double lifetime;
        
        Particle meson(0);
        
        std::function<double(double, double, double, double, double, double)> func;

        double on_shell_br=0;

        //THIS NEEDS TO BE UPDATED. 
        //NEED TO ADD HANDLING FOR:
        //1) Off-shell vs. on-shell decay
        //2) Rename the amplitudes to whatever I use in the file.
        //3) Chain decay of Dark_Matter_2
        //4) Any lingering references to Axion_Dark_Matter.
        //Can any of this be generalized?
        if(prodchoice=="pi0_decay"){
            if(MASS_PION_NEUTRAL<mass_dm1+mass_dm2){
                Vnum=0;
                DMGen = std::shared_ptr<DMGenerator>(new Do_Nothing_Gen());
                return true;
            }

            meson.name = "Pi0";
            meson.m = MASS_PION_NEUTRAL;
            meson.width=WIDTH_PION_NEUTRAL;
            Dist->set_mass(MASS_PION_NEUTRAL);
            lifetime=0;//Decay is effectively instantaneous.
            func = std::bind(&Inelastic_Dark_Matter::meson_decay_amplitude2,this,_1,_2,_3,_4,_5,_6);
            on_shell_br=brpi0toVgamma(mass_dp, 0, epsilon, alpha_D)*Inelastic_DM::Gamma_A_to_dm1_dm2(mass_dp,mass_dm1,mass_dm2,alpha_D)/A_width();
        }
        else if(prodchoice=="eta_decay"){
            if(MASS_ETA<mass_dm1+mass_dm2){
                Vnum=0;
                DMGen = std::shared_ptr<DMGenerator>(new Do_Nothing_Gen());
                return true;
            }
            meson.name = "eta";
            meson.m = MASS_ETA;
            meson.width=WIDTH_ETA;
            Dist->set_mass(MASS_ETA);
            lifetime=0;//Decay is effectively instantaneous.
            func = std::bind(&Inelastic_Dark_Matter::meson_decay_amplitude2,this,_1,_2,_3,_4,_5,_6);
            on_shell_br=bretatoVgamma(mass_dp, 0, epsilon, alpha_D)*Inelastic_DM::Gamma_A_to_dm1_dm2(mass_dp,mass_dm1,mass_dm2,alpha_D)/A_width();
        }
        else{
            std::cerr << "Invalid production channel" << prodchoice << "\n" << endl;
            throw -1;
        }
        std::function<double(double,double)> func_width = std::bind(func,_1,_2,meson.m,mass_dm1,mass_dm2,0);
 
        //Only calculates width to 2% accuracy.
        double Gamma_meson_3body=Three_Body_Decay_Space::integrate_decay_width(func_width,meson.m,mass_dm1,mass_dm2,0);
 
        cout << "Three Body Decay Width = " << Gamma_meson_3body << endl;

        std::shared_ptr<Three_Body_Decay_Gen> tmp_gen(new Three_Body_Decay_Gen(meson,dm1,dm2,photon,prodchoice,lifetime,Gamma_meson_3body,func,0));
        tmp_gen->d1=false;

        if(sig_choice!="DP_Signal_Decay"){
            //Switch to off-shell behavior
            if(tmp_gen->BranchingRatio()>off_shell_ratio*on_shell_br){
                tmp_gen->Burn_In(1000);
                DMGen = tmp_gen;

                if(sig_choice!="DM2_Signal_Decay"){
                    std::cerr << sig_choice << " not yet implemented\n";
                    throw -1;
                }

            }
            else{
                cout << "Turning on the on-shell decay of the Dark Photon with lifetime of " << hbar/A_width() << " seconds\n";
                std::shared_ptr<Two_Body_Decay_Gen> meson_decay_gen(new Two_Body_Decay_Gen(tmp_gen->BranchingRatio(),meson.m,meson.name,photon,dark_photon,hbar/meson.width));
 
                meson_decay_gen->d1=false;

                std::shared_ptr<Two_Body_Decay_Gen> tmp_gen2(new Two_Body_Decay_Gen(Inelastic_DM::Gamma_A_to_dm1_dm2(mass_dp, mass_dm1, mass_dm2,alpha_D)/A_width(),mass_dp,string("Dark_Photon"),dm1,dm2,hbar/A_width()));

                tmp_gen2->d1=false;//Don't care about DM1
                tmp_gen2->record_parent=false;
                meson_decay_gen->Toggle_Daughter_Decay(2,tmp_gen2);

                if(sig_choice!="DM2_Signal_Decay"){
                    std::cerr << sig_choice << " not yet implemented\n";
                    throw -1;
                }
                DMGen = meson_decay_gen;
            }       
        }

        Vnum=DMGen->BranchingRatio()*par.Protons_on_Target()*prodchan.Meson_Per_Pi0()*par.Pi0_per_POT();

        if(DMGen->Channel_Name()==""){
            DMGen->Set_Channel_Name(prodchoice);
        }
        
        return true;
    }
    else if(prodchoice=="proton_brem"){
        double part_mass;
        if(prodchoice=="proton_brem"){
            part_mass=mass_dp;
            if(sig_choice!="DP_Signal_Decay"&&A_width()>0){
                cout << "Turning on decay of the Dark Photon with lifetime of " << hbar/A_width() << " seconds\n";
                std::shared_ptr<Two_Body_Decay_Gen> invis_dec(new Two_Body_Decay_Gen(Inelastic_DM::Gamma_A_to_dm1_dm2(mass_dp,mass_dm1,mass_dm2,alpha_D)/A_width(),mass_dp,string("Dark_Photon"),dm1,dm2,hbar/A_width()));
                invis_dec->record_parent=false;
                DMGen = invis_dec;
            }
            else{
                std::shared_ptr<Do_Nothing_Gen> do_not(new Do_Nothing_Gen(string("Dark_Bremsstrahlung"),string("Dark_Photon")));
                DMGen = do_not;
            }
        }

        Dist->set_mass(part_mass);
        Vnum=par.Protons_on_Target()*prodchan.Num_per_pot();
        return true;
    }
    return false;
}

bool Inelastic_Dark_Matter::Prepare_Signal_Channel(Parameter& par){
    string sig_choice = par.Signal_Channel();
    if(sig_choice=="DM2_Signal_Decay"){
        sig_part_vec.push_back(string("Dark_Matter_2"));

        Particle dm2(mass_dm2);
        dm2.name = "Dark_Matter_2";
        dm2.width = dm2_width();
        Particle electron(MASS_ELECTRON);
        electron.name = "Decay_Electron";
        Particle positron(MASS_ELECTRON);
        positron.name = "Decay_Positron";
        Particle dm1(mass_dm1);
        dm1.name = "Decay_Dark_Matter_1";

        Particle muon(MASS_MUON);
        electron.name = "Decay_Muon";
        Particle antimuon(MASS_MUON);
        positron.name = "Decay_Antimuon";

        vector<std::shared_ptr<DMGenerator> > dec_vec;

        if(mass_dm2>mass_dm1+2*MASS_ELECTRON){
            std::function<double(double,double,double,double,double,double)> func = std::bind(&Inelastic_Dark_Matter::Amplitude2_dm2_to_lepton_lepton_dm1,this,_1,_2,_3,_4,_5,_6);

            //Set lifetime to zero, Signal_Decay handles the lifetime.
            std::shared_ptr<Three_Body_Decay_Gen> dp_3decay_gen(new Three_Body_Decay_Gen(dm2,electron,positron,dm1,string("Decay_Electron_Positron_DM1"),0,width_dm2_to_elec_elec_dm1,func));

            dp_3decay_gen->record_parent = false;            

            dec_vec.push_back(dp_3decay_gen);
        }
        
        if(mass_dm2>mass_dm1+2*MASS_MUON){

            std::function<double(double,double,double,double,double,double)> func2 = std::bind(&Inelastic_Dark_Matter::Amplitude2_dm2_to_lepton_lepton_dm1,this,_1,_2,_3,_4,_5,_6);

            std::shared_ptr<Three_Body_Decay_Gen> dp_3decay_gen_2(new Three_Body_Decay_Gen(dm2,muon,antimuon,dm1,string("Decay_Muon_Antimuon_DM1"),0,width_dm2_to_muon_muon_dm1,func2));

            dp_3decay_gen_2->record_parent = false;

            dec_vec.push_back(dp_3decay_gen_2);

        }

        if(dec_vec.size()==0){
            cout << "No valid decay channels for DM2!\n";
            return false;
        }

        double lifetime = hbar/dm2_width();

        std::shared_ptr<SignalDecay_2> sig_dec(new SignalDecay_2(lifetime, dec_vec));
        Sig_list.push_back(sig_dec);
        return true;
    }
    return false;
}
