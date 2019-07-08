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
#include "Proton_Brem_Distribution.h"

using std::string;
using std::vector;
using std::cout;
using std::endl;
using std::bind;
using std::function;

using namespace std::placeholders;

//The ratio at which off-shell mode is activated.
const double off_shell_ratio=1.3;

const string GM_form_factor_filename = "data/delta_production_form_factor.dat";

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
    out << mass_dp << " " << mass_dm1 << " " << mass_dm2 << " " << epsilon << " " << alpha_D << " " << dm2_width() << " " << A_width() << " ";
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

//This amplitude is averaged over initial spins (1/4)
double Inelastic_Dark_Matter::dm_e_to_dm_e_amp(double s, double t, double mass_dm_in, double mass_dm_out, double mR){
    return 1/pow(pow(mass_dp,2)-t,2)*32*pi*pi*alpha_D*alphaEM*pow(epsilon,2)*(2*(pow(mass_dm_in*mass_dm_out + pow(mR,2), 2) - (pow(mass_dm_in,2) + pow(mass_dm_out,2) + 2*pow(mR,2))*s + pow(s,2)) - (pow(mass_dm_in - mass_dm_out,2) - 2*s)*t + pow(t,2));
}

double Inelastic_Dark_Matter::dsigma_dm_e_to_dm_e(double E1lab, double E4, double mass_dm_in, double mass_dm_out, double mR){
    double s = two_to_two_scattering::s_lab(E1lab, mass_dm_in, mR);
    if(s<pow(mR+mass_dm_out,2)){
        return 0;
    }
    return 2*mR*1/64.0/pi/pow(mR,2)/(pow(E1lab,2)-pow(mass_dm_in,2))*dm_e_to_dm_e_amp(s,two_to_two_scattering::t_lab(E4,mR,mR),mass_dm_in,mass_dm_out,mR);
}

double Inelastic_Dark_Matter::dsigma_dm_e_to_dm_e_2(double E1lab, double t, double mass_dm_in, double mass_dm_out, double mR){
    double s = two_to_two_scattering::s_lab(E1lab, mass_dm_in, mR);
    return 1/64.0/pi/pow(mR,2)/(pow(E1lab,2)-pow(mass_dm_in,2))*dm_e_to_dm_e_amp(s,t,mass_dm_in,mass_dm_out,mR);
}

//m1, m2 can be either dm1 or dm2, whichever is inbound.
double Inelastic_Dark_Matter::amp_dm_N_to_dm_Delta(double s, double t, double m1, double m2, double mN, double mD){
    //1/4 is to average over initial spins.
/*    return 1/4.0*(8*pow(mD + mN,2)*pow(pi,2)*(2*(pow(mD,2) + pow(mN,2)) - t)*t*
     (pow(m1,4)*(-16*pow(mass_dp,2)*pow(mD,2) + pow(mD,4) - 2*pow(mD,2)*(pow(mN,2) - 5*t) + 
          pow(pow(mN,2) - t,2)) + 2*m1*m2*
        (36*pow(mass_dp,4)*pow(mD,2) - 2*pow(mass_dp,2)*
           (pow(mD,4) - 2*pow(mD,2)*(pow(mN,2) - 5*t) + pow(pow(mN,2) - t,2)) + 
          (pow(mD,4) - 2*pow(mD,2)*(pow(mN,2) - 5*t) + pow(pow(mN,2) - t,2))*t) + 
       pow(m2,4)*(4*pow(mass_dp,4) + pow(mD,4) - 2*pow(mD,2)*(pow(mN,2) - 5*t) + 
          pow(pow(mN,2) - t,2) - 4*pow(mass_dp,2)*(5*pow(mD,2) - pow(mN,2) + t)) + 
       4*pow(mass_dp,4)*(s*(-pow(mN,2) + s + t) + pow(mD,2)*(pow(mN,2) - s + 4*t)) - 
       pow(m1,2)*(20*pow(mass_dp,4)*pow(mD,2) + 
          (pow(mD,4) - 2*pow(mD,2)*(pow(mN,2) - 5*t) + pow(pow(mN,2) - t,2))*t + 
          4*pow(mass_dp,2)*(-pow(mD,4) + pow(mD,2)*(pow(mN,2) + s - 5*t) + s*(-pow(mN,2) + t))) + 
       pow(m2,2)*(-((pow(mD,4) - 2*pow(mD,2)*(pow(mN,2) - 5*t) + pow(pow(mN,2) - t,2))*
             t) - 4*pow(mass_dp,4)*(4*pow(mD,2) - pow(mN,2) + 2*s + t) + 
          pow(m1,2)*(-2*(pow(mD,4) - 2*pow(mD,2)*(pow(mN,2) - 5*t) + 
                pow(pow(mN,2) - t,2)) + 4*pow(mass_dp,2)*(9*pow(mD,2) - pow(mN,2) + t)) + 
          4*pow(mass_dp,2)*((pow(mN,2) - t)*(pow(mN,2) - s - t) + 
             pow(mD,2)*(-pow(mN,2) + s + 5*t))))*alpha_D*alphaEM*pow(epsilon,2)*pow(GM(-t),2))/
   (pow(mass_dp,4)*pow(mD,2)*pow(mN,2)*pow(pow(mass_dp,2) - t,2)*(pow(mD + mN,2) - t));
*/
    return (8.0*alpha_D*alphaEM*pow(epsilon,2)*pow(mD + mN,2)*pow(pi,2)*
     (2*pow(m1,4)*pow(mD,2) + 2*pow(m2,4)*pow(mN,2) - 
       pow(m2,2)*(pow(mD,2) - pow(mN,2) - t)*(pow(mD,2) + pow(mN,2) - 2*s - t) + 
       pow(m1,2)*(-2*pow(m2,2)*(pow(mD,2) + pow(mN,2) - t) + 
          (pow(mD,2) + pow(mN,2) - 2*s - t)*(pow(mD,2) - pow(mN,2) + t)) + 
       2*m1*m2*(pow(mD,4) + pow(pow(mN,2) - t,2) - 2*pow(mD,2)*(pow(mN,2) + t)) + 
       t*(pow(mD,4) + pow(mN,4) + 2*pow(s,2) + 2*s*t + pow(t,2) - 2*pow(mD,2)*(s + t) - 2*pow(mN,2)*(s + t)))*
     pow(GM(-t),2))/(pow(mN,2)*pow(pow(mass_dp,2) - t,2)*(pow(mD + mN,2) - t));
}

double Inelastic_Dark_Matter::GM(double q2){
    if(q2>MAX_Q2){
        if(!MAX_Q2_WARNING){
            MAX_Q2_WARNING=true;
            std::cerr << "First request of out of bounds q^2 for pion inelastic form factor.\n";
            std::cerr << "q^2 tried =" << q2 << "\nMAX_Q2 = " <<  MAX_Q2 << endl;
        }
        return 0.0;
    }
    return GM_form_factor->Interpolate(q2);
}

double Inelastic_Dark_Matter::dsigma_dm_N_to_dm_Delta(double E1lab, double E4, double mass_dm_in, double mass_dm_out, double mN, double mDelta){
    double s = two_to_two_scattering::s_lab(E1lab, mass_dm_in, mN);
    if(s<pow(mDelta+mass_dm_out,2)){
        return 0;
    }
    return -2*mN*1/64.0/pi/pow(mN,2)/(pow(E1lab,2)-pow(mass_dm_in,2))*amp_dm_N_to_dm_Delta(s,two_to_two_scattering::t_lab(E4,mN,mDelta),mass_dm_in,mass_dm_out,mN,mDelta);
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
                tmp_gen->set_Off_Shell(true);
                DMGen = tmp_gen;

                if(sig_choice!="DM2_Signal_Decay"&&dm2width>0){
                    std::cerr << "Need to add DM2 decay! Check Model file.\n";
                    throw -1;
                }
            }
            else{
                cout << "Turning on the on-shell decay of the Dark Photon with lifetime of " << hbar/A_width() << " seconds\n";
                std::shared_ptr<Two_Body_Decay_Gen> meson_decay_gen(new Two_Body_Decay_Gen(tmp_gen->BranchingRatio(),meson.m,meson.name,photon,dark_photon,hbar/meson.width));
 
                meson_decay_gen->d1=false;
                meson_decay_gen->d2=false;

                std::shared_ptr<Two_Body_Decay_Gen> tmp_gen2(new Two_Body_Decay_Gen(Inelastic_DM::Gamma_A_to_dm1_dm2(mass_dp, mass_dm1, mass_dm2,alpha_D)/A_width(),mass_dp,string("Dark_Photon"),dm1,dm2,hbar/A_width()));

                tmp_gen2->d1=false;//Don't care about DM1
                tmp_gen2->record_parent=true;
                meson_decay_gen->Toggle_Daughter_Decay(2,tmp_gen2);

                if(sig_choice!="DM2_Signal_Decay"&&dm2width>0){
                    std::cerr << "Need to add DM2 decay! Check Model file.\n";
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
    }//Need to add proton beam handling eventually.
    else if(prodchoice=="V_decay"){
//        double part_mass;
//        part_mass=mass_dp;        
        if(sig_choice!="DP_Signal_Decay"&&A_width()>0){
            cout << "Turning on decay of the Dark Photon with lifetime of " << hbar/A_width() << " seconds\n";
            std::shared_ptr<Two_Body_Decay_Gen> invis_dec(new Two_Body_Decay_Gen(Inelastic_DM::Gamma_A_to_dm1_dm2(mass_dp,mass_dm1,mass_dm2,alpha_D)/A_width(),mass_dp,string("Dark_Photon"),dm1,dm2,hbar/A_width()));
            invis_dec->record_parent=true   ;
            DMGen = invis_dec;
            if(sig_choice=="DM2_Signal_Decay"){
                invis_dec->d1=false;
            }
            //Option for decaying dm2 needs to be included for electron scatter.
        }
        else{
            std::shared_ptr<Do_Nothing_Gen> do_not(new Do_Nothing_Gen(string("Dark_Bremsstrahlung"),string("Dark_Photon")));
            DMGen = do_not;
        }
        if(proddist=="proton_brem"){
            if(prodchan.ptmax()<prodchan.ptmin() || prodchan.zmax() < 0 || prodchan.zmax()<prodchan.zmin() || prodchan.zmin() < 0){
                    std::cerr << "Invalid properties for production_distribution proton_brem." << endl;
                    return -1;
            }   
            std::shared_ptr<Proton_Brem_Distribution> pbd(new Proton_Brem_Distribution(par.Beam_Energy(),epsilon,mass_dp,prodchan.ptmax(),prodchan.zmax(),prodchan.zmin(),alpha_D,proddist,prodchan.ptmin()));
            Vnum = pbd->V_prod_rate()*par.Protons_on_Target();//V_decay should return a branching ratio for this. If it doesn't, I will have to change it in the future. prodchan.Num_per_pot()
            DMGen->Set_Channel_Name("Dark_Bremsstrahlung");
            Dist->Add_Dist(pbd);
            return true;
        }

        //Dist->set_mass(part_mass);
        Vnum=par.Protons_on_Target()*prodchan.Num_per_pot();
        return true;
    }
    else if(prodchoice==""){
    
    }
    return false;
}

bool Inelastic_Dark_Matter::Prepare_Signal_Channel(Parameter& par){
    string sig_choice = par.Signal_Channel();
    if(sig_choice=="Electron_Scatter" || sig_choice=="NCE_electron" or sig_choice=="Pion_Inelastic" or sig_choice == "Pion_Inelastic_Charged" or sig_choice=="Inelastic_Delta_to_Gamma"){
 
        // cout << "Beginning diagnostic\n";

        // double E1lab=1;
        // double s = pow(mass_dm1,2)+pow(MASS_ELECTRON,2)+2*E1lab*MASS_ELECTRON;

        // cout << "s=" << s << endl; 

        // cout << dm_e_to_dm_e_amp(s,-1e-6,mass_dm1,mass_dm2,MASS_ELECTRON) << endl;

        sig_part_vec.push_back(string("Dark_Matter_1"));
        sig_part_vec.push_back(string("Dark_Matter_2"));

        std::shared_ptr<Two_to_Two_Scatter> ttts(new Two_to_Two_Scatter());

        ttts->set_energy_limits(par.Min_Scatter_Energy(),par.Max_Scatter_Energy());
        ttts->set_angle_limits(par.Max_Angle(),par.Min_Angle());

        if(sig_choice=="Pion_Inelastic" or sig_choice == "Pion_Inelastic_Charged" or sig_choice=="Inelastic_Delta_to_Gamma"){

            GM_form_factor = std::shared_ptr<Linear_Interpolation>();
            MAX_Q2=load_form_factor(GM_form_factor_filename, GM_form_factor);

            double PNtot =(par.Get_Detector())->PNtot();
            double NNtot =(par.Get_Detector())->NNtot();

            Particle dm1_r(mass_dm1);
            dm1_r.name = "Recoil_Dark_Matter_1";
            Particle dm2_r(mass_dm2);
            dm2_r.name = "Recoil_Dark_Matter_2";
            Particle Delta(MASS_DELTA);
            Delta.name="Recoil_Delta";

            Particle pi0(MASS_PION_NEUTRAL);
            pi0.name = "pi0";
            Particle pion(MASS_PION_CHARGED);
            pi0.name = "pion";
            Particle photon(0);
            photon.name = "photon";
            //This could instead be a neutron, but I really don't want to handle Delta^+ vs Delta^0 vs Delta^++.
            Particle decay_nucleon(MASS_PROTON);
            decay_nucleon.name = "decay_nucleon";

            std::function<double(double,double)> dm1_p_to_dm2 = std::bind(&Inelastic_Dark_Matter::dsigma_dm_N_to_dm_Delta,this,_1,_2,mass_dm1,mass_dm2,MASS_PROTON,MASS_DELTA);

            ttts->Build_Channel(dm2_r, Delta, mass_dm1, MASS_PROTON, dm1_p_to_dm2, PNtot, "Dark_Matter_1", par.Max_DM_Energy(),par.EDM_RES());
            
            std::function<double(double,double)> dm2_p_to_dm1 = std::bind(&Inelastic_Dark_Matter::dsigma_dm_N_to_dm_Delta,this,_1,_2,mass_dm2,mass_dm1,MASS_PROTON,MASS_DELTA);

            ttts->Build_Channel(dm1_r, Delta, mass_dm2, MASS_PROTON, dm2_p_to_dm1, PNtot, "Dark_Matter_2", par.Max_DM_Energy(),par.EDM_RES());

            std::function<double(double,double)> dm1_n_to_dm2 = std::bind(&Inelastic_Dark_Matter::dsigma_dm_N_to_dm_Delta,this,_1,_2,mass_dm1,mass_dm2,MASS_NEUTRON,MASS_DELTA);

            ttts->Build_Channel(dm2_r, Delta, mass_dm1, MASS_NEUTRON, dm1_n_to_dm2, NNtot, "Dark_Matter_1", par.Max_DM_Energy(),par.EDM_RES());
            
            std::function<double(double,double)> dm2_n_to_dm1 = std::bind(&Inelastic_Dark_Matter::dsigma_dm_N_to_dm_Delta,this,_1,_2,mass_dm2,mass_dm1,MASS_NEUTRON,MASS_DELTA);

            ttts->Build_Channel(dm1_r, Delta, mass_dm2, MASS_NEUTRON, dm2_n_to_dm1, NNtot, "Dark_Matter_2", par.Max_DM_Energy(),par.EDM_RES());

            if(sig_choice=="Pion_Inelastic"){
                std::shared_ptr<Two_Body_Decay_Gen> pi0dec(new Two_Body_Decay_Gen(Delta_to_pi0, MASS_DELTA, "Recoil_Delta", pi0, decay_nucleon));
                ttts->add_decay("Recoil_Delta", pi0dec);
            }
            if(sig_choice=="Pion_Inelastic_Charged"){
                std::shared_ptr<Two_Body_Decay_Gen> piondec(new Two_Body_Decay_Gen(Delta_to_pion_charged, MASS_DELTA, "Recoil_Delta", pion, decay_nucleon));
                ttts->add_decay("Recoil_Delta", piondec);
            }
            if(sig_choice=="Inelastic_Delta_to_Gamma"){
                std::shared_ptr<Two_Body_Decay_Gen> photondec(new Two_Body_Decay_Gen(Delta_to_gamma, MASS_DELTA, "Recoil_Delta", photon, decay_nucleon));
                ttts->add_decay("Recoil_Delta", photondec);
            }
            Sig_list.push_back(ttts);
            return true;
        }
        else if(sig_choice=="Electron_Scatter" || sig_choice=="NCE_electron"){
            Particle electron(MASS_ELECTRON);
            electron.name="Electron";
            Particle dm1_r(mass_dm1);
            dm1_r.name = "Recoil_Dark_Matter_1";
            Particle dm2_r(mass_dm2);
            dm2_r.name = "Recoil_Dark_Matter_2";
            
            //Need to add decay of DM2 later on.
            /*if(dm2_width()>0){
                Particle dm1_f(mass_dm1);
                axion_f.name="Decay_dm1";
                Particle photon_f(0);
                photon_f.name = "Decay_Photon";
                std::shared_ptr<Two_Body_Decay_Gen> invis_dec(new Two_Body_Decay_Gen(A_to_a_gamma_width(mass_axion, mass_dp)/A_width(mass_axion, mass_dp),mass_dp,string("Recoil_Dark_Photon"),photon_f,axion_f,hbar/A_width(mass_axion, mass_dp)));
                invis_dec->record_parent=false;
                ttts->add_decay(string("Recoil_Dark_Photon"),invis_dec);
            }
            */
            //dm1+e->dm2+e
            std::shared_ptr<Linear_Interpolation> dm1_cross;
            std::shared_ptr<Linear_Interpolation> dm1_cross_max;
            

            std::function<double(double,double)> f_dm1_to_dm2 = std::bind(&Inelastic_Dark_Matter::dsigma_dm_e_to_dm_e,this,_1,_2,mass_dm1,mass_dm2,MASS_ELECTRON);
            
            std::function<double(double)> ER_min_dm1 = std::bind(&Two_to_Two_Scatter::scatmin,*ttts,_1,mass_dm1,MASS_ELECTRON,mass_dm2,MASS_ELECTRON);
            std::function<double(double)> ER_max_dm1 = std::bind(&Two_to_Two_Scatter::scatmax,*ttts,_1,mass_dm1,MASS_ELECTRON,mass_dm2,MASS_ELECTRON);

            Prepare_Cross_Section(f_dm1_to_dm2, ER_min_dm1, ER_max_dm1, dm1_cross,dm1_cross_max, mass_dm1, par.Max_DM_Energy(),par.EDM_RES());
            
            function<double(double)> dm1_cross_func = bind(&Linear_Interpolation::Interpolate,dm1_cross,_1);
            function<double(double)> dm1_cross_max_func = bind(&Linear_Interpolation::Interpolate,dm1_cross_max,_1);

            //dm2+e->dm1+e
            std::shared_ptr<Linear_Interpolation> dm2_cross;
            std::shared_ptr<Linear_Interpolation> dm2_cross_max;

            std::function<double(double,double)> f_dm2_to_dm1 = std::bind(&Inelastic_Dark_Matter::dsigma_dm_e_to_dm_e,this,_1,_2,mass_dm2,mass_dm1,MASS_ELECTRON);

            std::function<double(double)> ER_min_dm2 = std::bind(&Two_to_Two_Scatter::scatmin,*ttts,_1,mass_dm2,MASS_ELECTRON,mass_dm1,MASS_ELECTRON);
            std::function<double(double)> ER_max_dm2 = std::bind(&Two_to_Two_Scatter::scatmax,*ttts,_1,mass_dm2,MASS_ELECTRON,mass_dm1,MASS_ELECTRON);

            Prepare_Cross_Section(f_dm2_to_dm1, ER_min_dm2, ER_max_dm2,dm2_cross,dm2_cross_max, mass_dm2, par.Max_DM_Energy(),par.EDM_RES());

            function<double(double)> dm2_cross_func = bind(&Linear_Interpolation::Interpolate,*dm2_cross,_1);
            function<double(double)> dm2_cross_max_func = bind(&Linear_Interpolation::Interpolate,*dm2_cross_max,_1);


            double ENtot =(par.Get_Detector())->ENtot();
            ttts->add_channel(dm1_r, electron, MASS_ELECTRON, dm2_cross_func, dm2_cross_max_func, f_dm2_to_dm1, ENtot,"Dark_Matter_2");
            ttts->add_channel(dm2_r, electron, MASS_ELECTRON, dm1_cross_func, dm1_cross_max_func, f_dm1_to_dm2, ENtot,"Dark_Matter_1");

        }

        Sig_list.push_back(ttts);
        return true;
    }
    else if(sig_choice=="DM2_Signal_Decay"){
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
        muon.name = "Decay_Muon";
        Particle antimuon(MASS_MUON);
        antimuon.name = "Decay_Antimuon";

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
