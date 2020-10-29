#include "constants.h"
#include "branchingratios.h"
#include "Scatter.h"
#include <cmath>
#include "decay.h"
#include "DMgenerator.h"
#include "Model.h"
#include "Kinematics.h"
#include "SignalDecay.h"

using std::string;
using std::vector;
using std::cout;
using std::endl;
using std::bind;
using std::function;

using namespace std::placeholders;

//m12s=m12^2=(p1+p2)^2, s=m23^2=(p2+p3)^2
//2 -> axion
//3 -> Dark Photon
double Axion_Dark_Photon::meson_decay_amplitude2(double m12s, double m23s, double mmeson, double mgamma, double m2, double m3){
    return (pow(Gagpg,2)*(2*pow(m12s,2)*pow(m23s,2) - 2*pow(m2,2)*(m23s - m3*mmeson)*(m23s + m3*mmeson)*(m23s - pow(mmeson,2)) + pow(m23s - pow(m3,2),2)*
pow(m23s - pow(mmeson,2),2) + pow(m2,4)*(pow(m23s,2) + pow(mmeson,4)) + 2*m12s*m23s*((m23s - pow(m3,2))*(m23s - pow(mmeson,2))\
    - pow(m2,2)*(m23s + pow(mmeson,2))))*pow(alphaEM,2))/(8.*pow(PION_DECAY_CONSTANT,2)*pow(m23s,2)*pow(pi,2));
}

/*{
    return 1.0/pow(2*pi*PION_DECAY_CONSTANT*s,2)*pow(alphaEM*Gagpg,2)*(-2*pow(ma,2)*(m12s*pow(mmeson,2) + 
      pow(pow(mmeson,2) - s,2))*s + 
   s*(-((pow(mA,2) - s)*
         pow(pow(mmeson,2) - s,2)) + 
      pow(m12s,2)*s) + 
   pow(ma,4)*(2*pow(mmeson,4) - 
      2*pow(mmeson,2)*s + pow(s,2)));
    //(s*(m12s*m12s*s-pow(mA*mA-s,2))-2*ma*ma*s*(m12s*mpi0*mpi0+pow(mpi0*mpi0-s,2))+pow(ma,4)*(2*pow(mpi0,4)-2*mpi0*mpi0*s+s*s));
}*/

double Axion_Dark_Photon::eta_decay_amplitude2_c(double m12s, double s, double m0, double m1, double m2, double m3){
    //WIDTH_ETA/WIDTH_PION_NEUTRAL*bretato2gamma/brpi0to2gamma*
    return meson_decay_amplitude2(m12s,s,m0,m1,m2,m3);
}

double Axion_Dark_Photon::eta_decay_width_2(double m12, double CosTheta, double m0, double m1, double m2, double m3){
    //WIDTH_ETA/WIDTH_PION_NEUTRAL*bretato2gamma/brpi0to2gamma*
    return meson_decay_amplitude2(m12*m12,Three_Body_Decay_Space::Cos_Theta_to_m23s(m12*m12, CosTheta, m0, m1, m2, m3),m0,m1,m2,m3);
}

double Axion_Dark_Photon::pi0_decay_amplitude2(double m12s, double m23s, double mgamma, double ma, double mA){
    return meson_decay_amplitude2(m12s, m23s, MASS_PION_NEUTRAL, mgamma, ma, mA);
}

//Assumes mgamma=0;
/*
double Axion_Dark_Photon::pi0_decay_amplitude2(double m12s, double s, double mgamma, double ma, double mA){
    return 1.0/pow(2*pi*PION_DECAY_CONSTANT*s,2)*pow(alphaEM*Gagpg,2)*(-2*pow(ma,2)*(m12s*pow(mpi0,2) + 
      pow(pow(mpi0,2) - s,2))*s + 
   s*(-((pow(mA,2) - s)*
         pow(pow(mpi0,2) - s,2)) + 
      pow(m12s,2)*s) + 
   pow(ma,4)*(2*pow(mpi0,4) - 
      2*pow(mpi0,2)*s + pow(s,2)));
    //(s*(m12s*m12s*s-pow(mA*mA-s,2))-2*ma*ma*s*(m12s*mpi0*mpi0+pow(mpi0*mpi0-s,2))+pow(ma,4)*(2*pow(mpi0,4)-2*mpi0*mpi0*s+s*s));
}
*/
double Axion_Dark_Photon::pi0_decay_amplitude2_c(double m12s, double s, double m0, double m1, double m2, double m3){
    return pi0_decay_amplitude2(m12s,s,m1,m2,m3);
}

double Axion_Dark_Photon::pi0_decay_amplitude2_b(double m12, double CosTheta, double m0, double m1, double m2, double m3){
    return pi0_decay_amplitude2(m12*m12,Three_Body_Decay_Space::Cos_Theta_to_m23s(m12*m12, CosTheta, m0, m1, m2, m3),m1,m2,m3);
}

//m12s=m12^2, s= m23^2
double Axion_Dark_Photon::pi0_decay_width(double m12s, double s, double mgamma, double ma, double mA){
    return pi0_decay_amplitude2(m12s,s,mgamma,ma,mA)/(32*pow(2*pi*MASS_PION_NEUTRAL,3));
}

//This one works with Cos(theta) and m12 rather than m12^2 and m23.
double Axion_Dark_Photon::pi0_decay_width_2(double m12, double CosTheta, double m0, double m1, double m2, double m3){
    return pi0_decay_amplitude2(m12*m12,Three_Body_Decay_Space::Cos_Theta_to_m23s(m12*m12, CosTheta, m0, m1, m2, m3),m1,m2,m3);
}

//Warning, does not average over initial spins of DP!
//Assumes effectively massless axion! Would not be hard to update to include massive axion, later time.
double Axion_Dark_Photon::Axion_DP_electron_Amp(double s, double t, double mA, double mR){
    return -pow(Gagpg*G_ELEC,2)*(pow(mA,4)*(2*mR*mR+t)-2*mA*mA*t*(mR*mR+s+t)+t*(2*pow(mR*mR-s,2)+2*s*t+t*t))/2.0/t/t;
}


//Factor of 2*MASS_ELECTRON from change of variable to E4 from t, dt = 2 me dE4.
double Axion_Dark_Photon::dsigma_a_to_DP(double E1lab, double E4, double mA, double mR){
    double s = two_to_two_scattering::s_lab(E1lab, 0, mR);
    if(s<pow(mA+mR,2)){
        return 0;
    }
    return 2*MASS_ELECTRON*1.0/64.0/pi/pow(mR*E1lab,2)*Axion_DP_electron_Amp(s,two_to_two_scattering::t_lab(E4,mR,mR),mA,mR);
}

//Factor of 2*MASS_ELECTRON from change of variable to t, dt = 2 me dE4.
//Averages over initial dark photon spins, extra 1/3 factor.
double Axion_Dark_Photon::dsigma_DP_to_a(double E1lab, double E4, double mA, double mR){
    double s = two_to_two_scattering::s_lab(E1lab, mA, mR);
    if(s<pow(mR,2)){
        return 0;
    }
    return 2*mR*1.0/3.0/64.0/pi/s/(pow(E1lab,2)-mA*mA)*Axion_DP_electron_Amp(s,two_to_two_scattering::t_lab(E4,mR,mR),mA,mR);
}

//Decay of the DP
double Axion_Dark_Photon::A_to_a_gamma_width(double ma, double mA){
    if(mA<ma)
        return 0.0;
    return pow(Gagpg,2)*pow(mA*mA-ma*ma,3)/(96*pi*pow(mA,3));
}

//Decay of DP to axion+electron+positron, me2 is a dummy variable, should be the same as me.
//p1 -> elec, p2 -> positron, p3 -> axion
double Axion_Dark_Photon::A_to_a_elec_pos_amplitude(double m12s, double m23s, double mA, double me, double me2, double ma){
    if(mA<ma+2*MASS_ELECTRON)
        return 0.0;
    return (pow(G_ELEC,2)*pow(Gagpg,2)*(m12s*
        (pow(m12s,2) + 2*m12s*m23s + 2*pow(m23s,2) - 2*m12s*pow(ma,2) - 
          2*m23s*pow(ma,2) + pow(ma,4) - 2*(m12s + m23s)*pow(mA,2) + 
          pow(mA,4)) - 2*(-pow(pow(ma,2) - pow(mA,2),2) + 
          m12s*(2*m23s + pow(ma,2) + pow(mA,2)))*pow(me,2) + 
       2*m12s*pow(me,4)))/pow(m12s,2);
}

//Currently ignores DP->ae^+e^-.
double Axion_Dark_Photon::A_width(double ma, double mA){
    return A_to_a_gamma_width(ma,mA);
}
//End of decay of the DP

bool Axion_Dark_Photon::Set_Model_Parameters(Parameter& par){
    if(par.Query_Map("axion_mass",mass_axion)){
        mass_dp = par.MassDP();
        epsilon = par.Epsilon();
        par.Query_Map("gagg", Gagg);
        par.Query_Map("gagpg", Gagpg);
        par.Query_Map("gagpgp", Gagpgp);
        return true;
    }
    return false;
}

void Axion_Dark_Photon::Report(std::ostream& out){
    out << mass_axion << " " << mass_dp << " " << epsilon << " " << Gagg << " " << Gagpg << " " << Gagpgp << " " << A_width(mass_axion, mass_dp) << " ";
}

void Axion_Dark_Photon::Report_Model(){
    cout << "Mass Axion = " << mass_axion << " GeV" << endl;
    cout << "Mass Dark Photon = " << mass_dp << " GeV" << endl;
    cout << "epsilon = " << epsilon << endl;
    cout << "G_axion-photon-photon = " << Gagg << endl;
    cout << "G_axion-photon-dp = " << Gagpg << endl;
    cout << "G_axion-dp-dp = " << Gagpgp << endl;

}

bool Axion_Dark_Photon::Prepare_Production_Channel(std::string prodchoice, std::string proddist, production_channel& prodchan, std::shared_ptr<DMGenerator>& DMGen, std::shared_ptr<Distribution>& Dist, double& Vnum, Parameter& par){

    Particle photon(0);
    photon.name = "Photon";

    Particle axion(mass_axion);
    axion.name="Axion";
    
    Particle dark_photon(mass_dp);
    dark_photon.name="Dark_Photon";

    string sig_choice = par.Signal_Channel();
    
    if(prodchoice=="pi0_decay"||prodchoice=="eta_decay"){
        double lifetime;
        
        Particle meson(0);
        
        std::function<double(double, double, double, double, double, double)> func;

        if(prodchoice=="pi0_decay"){
            if(MASS_PION_NEUTRAL<mass_axion+mass_dp){
                Vnum=0;
                DMGen = std::shared_ptr<DMGenerator>(new Do_Nothing_Gen());
                return true;
            }

            meson.name = "Pi0";
            meson.m = MASS_PION_NEUTRAL;
            meson.width=WIDTH_PION_NEUTRAL;
            Dist->set_mass(MASS_PION_NEUTRAL);
            lifetime=0;//Decay is effectively instantaneous.
            func = std::bind(&Axion_Dark_Photon::pi0_decay_amplitude2_c,this,_1,_2,_3,_4,_5,_6);
        }
        else if(prodchoice=="eta_decay"){
            if(MASS_ETA<mass_axion+mass_dp){
                Vnum=0;
                DMGen = std::shared_ptr<DMGenerator>(new Do_Nothing_Gen());
                return true;
            }
            meson.name = "eta";
            meson.m = MASS_ETA;
            meson.width=WIDTH_ETA;
            Dist->set_mass(MASS_ETA);
            lifetime=0;//Decay is effectively instantaneous.
            func = std::bind(&Axion_Dark_Photon::eta_decay_amplitude2_c,this,_1,_2,_3,_4,_5,_6);
        }
        else{
            std::cerr << "Invalid production channel" << prodchoice << "\n" << endl;
        }


        std::shared_ptr<Three_Body_Decay_Gen> tmp_gen(new Three_Body_Decay_Gen(meson,photon,axion,dark_photon,prodchoice,lifetime,func));
        tmp_gen->d1=false;

        //Activate the decay of the DP unless we are looking for the decay of the DP
        
        if((sig_choice!="DP_Signal_Decay" and sig_choice!="DP_e_ep_a_Decay" and sig_choice!="DP_gamma_a_Decay")&&A_width(mass_axion, mass_dp)>0){
            cout << "Turning on decay of the Dark Photon with lifetime of " << hbar/A_width(mass_axion, mass_dp) << " seconds\n";
            std::shared_ptr<Two_Body_Decay_Gen> invis_dec(new Two_Body_Decay_Gen(A_to_a_gamma_width(mass_axion, mass_dp)/A_width(mass_axion, mass_dp),mass_dp,string("Dark_Photon"),photon,axion,hbar/A_width(mass_axion, mass_dp)));
            //Do not care about photons!
            invis_dec->d1=false;
            invis_dec->record_parent=false;
            tmp_gen->Toggle_Daughter_Decay(3,invis_dec);
        }
        //Probably going to have to wait to implement this properly, no time at the moment. Each of the generators will either do nothing, or perform a decay. 
        //Stopgap is everything does nothing for the time being, later on they'll have proper lifetimes.
        Vnum=tmp_gen->BranchingRatio()*par.Protons_on_Target()*prodchan.Meson_Per_Pi0()*par.Pi0_per_POT();
        if(tmp_gen->Channel_Name()==""){
            tmp_gen->Set_Channel_Name(prodchoice);
        }
        DMGen = tmp_gen;
        return true;
    }
    else if(prodchoice=="brem_dp" or prodchoice=="brem_axion"){
        double part_mass;
        
        if(prodchoice=="brem_dp"){
            part_mass=mass_dp;
            if(sig_choice!="DP_Signal_Decay"&&A_width(mass_axion, mass_dp)>0){
                cout << "Turning on decay of the Dark Photon with lifetime of " << hbar/A_width(mass_axion, mass_dp) << " seconds\n";
                std::shared_ptr<Two_Body_Decay_Gen> invis_dec(new Two_Body_Decay_Gen(A_to_a_gamma_width(mass_axion, mass_dp)/A_width(mass_axion, mass_dp),mass_dp,string("Dark_Photon"),photon,axion,hbar/A_width(mass_axion, mass_dp)));
                invis_dec->record_parent=false;
                DMGen = invis_dec;
            }
            else{
                std::shared_ptr<Do_Nothing_Gen> do_not(new Do_Nothing_Gen(string("Dark_Bremsstrahlung"),string("Dark_Photon")));
                DMGen = do_not;
            }
        }
        else if(prodchoice=="brem_axion"){
            part_mass=mass_axion;
            std::shared_ptr<Do_Nothing_Gen> do_not(new Do_Nothing_Gen(string("Dark_Bremsstrahlung"),string("Axion")));
            DMGen = do_not;
        }

        Dist->set_mass(part_mass);
        Vnum=par.Protons_on_Target()*prodchan.Num_per_pot();
        return true;
    }
    //Currently Reactor Production considers POT to be the number of V's produced.
    else if(prodchoice=="Reactor_Production"){
        if(sig_choice=="DP_gamma_a_Decay" or sig_choice=="DP_e_ep_a_Decay" or sig_choice=="DP_Signal_Decay"){
            std::shared_ptr<Do_Nothing_Gen> do_not(new Do_Nothing_Gen(string("Reactor_Production"),string("Dark_Photon")));
            DMGen = do_not;
            Dist->set_mass(mass_dp);
            Vnum=par.Protons_on_Target();
            return true;
        }
    }
    return false;
}

bool Axion_Dark_Photon::Prepare_Signal_Channel(Parameter& par){
    string sig_choice = par.Signal_Channel();
    if(sig_choice=="Electron_Scatter"){
        
        sig_part_vec.push_back(string("Axion"));
        sig_part_vec.push_back(string("Dark_Photon"));

        std::shared_ptr<Two_to_Two_Scatter> ttts(new Two_to_Two_Scatter());

        ttts->set_energy_limits(par.Min_Scatter_Energy(),par.Max_Scatter_Energy());
        ttts->set_angle_limits(par.Max_Angle(),par.Min_Angle());

        Particle electron(MASS_ELECTRON);
        electron.name="Electron";
        Particle axion_r(mass_axion);
        axion_r.name = "Recoil_Axion";
        Particle dark_photon_r(mass_dp);
        dark_photon_r.name = "Recoil_Dark_Photon";
        
        if(A_width(mass_axion, mass_dp)>0){
            Particle axion_f(mass_axion);
            axion_f.name="Decay_Axion";
            Particle photon_f(0);
            photon_f.name = "Decay_Photon";
            std::shared_ptr<Two_Body_Decay_Gen> invis_dec(new Two_Body_Decay_Gen(A_to_a_gamma_width(mass_axion, mass_dp)/A_width(mass_axion, mass_dp),mass_dp,string("Recoil_Dark_Photon"),photon_f,axion_f,hbar/A_width(mass_axion, mass_dp)));
            invis_dec->record_parent=false;
            ttts->add_decay(string("Recoil_Dark_Photon"),invis_dec);
        }

        //axion+e->dp+e
        std::shared_ptr<Linear_Interpolation> axion_cross;
        std::shared_ptr<Linear_Interpolation> axion_cross_max;
        
        //double s_test = two_to_two_scattering::s_lab(0.07,0,MASS_ELECTRON);
        //double t_test = two_to_two_scattering::t_lab(0.0318287,MASS_ELECTRON,MASS_ELECTRON);
        //cout << "s=" << s_test << " t=" << t_test << " ";

        //cout << "Amplitude = " << Axion_DP_electron_Amp(two_to_two_scattering::s_lab(0.07,0,MASS_ELECTRON),t_test, mass_dp, MASS_ELECTRON) << endl;

        //cout << "disgma=" << dsigma_a_to_DP(0.07, 0.0318287, mass_dp, MASS_ELECTRON) << endl;

        std::function<double(double,double)> f_a_to_dp = std::bind(&Axion_Dark_Photon::dsigma_a_to_DP,this,_1,_2,mass_dp,MASS_ELECTRON);
        std::function<double(double)> ER_min_ax = std::bind(&Two_to_Two_Scatter::scatmin,*ttts,_1,mass_axion,MASS_ELECTRON,mass_dp,MASS_ELECTRON);
        std::function<double(double)> ER_max_ax = std::bind(&Two_to_Two_Scatter::scatmax,*ttts,_1,mass_axion,MASS_ELECTRON,mass_dp,MASS_ELECTRON);
        
        Prepare_Cross_Section(f_a_to_dp, ER_min_ax, ER_max_ax, axion_cross,axion_cross_max, mass_axion, par.Max_DM_Energy(),par.EDM_RES());
        
        function<double(double)> axion_cross_func = bind(&Linear_Interpolation::Interpolate,axion_cross,_1);
        function<double(double)> axion_cross_max_func = bind(&Linear_Interpolation::Interpolate,axion_cross_max,_1);

        //dp+e->axion+e
        std::shared_ptr<Linear_Interpolation> dp_cross;
        std::shared_ptr<Linear_Interpolation> dp_cross_max;

        std::function<double(double,double)> f_dp_to_a = std::bind(&Axion_Dark_Photon::dsigma_DP_to_a,this,_1,_2,mass_dp,MASS_ELECTRON);
        std::function<double(double)> ER_min_dp = std::bind(&Two_to_Two_Scatter::scatmin,*ttts,_1,mass_dp,MASS_ELECTRON,mass_axion,MASS_ELECTRON);
        std::function<double(double)> ER_max_dp = std::bind(&Two_to_Two_Scatter::scatmax,*ttts,_1,mass_dp,MASS_ELECTRON,mass_axion,MASS_ELECTRON);

        Prepare_Cross_Section(f_dp_to_a, ER_min_dp, ER_max_dp,dp_cross,dp_cross_max, mass_dp, par.Max_DM_Energy(),par.EDM_RES());

        function<double(double)> dp_cross_func = bind(&Linear_Interpolation::Interpolate,*dp_cross,_1);
        function<double(double)> dp_cross_max_func = bind(&Linear_Interpolation::Interpolate,*dp_cross_max,_1);


        double ENtot =(par.Get_Detector())->ENtot();
        ttts->add_channel(dark_photon_r, electron, MASS_ELECTRON, axion_cross_func, axion_cross_max_func, f_a_to_dp, ENtot,"Axion");
        ttts->add_channel(axion_r, electron, MASS_ELECTRON, dp_cross_func, dp_cross_max_func, f_dp_to_a, ENtot,"Dark_Photon");

        Sig_list.push_back(ttts);

        return true;
    }//Searching for DP in signal decays
    else if(sig_choice=="DP_e_ep_a_Decay"){

        sig_part_vec.push_back(string("Dark_Photon"));

        Particle dark_photon(mass_dp);
        dark_photon.name = "Dark_Photon";
        dark_photon.width = A_width(mass_axion,mass_dp);
        Particle electron(MASS_ELECTRON);
        electron.name = "Decay_Electron";
        Particle positron(MASS_ELECTRON);
        positron.name = "Decay_Positron";
        Particle axion(mass_axion);
        axion.name = "Decay_Axion";

        std::function<double(double,double,double,double,double,double)> func = std::bind(&Axion_Dark_Photon::A_to_a_elec_pos_amplitude,this,_1,_2,_3,_4,_5,_6);

        std::shared_ptr<Three_Body_Decay_Gen> dp_3decay_gen(new Three_Body_Decay_Gen(dark_photon,electron,positron,axion,string("Decay_Electron_Positron_Axion"),0,func));

        dp_3decay_gen->record_parent = false;

        vector<std::shared_ptr<DMGenerator> > dec_vec;
        dec_vec.push_back(dp_3decay_gen);

        double lifetime = hbar/A_width(mass_axion, mass_dp);

        std::shared_ptr<SignalDecay_2> sig_dec(new SignalDecay_2(lifetime, dec_vec));
        Sig_list.push_back(sig_dec);
        return true;
    }
    else if(sig_choice=="DP_gamma_a_Decay"){

        sig_part_vec.push_back(string("Dark_Photon"));

        Particle dark_photon(mass_dp);
        dark_photon.name = "Dark_Photon";
        dark_photon.width = A_width(mass_axion,mass_dp);
        Particle photon(0);
        photon.name = "Decay_Photon";
        Particle axion(mass_axion);
        axion.name = "Decay_Axion";

        std::shared_ptr<Two_Body_Decay_Gen> dp_2decay_gen(new Two_Body_Decay_Gen(A_to_a_gamma_width(mass_axion, mass_dp)/A_width(mass_axion, mass_dp),mass_dp, string("Dark_Photon"), photon, axion,0));//Set lifetime to zero to force immediate decay!

        dp_2decay_gen->record_parent = false;

        vector<std::shared_ptr<DMGenerator> > dec_vec;
        dec_vec.push_back(dp_2decay_gen);

        double lifetime = hbar/dark_photon.width;

        std::shared_ptr<SignalDecay_2> sig_dec(new SignalDecay_2(lifetime, dec_vec));
        Sig_list.push_back(sig_dec);
        return true;
    }
    return false;
}


/*
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
}*/
