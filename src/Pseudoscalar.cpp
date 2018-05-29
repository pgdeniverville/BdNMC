#include "Model.h"
#include "constants.h"
#include "Kinematics.h"
#include "Drell_Yan_Gen.h"
#include "branchingratios.h"
#include "DMgenerator.h"
#include "Proton_Brem.h"

#include <memory>
#include <math.h>

using std::cerr; using std::vector;
using std::string; using std::function;
using std::max; using std::min;
using std::cout; using std::endl;

using std::shared_ptr;

const double me = MASS_ELECTRON;
//const double gq = G_ELEC; 

const string g_chi_key="g_chi";
const string g_quark_key="g_quark";
const string g_electron_key = "g_pseudoscalar_electron";
const string dark_mediator_mass_key="dark_mediator_mass";
const string dark_matter_mass_key="dark_matter_mass";

using namespace std::placeholders;

using namespace elastic_scattering;

bool Pseudoscalar::Set_Model_Parameters(Parameter& par){
    if(par.Query_Map(g_chi_key,gchi)&&(par.Query_Map(dark_matter_mass_key,mchi))&&par.Query_Map(dark_mediator_mass_key, ma)&&par.Query_Map(g_quark_key,gq)&&par.Query_Map(g_electron_key,gae)){
        return true;
    }
    return false;
}

bool Pseudoscalar::Prepare_Signal_Channel(Parameter& par){
    string sig_choice=par.Signal_Channel();//This is eventually going to get looped.
    
    if(sig_choice=="NCE_electron"){
        std::shared_ptr<Elastic_Scatter> els(new Elastic_Scatter);
        
        Particle electron(me); 
        electron.name = "Electron";
        function<double(double, double)> dsig = bind(&Pseudoscalar::dsigma_dEf_electron,this,_1,_2);
        function<double(double)> sig_tot = bind(&Pseudoscalar::sigma_tot_electron,this,_1);
        function<double(double)> sig_max =  bind(&Pseudoscalar::dsig_max,this,_1);
        double ENtot = (par.Get_Detector())->ENtot();
        els->add_channel(electron,sig_tot,sig_max,dsig,ENtot);
        
        Sig_list.push_back(els);
        return true;
    }
    return false;
}

//Production channel comes from Parameter.h!
bool Pseudoscalar::Prepare_Production_Channel(std::string prodchoice, std::string proddist, production_channel& prodchan, std::shared_ptr<DMGenerator>& DMGen, std::shared_ptr<Distribution>& Dist, double& Vnum, Parameter& par){
    if(prodchoice=="Drell_Yan"){
//        cout << "ma=" << ma << " mx=" << mchi << endl;
//        cout << dsigma_dEk_qq_to_chichi(2,8.9,0.3,0.15,1,MASS_PROTON) << endl;
        //cout <<  "Test sigma for EA=20, x=0.3, y=0.15, MASS=0.938 : " << dsigma_hat_dt_qq_to_chi_chi(20, 0, 0.3, 0.15, 1, 0.938) << endl;
        function<double(double, double, double, double, double)> dsig_dEk = std::bind(&Pseudoscalar::dsigma_dEk_qq_to_chichi,this,_1,par.Beam_Energy(),_2,_3,_4,_5);
        std::shared_ptr<Drell_Yan_Gen> tmp_gen(new Drell_Yan_Gen(mchi,par.Beam_Energy(),dsig_dEk,par.Target_P_Num(),par.Target_N_Num(),prodchan));
        Vnum=tmp_gen->Interaction_Cross_Section()*par.Protons_on_Target()*par.Target_Length()*par.Target_N_Dens()*convGeV2cm2*m_to_cm;
        cout << "Compare " << Vnum << " with V_num estimate " << tmp_gen->Sig_P()*par.Protons_on_Target()/par.P_Cross()*convGeV2cm2 << endl;
        cout << "Compare " << Vnum << " with V_num estimate " << tmp_gen->Interaction_Cross_Section()*par.Protons_on_Target()/par.P_Cross()/8.0*convGeV2cm2 << endl;
        cout << tmp_gen->Sig_P() << " " << par.P_Cross() << " " << tmp_gen->Interaction_Cross_Section() << " " << par.Target_N_Dens() << endl;
        cout << "Note overestimate due to 1.75 interaction lengths of material.\n";
        DMGen = tmp_gen;
        return true;
    }
    //This might need to be proddist?
    else if(prodchoice=="Proton_Bremsstrahlung"){
        //Not sure if this is the right coupling.
        function<double(double, double)> dsig_dpt2dz = std::bind(brem_split_pseudoscalar,_1,_2,ma,gq*proton_form_factor(ma*ma));
        
        std::shared_ptr<DMGenerator> dmgen;
        if(ma<mchi*2){
            double branching_ratio_to_neutrinos=1;
            Particle neutrino(0);
            neutrino.name = "Neutrino";
            dmgen = shared_ptr<DMGenerator>(new Two_Body_Decay_Gen(branching_ratio_to_neutrinos,ma,"Pseudoscalar",neutrino, neutrino, Gamma_pseudoscalar_to_2fermion(gae,ma,0)));
            sig_part_name = "Neutrino";
        }
        else{
            double branching_ratio_to_invisible=1;
            Particle dm(mchi);
            dm.name = "DM";
            dmgen = shared_ptr<DMGenerator>(new Two_Body_Decay_Gen(branching_ratio_to_invisible,ma,"Pseudoscalar",dm, dm));
            sig_part_name = "DM";
        }
        Particle med(ma);
        med.name = "Pseudoscalar";
        string prod_chan_name = "Proton_Bremsstrahlung_Pseudoscalar";
        std::shared_ptr<Proton_Brem> tmp_gen(new Proton_Brem(par.Beam_Energy(),dsig_dpt2dz,med,prodchan.ptmax(),prodchan.zmax(),prodchan.zmin(),prod_chan_name, dmgen, prodchan.ptmin()));
        Vnum=tmp_gen->BranchingRatio()*par.Protons_on_Target();
        DMGen = tmp_gen;
        return true;
    }
    return false;
}
//Placeholder, no idea what it looks like!
double Pseudoscalar::proton_form_factor(double q2){
    return 1.0;
}

//Differential scattering cross section for chi+e->chi+e.
//gq is electron charge?
double Pseudoscalar::dsigma_dEf_electron(double Ei, double Ef){
    return 1.0/8.0/pi*pow(gchi*gae,2)/(Ei*Ei-mchi*mchi)*me*pow(Ef-me,2)/pow(2*me*me -2*me*Ef-ma*ma,2);
}

//Integrated version of dsigma_dEf_electron.
double Pseudoscalar::sigma_Ef_electron(double Ei, double Ef){
    return -(pow(gchi*gae,2)*(-2*Ef*me + 2*pow(me,2) + pow(ma,4)/(pow(ma,2) + 2*(Ef - me)*me) + 2*pow(ma,2)*log(pow(ma,2) + 2*(Ef - me)*me)))/(64.*(pow(Ei,2) - pow(mchi,2))*pow(me,2)*pi);
}

double Pseudoscalar::dsig_max(double Ei){
    return dsigma_dEf_electron(Ei, max(scat_min, E2fMin(Ei,mchi,me)));
}

double Pseudoscalar::sigma_tot_electron(double Ei){
    return sigma_Ef_electron(Ei,min(scat_max,E2fMax(Ei,mchi,me)))-sigma_Ef_electron(Ei,max(scat_min,E2fMin(Ei,mchi,me)));
}

//gf has no effect!
double Pseudoscalar::dsigma_dEk_qq_to_chichi(double Ek, double EA, double x, double y, double gf, double MASS){
    double s_hat = annihilation_to_pair::shat(EA*x,x*MASS_PROTON,y*MASS);
    //return 1/3.0*1.0/16/pi*pow(gq*gchi,2)/pow(s_hat-ma*ma,2)*s_hat/(x*EA);
    //return (pow(gf,2)*pow(gchi,2)*s_hat*(-s_hat + pow(MASS_PROTON*x - MASS*y,2)))/(96.*MASS*(pow(E1cm_hat,2) - pow(x*MASS_PROTON,2))*pi*pow(pow(ma,2) - s_hat,2)*pow(x,2)*y);
    return -pow(gchi*gq,2)*s_hat*(pow(x*MASS_PROTON-y*MASS,2)-s_hat)/(96*pi*y*pow(x,2)*(pow(EA,2)-pow(MASS_PROTON,2))*MASS)/pow(ma*ma-s_hat,2);
}

double Pseudoscalar::dsigma_hat_dt_qq_to_chi_chi(double EA, double t, double x, double y, double qf, double MASS){
    double s_hat=annihilation_to_pair::shat(EA*x,x*MASS_PROTON,y*MASS);
    //double E1cm_hat = (s_hat+pow(x*MASS_PROTON,2)-pow(y*MASS,2))/(2*sqrt(s_hat));
    return pow(gq*gchi,2)*s_hat*(s_hat-pow(MASS_PROTON*x-MASS*y,2))/192.0/pi/pow(MASS,2)/pow(x*y,2)/(pow(EA,2)-pow(MASS_PROTON,2))/pow(ma*ma-s_hat,2);
}

double Pseudoscalar::sigma_hat_tot_qq_to_chi_chi(double EA, double x, double y, double qf, double MASS){
    double t0=annihilation_to_pair::t0(EA*x,MASS_PROTON*x,MASS*y,mchi);
    double t1=annihilation_to_pair::t1(EA*x,MASS_PROTON*x,MASS*y,mchi);
    if(t0<t1)
        return 0;
    //differential equation is independent of t.
    return (t0-t1)*dsigma_hat_dt_qq_to_chi_chi(EA, 0, x, y, qf, MASS);
}

void Pseudoscalar::Report(std::ostream& out){
    out << ma << " " << mchi << " " << gchi << " " << gq << " " << gae << " ";
/*    cout << "Self-Check!" << endl;
    cout << dsigma_dEf_electron(1, 0.4) << endl;
    cout << sigma_Ef_electron(1, 0.4) << endl;
    cout << sigma_tot_electron(1) << endl;
    cout << "Efmax,min=" << min(scat_max,E2fMax(1,mchi,me)) << " , " << max(scat_min,E2fMin(1,mchi,me)) << endl;
    cout << "gq=" << gq << " gchi=" << gchi << " ma=" << ma << " me=" << me << " mx" << mchi << endl;*/
}

void Pseudoscalar::Report_Model(){
    cout << "Dark pseudoscalar mediator mass = " << ma << " GeV" << endl;  
    cout << "Dark matter mass = " << mchi << " GeV" << endl; 
    cout << "g_chi = " << gchi << endl;
    cout << "g_electron = " << gae << endl;
    cout << "g_quark = " << gq << endl; 
    //    cout << "kappa = " << kappa << endl;
}
