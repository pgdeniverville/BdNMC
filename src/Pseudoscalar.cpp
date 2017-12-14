#include "Model.h"
#include "constants.h"
#include "Kinematics.h"
#include "Drell_Yan_Gen.h"

#include <math.h>

using std::cerr; using std::vector;
using std::string; using std::function;
using std::max; using std::min;
using std::cout; using std::endl;

const double me = MASS_ELECTRON;
const double gq = G_ELEC; 

const string g_chi_key="g_chi";
const string dark_mediator_mass_key="dark_mediator_mass";
const string dark_matter_mass_key="dark_matter_mass";

using namespace std::placeholders;

using namespace elastic_scattering;

bool Pseudoscalar::Set_Model_Parameters(Parameter& par){
    if(par.Query_Map(g_chi_key,gchi)&&(par.Query_Map(dark_matter_mass_key,mchi))&&par.Query_Map(dark_mediator_mass_key, ma)){
        return true;
    }
    return false;
}

bool Pseudoscalar::Prepare_Signal_Channel(Parameter& par){
    string sig_choice=par.Signal_Channel();//This is eventually going to get looped.
    
    if(sig_choice=="NCE_electron"){
        std::shared_ptr<Elastic_Scatter> els(new Elastic_Scatter);
        
        Particle electron(me); 
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

bool Pseudoscalar::Prepare_Production_Channel(std::string prodchoice, std::string proddist, production_channel& prodchan, std::shared_ptr<DMGenerator> DMGen, std::shared_ptr<Distribution>, double& Vnum, Parameter& par){
    if(prodchoice=="Drell_Yan"){
        function<double(double, double, double, double, double)> dsig_dEk = std::bind(&Pseudoscalar::dsigma_dEk_qq_to_chichi,this,_1,par.Beam_Energy(),_2,_3,_4,_5);
        std::shared_ptr<Drell_Yan_Gen> tmp_gen(new Drell_Yan_Gen(mchi,par.Beam_Energy(),dsig_dEk,par.Target_P_Num(),par.Target_N_Num(),prodchan.Proton_PDF_File(),prodchan.Neutron_PDF_File(),prodchoice,1e-3));
        Vnum=tmp_gen->Interaction_Probability_Density()*par.Protons_on_Target()*par.Target_Length();
        cout << "Compare " << Vnum << " with V_num estimate " << tmp_gen->Sig_P()*par.Protons_on_Target()/par.P_Cross()*convGeV2cm2*(par.Target_N_Num()+par.Target_P_Num()) << endl;
        DMGen = tmp_gen;
        return true;
    }
    return false;
}

//Differential scattering cross section for chi+e->chi+e.
//gq is electron charge?
double Pseudoscalar::dsigma_dEf_electron(double Ei, double Ef){
    return 1.0/8.0/pi*pow(gchi*gq,2)/(Ei*Ei-mchi*mchi)*pow(Ef-me,2)/pow(2*me*me -2*me*Ef-ma*ma,2);
}

//Integrated version of dsigma_dEf_electron.
double Pseudoscalar::sigma_Ef_electron(double Ei, double Ef){
    return -(pow(gchi,2)*pow(gq,2)*(-2*Ef*me + 2*pow(me,2) + pow(ma,4)/(pow(ma,2) + 2*(Ef - me)*me) + 2*pow(ma,2)*log(pow(ma,2) + 2*(Ef - me)*me)))/(64.*(pow(Ei,2) - pow(mchi,2))*pow(me,2)*pi);
}

double Pseudoscalar::dsig_max(double Ei){
    return dsigma_dEf_electron(Ei, max(scat_min, E2fMin(Ei,mchi,me)));
}

double Pseudoscalar::sigma_tot_electron(double Ei){
    return sigma_Ef_electron(Ei,min(scat_max,E2fMin(Ei,mchi,me)))-sigma_Ef_electron(Ei,max(scat_min, E2fMin(Ei,mchi,me)));
}

double Pseudoscalar::dsigma_dEk_qq_to_chichi(double Ek, double EA, double x, double y, double gf, double MASS){
    double s_hat = pow(x*MASS_PROTON,2) + pow(y*MASS,2) + 2*x*y*EA*MASS;
    return 1.0/16/pi*pow(gf*gchi,2)/(s_hat-ma*ma)*s_hat/(x*EA);
}