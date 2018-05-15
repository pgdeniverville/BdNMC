#include "Proton_Brem.h"
#include "constants.h"
#include "Random.h"
#include "Integrator.h"

using std::bind;
using std::cout;
using std::endl;
using namespace std::placeholders;

Proton_Brem::Proton_Brem(double Beam_E, std::function<double(double, double)> splitting_function, double m_mediator, double ptmax, double zmax, double zmin, std::string &mode, std::shared_ptr<DMgenerator> V_dec, double ptmin){

    MA=m_mediator;
    PTMIN=ptmin;
    PTMAX=ptmax;
    ZMIN=zmin;
    ZMAX=zmax;
    dsig = splitting_function;
    Beam_Energy=Beam_E;
    V_decay = V_dec;

    calc_V_prod_rate();

}

bool GenDM(std::list<Particle>& vec, std::function<double(Particle&)> det_int, Particle& part){


}
 

void Proton_Brem::calc_V_prod_rate(){
    std::function<double(double, double)> func = std::bind(&Proton_Brem_Distribution::d2N_proton_brem_to_V,this,_1,_2);
    if(ZMIN>ZMAX){
        vprodrate=0;
    }
    else
        vprodrate = SimpsonCubature(func,ZMIN,ZMAX,100,PTMIN,PTMAX,100);//Will need to tweak this. Hopefully come up with a more general algorithm at some point.

    max_prod=d2N_proton_brem_to_V(ZMIN,PTMIN);
};

double Proton_Brem::d2N_proton_brem_to_out(double z, double pt2){
    return sigmapp(2*mp*(Beam_Energy-sqrt(MA*MA+pt2+z*z*(pow(Beam_Energy,2)-mp*mp))))/sigmapp(2*mp*Beam_Energy)*dsig(z,pt2);
}

double Proton_Brem::sigmapp(double s){
    return Hpp*pow(log(s/sppM),2)+Ppp+R1pp*pow(s/sppM,-eta1pp)-R2pp*pow(s/sppM,-eta2pp);
}

void Proton_Brem::set_fit_parameters(production_channel &par){
    par.query_dist_param("rD",rD);
    // par.query_dist_param("mD2",mD2);
    par.query_dist_param("H",Hpp);
    par.query_dist_param("M",Mpp);
    par.query_dist_param("eta1",eta1pp);
    par.query_dist_param("eta2",eta2pp);
    par.query_dist_param("R1pp",R1pp);
    par.query_dist_param("R2pp",R2pp);
    par.query_dist_param("Ppp",Ppp);
    sppM = pow(2*mp+Mpp,2);
    calc_V_prod_rate();
    max_prod = d2N_proton_brem_to_out(ZMIN,PTMIN);
}

void Proton_Brem::sample_particle(Particle &part){
    double mom, theta, phi;
    sample_momentum(mom, theta, phi);
    part.ThreeMomentumPolar(mom, theta, phi);
}

void Proton_Brem::sample_momentum(double &pmom, double &theta, double &phi){
    double z, pt, hold;
    while(true){
        z = Random::Flat(ZMIN,ZMAX);
        pt = Random::Flat(PTMIN,PTMAX);
        if((hold=d2N_proton_brem_to_out(z,pt))>Random::Flat()*max_prod){
            if(max_prod<hold)
                max_prod=hold;
            break;
        }
    }
    phi = Random::Flat(0,2*pi);
    double pz = sqrt(pow(Beam_Energy,2)-mp*mp)*z;
    theta=atan2(pt,pz);
    pmom = sqrt(pz*pz+pt*pt);
}