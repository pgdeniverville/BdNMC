#include "branchingratios.h"
#include "Integrator.h"
#include "Random.h"
#include "Kinematics.h"
#include <functional>
#include <cmath>
#include <iostream>
#include <iomanip>

using std::cout; using std::endl;

const double mpion = 0.1349766;
const double mpim = 0.139570;
const double meta = 0.547853;
const double momega = 0.782;
const double mrho = 0.770;
const double mphi = 1.020;
const double metap = 0.95778;
const double Br_etap_to_gamma_gamma=0.022;
const double Br_etap_to_gamma_rho=0.291;
const double Br_etap_to_gamma_omega=0.0275;
const double Br_phi_to_e_e=2.9e-4;
const double Br_rho_to_e_e = 4.72e-5;
const double brpi0to2gamma = 0.98823;
const double bretato2gamma = 0.3941;
const double Br_omega_to_e_e = 7.28e-5;
const double pi = 3.14159;
const double melec = 0.000511;
const double mmuon = 0.1057;
const double alphaem = 1/137.036;
const double etafactor = 0.61;

/********************
 *KINETIC COUPLING
 ********************/



double brpi0toVgamma(double mv, double mx, double kappa, double alphaD){
    if(mv>mpion)
        return 0;
    return 2*pow(kappa,2)*pow(1-mv*mv/pow(mpion,2),3)*brpi0to2gamma;
}

double bretatoVgamma(double mv, double mx, double kappa, double alphaD){
    if(mv>meta)
        return 0;
    return 2*pow(kappa,2)*pow(1-mv*mv/pow(meta,2),3)*bretato2gamma;
}

double rhofunc(double x){
    if(x<0.77)
        return exp(-14*(0.77-x));
    else
        return exp(20*(0.77-x));
}


double brrho_to_V_1(double mv, double mx, double kappa, double alphaD){
    return 3.0*pow(kappa,2.0)/alphaem*Br_rho_to_e_e*pow(mrho,3)*GammaV_to_dm_dm(mv,mx,kappa,alphaD)/
        (pow(pow(mrho,2)-mv*mv,2)+pow(mv*GammaV(mv,mx,kappa,alphaD),2));
}

double brrho_to_V(double mv, double mx, double kappa, double alphaD){
    if(mv>1.6)
        return 0.0;
    return brrho_to_V_1(0.77, mx, kappa, alphaD)*rhofunc(mv);
}

double bromega_to_V(double mv, double mx, double kappa, double alphaD){
    if(mv-momega>0.1||momega-mv>0.1)
        return 0.0;
    return 3.0*pow(kappa,2)/alphaem*Br_omega_to_e_e*pow(momega,3)*GammaV_to_dm_dm(mv,mx,kappa,alphaD)/
        (pow(pow(momega,2)-mv*mv,2)+pow(mv*GammaV(mv,mx,kappa,alphaD),2));
}

//I could give this the off-shell treatment, but probably not worth it.
double br_eta_prime_to_V(double mv, double mx, double kappa, double alphaD){
    if(mv>=metap)
        return 0;
    double tot = 2*pow(1-pow(mv/metap,2),3)*Br_etap_to_gamma_gamma;
    if(mv<metap-mrho)
        tot+=pow(1-pow(mv/(metap-mrho),2),1.5)*Br_etap_to_gamma_rho;
    if(mv<metap-momega)
        tot+=pow(1-pow(mv/(metap-momega),2),1.5)*Br_etap_to_gamma_omega;
    return tot*pow(kappa,2);
}

double brphi_to_V(double mv, double mx, double kappa, double alphaD){
    if(mv-mphi>0.15||mphi-mv>0.15)
        return 0.0;
    return 3.0*pow(kappa,2)/alphaem*Br_phi_to_e_e*pow(mphi,3)*GammaV_to_dm_dm(mv,mx,kappa,alphaD)/
        (pow(pow(mphi,2)-mv*mv,2)+pow(mv*GammaV(mv,mx,kappa,alphaD),2));
}

double brmasstoVgamma(double mass, double mv, double mx, double kappa, double alphaD){
    if(mv>mass)
        return 0;
    return 2*pow(kappa,2)*pow(1-mv*mv/pow(mass,2),0.5);//scaling less dubious
}

double GammaV(double mv, double mx, double kappa, double alphaD){
    double term = 0;
    if(mv>2*mx){
        term = alphaD*(mv*mv-4*mx*mx)*sqrt(mv*mv/4-mx*mx);
    }
    if(mv>2*melec){
        term += 4*pow(kappa,2)*alphaem*(2*pow(melec,2)+mv*mv)*sqrt(mv*mv/4-pow(melec,2));
    }
    if(mv>2*mmuon){
        term += 4*pow(kappa,2)*alphaem*(2*pow(mmuon,2)+mv*mv)*sqrt(mv*mv/4-pow(mmuon,2));//I need to multiply this by the RRatio, but it's a minor overall effect when V->DM+DM is available;
    }
    return 1.0/(6.0*mv*mv)*(term);
}

double GammaV_to_dm_dm(double mv, double mx, double kappa, double alphaD){
    if(mv<=2*mx)
        return 0;
    return alphaD*(mv*mv-4*mx*mx)*sqrt(mv*mv/4.0-mx*mx)/(6.0*mv*mv);
}

double brV_to_dm_dm(double mv, double mx, double kappa, double alphaD){
    if(mv<=2*mx)
        return 0;
    return GammaV_to_dm_dm(mv, mx, kappa, alphaD)/GammaV(mv, mx, kappa, alphaD);
}

double d2brpi0_to_gamma_dm_dm(double mv, double mx, double kappa, double alphaD, double s, double theta){
    return sqrt(1-4*mx*mx/s)*pow(pow(mpion,2)-s,3)*(s-4*mx*mx)*alphaD*pow(kappa,2)*pow(sin(theta),3)/(8*pow(mpion,6)*pi*(pow(s-mv*mv,2)+pow(mv*GammaV(mv,mx,kappa,alphaD),2)));
}

double dbrpi0_to_gamma_dm_dm(double mv, double mx, double kappa, double alphaD, double s){
    return sqrt(1.0-4.0*mx*mx/s)*pow(pow(mpion,2)-s,3)*(s-4.0*mx*mx)*alphaD*pow(kappa,2)/(6.0*pow(mpion,6)*pi*(pow(s-mv*mv,2)+pow(mv*GammaV(mv,mx,kappa,alphaD),2)));
}

double brpi0_to_gamma_dm_dm(double mv, double mx, double kappa, double alphaD){
    using namespace std::placeholders;
    auto dbr = std::bind(dbrpi0_to_gamma_dm_dm,mv,mx,kappa,alphaD,_1);
    return DoubleExponential_adapt(dbr,4*mx*mx,pow(mpion,2),100,0.1,1e-4);
}

double d2breta_to_gamma_dm_dm(double mv, double mx, double kappa, double alphaD, double s, double theta){
    return bretato2gamma/brpi0to2gamma*sqrt(1-4*mx*mx/s)*pow(pow(meta,2)-s,3)*(s-4*mx*mx)*alphaD*pow(kappa,2)*pow(sin(theta),3)/(8*pow(meta,6)*pi*(pow(s-mv*mv,2)+pow(mv*GammaV(mv,mx,kappa,alphaD),2)));
}

double dbreta_to_gamma_dm_dm(double mv, double mx, double kappa, double alphaD, double s){
    return bretato2gamma/brpi0to2gamma*sqrt(1.0-4.0*mx*mx/s)*pow(pow(meta,2)-s,3)*(s-4.0*mx*mx)*alphaD*pow(kappa,2)/(6.0*pow(meta,6)*pi*(pow(s-mv*mv,2)+pow(mv*GammaV(mv,mx,kappa,alphaD),2)));
}

double breta_to_gamma_dm_dm(double mv, double mx, double kappa, double alphaD){
    using namespace std::placeholders;
    auto dbr = std::bind(dbreta_to_gamma_dm_dm,mv,mx,kappa,alphaD,_1);
	return DoubleExponential_adapt(dbr,4*mx*mx,pow(meta,2),100,0.1,1e-4);
}

//This is a nebulous placeholder function. I don't know how well it will work!
double brmasstoV(double mass, double mv, double mx, double kappa, double alphaD){
    if(mv>mass){
        return 0;
    }
    return pow(kappa,2)*pow(1-mv*mv/pow(mass,2),3)*0.5;
}

double d2brmass_to_dm_dm(double mmeson, double mv, double mx, double kappa, double alphaD, double s, double theta){
//   return alphaD*pow(kappa,2)*sqrt(1-4*mx*mx/s)*(1-s/pow(mpion,2))*(pow(pow(mpion,2)-s,2)-16*pow(TriangleFunction(mpion,s,0)*TriangleFunction(s,mx,mx)*cos(theta),2)*s)*Sin(theta)/(8*mpion^4*pi*(pow(s-mv*mv,2)+pow(mv*GammaV(mv,mx,kappa,alphaD),2)))
    return sqrt(1-4*mx*mx/s)*pow(pow(mmeson,2)-s,3)*(s-4*mx*mx)*alphaD*pow(kappa,2)*pow(sin(theta),3)/(8*pow(mmeson,6)*pi*(pow(s-mv*mv,2)+pow(mv*GammaV(mv,mx,kappa,alphaD),2)));
}

double dbrmass_to_dm_dm(double mmeson, double mv, double mx, double kappa, double alphaD, double s){
    return sqrt(1.0-4.0*mx*mx/s)*pow(pow(mmeson,2)-s,3)*(s-4.0*mx*mx)*alphaD*pow(kappa,2)/(6.0*pow(mmeson,6)*pi*(pow(s-mv*mv,2)+pow(mv*GammaV(mv,mx,kappa,alphaD),2)));
}

double brmass_to_dm_dm(double mmeson, double mv, double mx, double kappa, double alphaD){
    using namespace std::placeholders;
    auto dbr = std::bind(dbrmass_to_dm_dm,mmeson, mv,mx,kappa,alphaD,_1);
    return DoubleExponential_adapt(dbr,4*mx*mx,pow(mmeson,2),100,0.1,1e-4);
}

/********************
 *BARYONIC COUPLING
 ********************/


//Not implemented for masses above the three pion threshold.

double GammaVB_to_dm_dm(double mv, double mx, double kappa, double alphaD){
    if(mv>2*mx){
         return alphaD*(mv*mv-4*mx*mx)*sqrt(mv*mv/4-mx*mx)/(6*mv*mv);
    }
    return 0; 
}

double GammaVB(double mv, double mx, double kappa, double alphaD){
    return GammaVB_to_dm_dm(mv, mx, kappa, alphaD);
}

double brVB_to_dm_dm(double mv, double mx, double kappa, double alphaD){
    if(mv<=2*mx)
        return 0;
    return GammaVB_to_dm_dm(mv, mx, kappa, alphaD)/GammaVB(mv, mx, kappa, alphaD);
}

double brpi0toVBgamma(double mv, double mx, double kappa, double alphaD){
    if(mv>mpion)
        return 0;
    return 2*pow(sqrt(alphaD/alphaem)-kappa,2)*pow(1-mv*mv/pow(mpion,2),3)*brpi0to2gamma;
}

double d2brpi0_to_gamma_dm_dm_baryonic(double mv, double mx, double kappa, double alphaD, double s, double theta){
    return sqrt(1-4*mx*mx/s)*pow(pow(mpion,2)-s,3)*(s-4*mx*mx)*pow(alphaD,1)*pow(sqrt(alphaD/alphaem)-kappa,2)*pow(sin(theta),3)/(8*pow(mpion,6)*pi*(pow(s-mv*mv,2)+pow(mv*GammaVB(mv,mx,kappa,alphaD),2)));
}

double dbrpi0_to_gamma_dm_dm_baryonic(double mv, double mx, double kappa, double alphaD, double s){
    return sqrt(1.0-4.0*mx*mx/s)*pow(pow(mpion,2)-s,3)*(s-4.0*mx*mx)*pow(alphaD,1)*pow(sqrt(alphaD/alphaem)-kappa,2)/(6.0*pow(mpion,6)*pi*(pow(s-mv*mv,2)+pow(mv*GammaVB(mv,mx,kappa,alphaD),2)));
}

double brpi0_to_gamma_dm_dm_baryonic(double mv, double mx, double kappa, double alphaD){
    using namespace std::placeholders;
    auto dbr = std::bind(dbrpi0_to_gamma_dm_dm_baryonic,mv,mx,kappa,alphaD,_1);
    if(mv>2*mx&&mv<mpion){
        double up = mv*mv+5*mv*GammaVB(mv, mx, kappa, alphaD);
        double down = mv*mv-5*mv*GammaVB(mv, mx, kappa, alphaD);
        if(down<4*mx*mx)
            down = 4*mx*mx;

        if(up>pow(mpion,2))
            up = pow(mpion,2);

        double low = DoubleExponential_adapt(dbr,4*mx*mx,down, 2 ,0.1,1e-4);
        double mid = DoubleExponential_adapt(dbr,down, up, 2 ,0.1,1e-5);
        double high = DoubleExponential_adapt(dbr, up, pow(mpion,2), 2 ,0.1,1e-4);
        return low+mid+high;
    }
    else{
        return DoubleExponential_adapt(dbr,4*mx*mx,pow(mpion,2),100,0.1,1e-4);
    }
}

double bretatoVBgamma(double mv, double mx, double kappa, double alphaD){
    if(mv>meta)
        return 0;
    return 2*pow(etafactor*sqrt(alphaD/alphaem)-kappa,2)*pow(1-mv*mv/pow(meta,2),3)*bretato2gamma;
}

double d2breta_to_gamma_dm_dm_baryonic(double mv, double mx, double kappa, double alphaD, double s, double theta){
    return bretato2gamma*sqrt(1-4*mx*mx/s)*pow(pow(meta,2)-s,3)*(s-4*mx*mx)*pow(alphaD,1)*pow(etafactor*sqrt(alphaD/alphaem)-kappa,2)*pow(sin(theta),3)/(8*pow(meta,6)*pi*(pow(s-mv*mv,2)+pow(mv*GammaVB(mv,mx,kappa,alphaD),2)));
}

double dbreta_to_gamma_dm_dm_baryonic(double mv, double mx, double kappa, double alphaD, double s){
    return bretato2gamma*sqrt(1.0-4.0*mx*mx/s)*pow(pow(meta,2)-s,3)*(s-4.0*mx*mx)*pow(alphaD,1)*pow(etafactor*sqrt(alphaD/alphaem)-kappa,2)/(6.0*pow(meta,6)*pi*(pow(s-mv*mv,2)+pow(mv*GammaVB(mv,mx,kappa,alphaD),2)));
}

double breta_to_gamma_dm_dm_baryonic(double mv, double mx, double kappa, double alphaD){
    using namespace std::placeholders;
    auto dbr = std::bind(dbreta_to_gamma_dm_dm_baryonic,mv,mx,kappa,alphaD,_1);
    if(mv>2*mx&&mv<meta){
        double up = mv*mv+5*mv*GammaVB(mv, mx, kappa, alphaD);
        double down = mv*mv-5*mv*GammaVB(mv, mx, kappa, alphaD);
        if(down<4*mx*mx)
            down = 4*mx*mx;

        if(up>pow(meta,2))
            up = pow(meta,2);
        
        double low = DoubleExponential_adapt(dbr,4*mx*mx,down, 100 ,0.1,1e-4);
        double mid = DoubleExponential_adapt(dbr,down, up, 100 ,0.1,1e-4);
        double high = DoubleExponential_adapt(dbr, up, pow(meta,2), 100 ,0.1,1e-4);
        return low+mid+high;
    }
    else{
        return DoubleExponential_adapt(dbr,4*mx*mx,pow(meta,2),100,0.1,1e-4);
    }
}

double bromega_to_Vb(double mv, double mx, double kappa, double alphaD){

    return 0.25*pow(2*sqrt(alphaD/alphaem)-kappa,2)*alphaD/alphaem*Br_omega_to_e_e*pow(momega,4)/ \
        (pow(pow(momega,2)-pow(mv,2),2)+1e-4*pow(momega,4))*pow(1-4*pow(mx/mv,2),1.5);//Where does that 1e-4 come from? Must be the width of Omega.
}

double brphi_to_Vb(double mv, double mx, double kappa, double alphaD){
	if(mv-mphi>0.25||mphi-mv>0.25)
        return 0.0;
    return 0.25*pow(-sqrt(alphaD/alphaem)-kappa,2)*alphaD/alphaem*Br_phi_to_e_e*pow(momega,4)/ \
        (pow(pow(mphi,2)-pow(mv,2),2)+1e-4*pow(mphi,4))*pow(1-4*pow(mx/mv,2),1.5);
}
