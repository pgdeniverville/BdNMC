//#include "Event.h"
#include "DMNscattering.h"

#include <iostream>
#include "Kinematics.h"

const double pi = 3.14159265359;
const double alphaEM = 1.0/137.035999074;
//All masses are defined in terms of GeV.
const double mp = 0.938272;
const double mn = 0.939565;

//unit conversions
const double convGeV2cm2 = 3.89e-28;
const double convmcm = 100.0;

//Form factor fit parameters from hep-ex/0308005.
const double aEN = 0.942; const double bEN = 4.61; const double mup = 2.793;
const double mun = -1.913; const double MV2 = 0.71;
const double aEP [7] = { 1 , 3.253 , 1.422 , 0.08582 , 0.3318 , -0.09371, 0.01076};
const double aMP [7] = { 1 , 3.104 , 1.428 , 0.1112 , -0.006981 , 0.0003705, -0.7063 * pow(10,-5)};
const double aMN [7] = { 1, 3.043, 0.8548 , 0.6806, -0.1287, 0.008912, 0};
const double GEP0 = 1; const double GMP0 = mup; const double GMN0 = mun;

const double suppression_constant = 1.0;

double A(double, double, double, const double);
double B(double, double, double, const double);
double C(double, double ,double, const double);
double F1P(double);
double F2P(double);
double F1N(double);
double F2N(double);
double GEPact(double);
double GEP03(double);
double GEN03(double);
double GMP03(double);
double GMN03(double);
double CoherentFormFactor(double q);


using std::cout; using std::endl;

/*
 * Returns the effective charge radius of a nucleon with number of nucleons A. 
 */

const double R_a=0.52;
const double R_s=0.9;

double RadiusFunction(double A){
    return sqrt(pow(1.23*pow(A,1.0/3.0)-0.6,2)+7.0/3.0*pow(pi*R_a,2)-5*pow(R_s,2));
}

double BesselJ1(double x){
    return (sin(x)/x-cos(x))/x;
}

double CoherentFormFactor(double q, double A){
    if(q==0)
        return 1.0;
    return 3.0/(q*RadiusFunction(A))*BesselJ1(q*RadiusFunction(A))*exp(-suppression_constant*q*q/2.0);
}
/* Returns the differential scattering cross section between a proton and 
 * a dark matter particle for initial dark matter energy E and final 
 * dark matter energy Edm. All energies in GeV.
 * I may need to implement a more general equation that contains nucleon angle. 
*/



double dsigmadEdmP(double E, double Edm,  double mdm, double mV, double alphaprime, double kappa)
{
    return alphaprime * alphaEM * pow(kappa,2) * 4 * pi * \
        (pow(F1P(2*mp*(E-Edm)),2) * A(E, Edm, mdm, mp) + \
        F1P(2*mp*(E-Edm)) * F2P(2*mp*(E-Edm)) * C(E, Edm, mdm, mp) - \
        1.0/4.0 * pow(F2P(2*mp*(E-Edm)),2) * B(E, Edm, mdm, mp)) \
        /((pow((pow(mV,2) + 2*mp * (E - Edm)),2))*(E*E - pow(mdm, 2)));
}

double dsigmadEdmP_coherent(double E, double Edm,  double mdm, double mV, double alphaprime, double kappa, double A, double Z)
{
    double mass = (A-Z)*mn+Z*mp;
    return pow(Z,2)*pow(CoherentFormFactor(sqrt(2*mass*(E-Edm)),A),2)*dsigmadEdmP(E, Edm, mdm, mV, alphaprime, kappa);
}

double dsigmadEdmN(double E, double Edm,  double mdm, double mV, double alphaprime, double kappa)
{
    return alphaprime * alphaEM * pow(kappa,2) * 4 * pi * \
        (pow(F1N(2*mn*(E-Edm)),2) * A(E, Edm, mdm, mn) + \
        F1N(2*mn*(E-Edm)) * F2N(2*mn*(E-Edm)) * C(E, Edm, mdm, mn) - \
        1.0/4.0 * pow(F2N(2*mn*(E-Edm)),2) * B(E, Edm, mdm, mn)) \
        /((pow((pow(mV,2) + 2*mn * (E - Edm)),2))*(E*E - pow(mdm, 2)));
}

/*Coefficients in the scattering cross section.*/

double A(double E, double Edm, double mdm, const double mN){
    return 2*mN*E*Edm-pow(mdm, 2)*(E-Edm);
}

double B(double E, double Edm, double mdm, const double mN){
    return (Edm-E)*(pow((Edm+E),2)+2*mN*(Edm-E)-4*pow(mdm,2));
}

double C(double E, double Edm, double mdm, const double mN){
    return (Edm - E) * (mN*(E-Edm)+2*pow(mdm,2));
}

/*BBA-03 form factors for neutrino-nucleon scattering from hep-ex/0308005.*/

double F1P(double q2){
   return (GEP03(q2) + q2/(4*pow(mp,2)) * GMP03(q2))/(1+q2/(4 * pow(mp,2) )); 
}

double F2P(double q2){
    return (GMP03(q2)-GEP03(q2))/(1+q2/(4*pow(mp,2)));
}

double F1N(double q2){
    return (GEN03(q2) + q2/(4*pow(mn,2)) * GMN03(q2))/(1+q2/(4 * pow(mn,2) )); 
}

double F2N(double q2){
    return (GMN03(q2)-GEN03(q2))/(1+q2/(4*pow(mn,2)));
}

double GEPact(double q2){
    return GEP0/(1 + aEP[1]*q2 + aEP[2]*pow(q2,2) + aEP[3]*pow(q2,3)+ aEP[4]*pow(q2,4)+\
        aEP[5]*pow(q2,5)+ aEP[6]*pow(q2,6));
}
//There is a change in the fit at q^2 = 6 GeV^2.
double GEP03(double q2){
    if(q2<6){
        return GEPact(q2);
    } else{
        return GMP03(q2)*GEPact(6)/GMP03(6); 
    }
}

double GMP03(double q2){
    return GMP0/(1 + aMP[1]*q2 + aMP[2]*pow(q2,2) + aMP[3]*pow(q2,3)+ aMP[4]*pow(q2,4)+\
        aMP[5]*pow(q2,5)+ aMP[6]*pow(q2,6));
}

double GEN03(double q2){
    return -mun * aEN * q2/(4*mn*mn) / (1 + bEN*q2/(4 * mn*mn)) / pow((1+ q2/MV2),2);
}

double GMN03(double q2){
    return GMN0/(1 + aMN[1]*q2 + aMN[2]*pow(q2,2) + aMN[3]*pow(q2,3)+ aMN[4]*pow(q2,4)+\
        aMN[5]*pow(q2,5)+ aMN[6]*pow(q2,6));
}

double nucleon_recoil_energy(double Edm, double mdm, double mN, double theta){
	return mN*(3*Edm*Edm+4*Edm*mN+2*mN*mN-mdm*mdm+(Edm*Edm-mdm*mdm)*cos(2*theta))/(Edm*Edm+4*Edm*mN+2*mN*mN+mdm*mdm+(mdm*mdm-Edm*Edm)*cos(2*theta));
}

/*double Efmin(double En, double mdm, double m)
{
    return En + m - m*(pow(En+m,2)+(pow(En,2)-pow(mdm,2)))/(pow(En+m,2)-(pow(En,2)-pow(mdm,2)));
}*/

double Efmin(double Edm, double mdm, double mN){
	return (2*mN*mdm*mdm+Edm*(mN*mN+mdm*mdm))/(2*Edm*mN+mN*mN+mdm*mdm);	
}

//Calculates the angle of the outgoing nucleon, relative to the incoming dark matter trajectory.
//I should find a way to replace this with atan2
double Ef_to_N_Theta(double En, double Ef, double mdm, double m){
    return acos((En*(En+m-Ef)-m*Ef)/sqrt(En*En-mdm*mdm)/sqrt(pow(En+m-Ef,2)-m*m));
}

/*
bool probscatter (double &pMax, double mv, double kap, double alD, detector det, Particle DM, Particle Nucleon, Linear_Interpolation LIN, Linear_Interpolation LIP){
    
    using namespace std::placeholders;

    double LXDet = convmcm*det.Ldet(DM);
    double XDMp = LIP.Interpolate(DM.E)*npCH2;
    double XDMn = LIN.Interpolate(DM.E)*nnCH2;
    double prob;
    if((prob = LXDet*convGeV2cm2*(XDMp+XDMn)) > pMax)
        pMax = prob;
    if(prob > pMax*Random::Flat(0,1)){
        if(XDMp/(XDMn+XDMp) <= Random::Flat(0,1)){
            auto fn = std::bind(dsigmadEdmP,DM.E,_1,DM.m,mv,alD,kap);//this might not work.
            std::function<double(double)> Xsec = fn;
            Nucleon.m = mp;
            Nucleon.name = "proton";
            scatterevent(DM, Nucleon, Xsec);
        }
        else{
            auto fn = std::bind(dsigmadEdmN,DM.E,_1,DM.m,mv,alD,kap);
            std::function <double(double)> Xsec = fn;
            Nucleon.m = mn;
            Nucleon.name = "neutron";
            scatterevent(DM, Nucleon, Xsec);
        }
        Event coordinates(0);
        coordinates.FourPosition(DM, det);
        Nucleon.coords = & coordinates;
        cout << Nucleon.name << " " << Nucleon.px << " " << Nucleon.coords->x << endl;;
        Nucleon.EVENT_SET = true;
        return true;
    }
    else
        return false;
}*/
/*
bool probscatter_coherent (double &pMax, const double mv, const double kap, const double alD, const std::shared_ptr<detector>& det, const int mat_index, Particle &DM, Linear_Interpolation &LIP){

    using namespace std::placeholders;

    double LXDet = convmcm*det->Ldet(DM);
    double XDMc = LIP.Interpolate(DM.E)*(det->get_nDensity(mat_index));
    double prob;
    if((prob = LXDet*convGeV2cm2*(XDMc)) > pMax*Random::Flat(0,1)){
        if(prob>pMax)
            pMax = prob;
        return true;
    }
    else
        return false;
}

bool probscatter_coherent (double &pMax, const double mv, const double kap, const double alD,  const std::shared_ptr<detector>& det, const int mat_index, Particle &DM, Particle *Nucleus, Linear_Interpolation &LIP){

    using namespace std::placeholders;
    double LXDet = convmcm*det->Ldet(DM);
    double XDMc = LIP.Interpolate(DM.E)*(det->get_nDensity(mat_index));
    double prob = LXDet*convGeV2cm2*(XDMc);
    if(prob > pMax*Random::Flat(0,1)){
        if(prob>pMax)
            pMax = prob;
        auto fn = std::bind(dsigmadEdmP,DM.E,_1,DM.m,mv,alD,kap);//this might not work.
        std::function<double(double)> Xsec = fn;
        scatterevent(DM, *Nucleus, Xsec);
        return true;
    }
    else
        return false;
}*/
