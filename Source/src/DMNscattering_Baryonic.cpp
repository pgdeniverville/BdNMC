//#include "Event.h"
#include "DMNscattering_Baryonic.h"

#include <iostream>
#include "Kinematics.h"
#include "Random.h"
#include "constants.h"

//All masses are defined in terms of GeV.
const double mp = MASS_PROTON;
const double mn = MASS_NEUTRON;

//unit conversions
//const double convmcm = 100.0;

//Form factor fit parameters from hep-ex/0308005.
const double aEN = 0.942; const double bEN = 4.61; const double mup = 2.793;
const double mun = -1.913; const double MV2 = 0.71;
const double aEP [7] = { 1 , 3.253 , 1.422 , 0.08582 , 0.3318 , -0.09371, 0.01076};
const double aMP [7] = { 1 , 3.104 , 1.428 , 0.1112 , -0.006981 , 0.0003705, -0.7063 * pow(10,-5)};
const double aMN [7] = { 1, 3.043, 0.8548 , 0.6806, -0.1287, 0.008912, 0};
const double GEP0 = 1; const double GMP0 = mup; const double GMN0 = mun;

const double suppression_constant = 0.9;


using std::cout; using std::endl;

/*
 * Returns the effective charge radius of a nucleon with number of nucleons A. 
 */

const double R_a=0.52;
const double R_s=0.9;



double DMNscattering_Baryonic::RadiusFunction(double A){
    return sqrt(pow(1.23*pow(A,1.0/3.0)-0.6,2)+7.0/3.0*pow(pi*R_a,2)-5*pow(R_s,2));
}

double DMNscattering_Baryonic::BesselJ1(double x){
    return (sin(x)/x-cos(x))/x;
}

double DMNscattering_Baryonic::CoherentFormFactor(double q, double A){
    if(q==0)
        return 1.0;
    return 3.0/(q*RadiusFunction(A))*BesselJ1(q*RadiusFunction(A))*exp(-suppression_constant*q);
}


/*Coefficients in the scattering cross section.*/

double DMNscattering_Baryonic::A(double E, double Edm, double mdm, const double mN){
    return 2*mN*E*Edm-pow(mdm, 2)*(E-Edm);
}

double DMNscattering_Baryonic::B(double E, double Edm, double mdm, const double mN){
    return (Edm-E)*(pow((Edm+E),2)+2*mN*(Edm-E)-4*pow(mdm,2));
}

double DMNscattering_Baryonic::C(double E, double Edm, double mdm, const double mN){
    return (Edm - E) * (mN*(E-Edm)+2*pow(mdm,2));
}

/*BBA-03 form factors for neutrino-nucleon scattering from hep-ex/0308005.*/
//Cut off at q2>10 GeV^2. These form factors are not well callibrated beyond this point.
double DMNscattering_Baryonic::F1P(double q2){
   if(q2>10)
	   return 0;
   else
		return (GEP03(q2) + q2/(4*pow(mp,2)) * GMP03(q2))/(1+q2/(4 * pow(mp,2) )); 
}

double DMNscattering_Baryonic::F2P(double q2){
    if(q2>10)
	   return 0;
   else
	return (GMP03(q2)-GEP03(q2))/(1+q2/(4*pow(mp,2)));
}

double DMNscattering_Baryonic::F1N(double q2){
    if(q2>10)
	   return 0;
   else
	return (GEN03(q2) + q2/(4*pow(mn,2)) * GMN03(q2))/(1+q2/(4 * pow(mn,2) )); 
}

double DMNscattering_Baryonic::F2N(double q2){
    if(q2>10)
	   return 0;
   else
	return (GMN03(q2)-GEN03(q2))/(1+q2/(4*pow(mn,2)));
}

double DMNscattering_Baryonic::Fs1(double q2){
    return F1P(q2)+F1N(q2);
}

double DMNscattering_Baryonic::Fs2(double q2){
    return F2P(q2)+F2N(q2);
}

double DMNscattering_Baryonic::Fv1(double q2){
    return F1P(q2)-F1N(q2);
}

double DMNscattering_Baryonic::Fv2(double q2){
    return F2P(q2)-F2N(q2);
}

double DMNscattering_Baryonic::FS1(double q2){
    return 0;
}

double DMNscattering_Baryonic::FS2(double q2){
   return 0;
}

double DMNscattering_Baryonic::gu(double alD, double kappa){
    return sqrt(4*pi*alD)/3-2*kappa*sqrt(4*pi*alphaEM)/3;
}

double DMNscattering_Baryonic::gd(double alD, double kappa){
    return sqrt(4*pi*alD)/3+kappa*sqrt(4*pi*alphaEM)/3;
}

double DMNscattering_Baryonic::F1Pb(double q2, double alD, double kappa){
   return  0.5*(gu(alD, kappa)-gd(alD, kappa))*Fv1(q2) + 1.5*(gu(alD, kappa)+gd(alD, kappa))*Fs1(q2)+(gu(alD, kappa)+2*gd(alD, kappa))*FS1(q2);
}

double DMNscattering_Baryonic::F2Pb(double q2, double alD, double kappa){
   return  0.5*(gu(alD, kappa)-gd(alD, kappa))*Fv2(q2) + 1.5*(gu(alD, kappa)+gd(alD, kappa))*Fs2(q2)+(gu(alD, kappa)+2*gd(alD, kappa))*FS2(q2);
}

double DMNscattering_Baryonic::F1Nb(double q2, double alD, double kappa){
   return  -0.5*(gu(alD, kappa)-gd(alD, kappa))*Fv1(q2) + 1.5*(gu(alD, kappa)+gd(alD, kappa))*Fs1(q2)+(gu(alD, kappa)+2*gd(alD, kappa))*FS1(q2);
}

double DMNscattering_Baryonic::F2Nb(double q2, double alD, double kappa){
   return  -0.5*(gu(alD, kappa)-gd(alD, kappa))*Fv2(q2) + 1.5*(gu(alD, kappa)+gd(alD, kappa))*Fs2(q2)+(gu(alD, kappa)+2*gd(alD, kappa))*FS2(q2);
}

double DMNscattering_Baryonic::tau(double q2, double m){
    return q2/(4*m*m);
}

double DMNscattering_Baryonic::GEPact(double q2){
    return GEP0/(1 + aEP[1]*q2 + aEP[2]*pow(q2,2) + aEP[3]*pow(q2,3)+ aEP[4]*pow(q2,4)+\
        aEP[5]*pow(q2,5)+ aEP[6]*pow(q2,6));
}
//There is a change in the fit at q^2 = 6 GeV^2.
double DMNscattering_Baryonic::GEP03(double q2){
    if(q2<6.0){
        return GEPact(q2);
    } else{
        return GMP03(q2)*GEPact(6)/GMP03(6); 
    }
}

double DMNscattering_Baryonic::GMP03(double q2){
    return GMP0/(1 + aMP[1]*q2 + aMP[2]*pow(q2,2) + aMP[3]*pow(q2,3)+ aMP[4]*pow(q2,4)+\
        aMP[5]*pow(q2,5)+ aMP[6]*pow(q2,6));
}

double DMNscattering_Baryonic::GEN03(double q2){
    return -mun * aEN * q2/(4.0*mn*mn) / (1.0 + bEN*q2/(4.0 * mn*mn)) / pow((1.0+ q2/MV2),2);
}

double DMNscattering_Baryonic::GMN03(double q2){
    return GMN0/(1 + aMN[1]*q2 + aMN[2]*pow(q2,2) + aMN[3]*pow(q2,3)+ aMN[4]*pow(q2,4)+\
        aMN[5]*pow(q2,5)+ aMN[6]*pow(q2,6));
}


/* Returns the differential scattering cross section between a proton and 
 * a dark matter particle for initial dark matter energy E and final 
 * dark matter energy Edm. All energies in GeV.
 * I may need to implement a more general equation that contains nucleon angle. 
*/
double DMNscattering_Baryonic::dsigmadEdmP(double E, double Edm,  double mdm, double mV, double alD, double kappa)
{
    return alD * (pow(F1Pb(2*mp*(E-Edm),alD,kappa),2) * A(E, Edm, mdm, mp) + \
        F1Pb(2*mp*(E-Edm),alD,kappa) * F2Pb(2*mp*(E-Edm),alD,kappa) * C(E, Edm, mdm, mp) - \
        1.0/4.0 * pow(F2Pb(2*mp*(E-Edm),alD,kappa),2) * B(E, Edm, mdm, mp)) \
        /((pow((pow(mV,2) + 2*mp * (E - Edm)),2))*(E*E - pow(mdm, 2)));
}

double DMNscattering_Baryonic::dsigmadEdmP_coherent(double E, double Edm,  double mdm, double mV, double alphaprime, double kappa, double A, double Z)
{
    double mass = (A-Z)*mn+Z*mp;
    return pow(A,2)*pow(CoherentFormFactor(sqrt(2*mass*(E-Edm))/GeVtofm,A),2)*dsigmadEdmP(E, Edm, mdm, mV, alphaprime, kappa);
}

double DMNscattering_Baryonic::dsigmadEdmN(double E, double Edm,  double mdm, double mV, double alD, double kappa)
{
    return alD * (pow(F1Nb(2*mn*(E-Edm),alD,kappa),2) * A(E, Edm, mdm, mn) + \
        F1Nb(2*mn*(E-Edm),alD,kappa) * F2Nb(2*mn*(E-Edm),alD,kappa) * C(E, Edm, mdm, mn) - \
        1.0/4.0 * pow(F2Nb(2*mn*(E-Edm),alD,kappa),2) * B(E, Edm, mdm, mn)) \
        /((pow((pow(mV,2) + 2*mn * (E - Edm)),2))*(E*E - pow(mdm, 2)));
}

double DMNscattering_Baryonic::Efmin(double En, double mdm, double m)
{
    return En + m - m*(pow(En+m,2)+(pow(En,2)-pow(mdm,2)))/(pow(En+m,2)-(pow(En,2)-pow(mdm,2)));
}

//Calculates the angle of the outgoing nucleon, relative to the incoming dark matter trajectory.
double DMNscattering_Baryonic::Ef_to_N_Theta(double En, double Ef, double mdm, double m){
    return acos((En*(En+m-Ef)-m*Ef)/sqrt(En*En-mdm*mdm)/sqrt(pow(En+m-Ef,2)-m*m));
}
