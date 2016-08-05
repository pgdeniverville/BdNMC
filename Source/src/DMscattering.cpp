#include "DMscattering.h"
#include <cmath>
#include <stdlib.h>
#include <iostream>
#include "Random.h"
#include "Kinematics.h"
#include "constants.h"

using std::cout; using std::endl;

const double Me = MASS_ELECTRON;
const double Pi = pi;
const double convmcm = 100.0;

// Electron angle as a function of 
// electron energy and dark matter energy and mass  
double ThetaEe (double Ee, double EDM, double MDM) {
	return(acos(sqrt((Ee-Me)/(Ee+Me))*(EDM+Me)/sqrt(EDM*EDM-MDM*MDM)));
}
// Electron energy as a function of 
// electron scattering angle and dark matter energy
double EeTheta (double EDM, double Thetael, double MDM) {
	double rEeTheta;
	double EeThetaN, EeThetaD;
	EeThetaN = Me*((EDM+Me)*(EDM+Me)+(EDM*EDM-MDM*MDM)*cos(Thetael)*cos(Thetael));
	EeThetaD = (EDM+Me)*(EDM+Me)-(EDM*EDM-MDM*MDM)*cos(Thetael)*cos(Thetael);
	rEeTheta = EeThetaN/EeThetaD;
	return(rEeTheta);
}
// Maximum electron energy as a function of 
// dark matter energy and mass
double EeTMax (double EDM, double MDM) {
	double rEeTMax;
	double ThetaelMax = 0.0;
	rEeTMax = EeTheta(EDM,ThetaelMax,MDM);
	return(rEeTMax);
}
// Minimum electron energy as a function of 
// dark matter energy and mass
double EeTMin (double EDM, double MDM) {
	double rEeTMin;
	double ThetaelMin = Pi/2.0;
	rEeTMin = EeTheta(EDM,ThetaelMin,MDM);
	return(rEeTMin);
}
// Function F1 
// dsigma/dEe = 4*Pi*kappa*kappa*alpha*alphaD*F1
double F1 (double Ee, double EDM, double MDM, double MDP) {
	double rF1;
	double F1N, F1D;
	F1N = 2.0*Me*EDM*EDM-(2.0*Me*EDM+MDM*MDM)*(Ee-Me);
	F1D = (EDM*EDM-MDM*MDM)*pow(MDP*MDP+2*Me*Ee-2*Me*Me,2);
	rF1 = F1N/F1D;
	return(rF1);
}
//  differential DM - electron scattering cross section dsigma/dEe
double dsigmadEe (double Ee, double EDM, double MDM, double MDP, double kappa, double alphaD) {
	double rdsig;
	double coef;
	coef = 4*Pi*kappa*kappa*alphaEM*alphaD;
	rdsig = coef*F1(Ee,EDM,MDM,MDP);
	return(rdsig);
}

double dsigmadEe_scaled (double Ee, double EDM, double MDM, double MDP, double kappa, double alphaD) {
    return dsigmadEe (Ee, EDM, MDM, MDP, kappa, alphaD)*convGeV2cm2*convmcm;
}

// Function F2(Ee) 
// Total DM - electron scattering cross section equals
// sigma =  4*Pi*kappa*kappa*alpha*alphaD*( F2(EeMax)- F2(EeMin) ) 
double F2 (double Ee, double EDM, double MDM, double MDP) {
	double rF2;
	double F2N1, F2N2, F2D;
	F2N1 = (4*EDM*EDM*Me*Me+2*EDM*Me*MDP*MDP + MDP*MDP*MDM*MDM)/(2*Ee*Me-2*Me*Me+MDP*MDP);
	F2N2 = (2*EDM*Me+MDM*MDM)*log(2*Ee*Me-2*Me*Me+MDP*MDP);
	F2D  = 4*Me*Me*(EDM*EDM-MDM*MDM);
	rF2  = -(F2N1+F2N2)/F2D;
	return(rF2);
}
// DM - electron scattering total cross section sigma
double sigma (double EDM, double MDM, double MDP, double kappa, double alphaD) {
	double rsig;
	double coef;
	coef = 4*Pi*kappa*kappa*alphaEM*alphaD;
	double EeMaxA, EeMinA;
	EeMaxA = EeTMax(EDM,MDM);
	EeMinA = EeTMin(EDM,MDM);
	rsig = coef*(F2(EeMaxA,EDM,MDM,MDP)-F2(EeMinA,EDM,MDM,MDP));
	return(rsig);
}


double sigma2 (double EDM, double MDM, double MDP, double kappa, double alphaD, double Emax, double Emin) {
	if(Emax<Emin)
		return 0;
	return(4*Pi*kappa*kappa*alphaEM*alphaD*(F2(Emax,EDM,MDM,MDP)-F2(Emin,EDM,MDM,MDP)));
}
