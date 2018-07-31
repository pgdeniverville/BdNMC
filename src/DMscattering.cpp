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
    //cout << Ee << " " << EDM << " " << MDM << " " << MDP << " " << kappa << " " << alphaD << endl;
	//cout << 4*Pi*kappa*kappa*alphaEM*alphaD*(2*EDM*Me*(-Ee+EDM+Me)+MDM*MDM*(Me-Ee))/(EDM-MDM)/(EDM+MDM)/pow(2*Me*(Ee-Me)+MDP*MDP,2) << endl;
	return(4*Pi*kappa*kappa*alphaEM*alphaD*(2*EDM*Me*(-Ee+EDM+Me)+MDM*MDM*(Me-Ee))/(EDM-MDM)/(EDM+MDM)/pow(2*Me*(Ee-Me)+MDP*MDP,2));
}

double dsigmadEe_scaled (double Ee, double EDM, double MDM, double MDP, double kappa, double alphaD) {
    return dsigmadEe (Ee, EDM, MDM, MDP, kappa, alphaD)*convGeV2cm2*convmcm;
}

// Function F2(Ee) 
// Total DM - electron scattering cross section equals
// sigma =  4*Pi*kappa*kappa*alpha*alphaD*( F2(EeMax)- F2(EeMin) ) 
double F2 (double Ee, double EDM, double MDM, double MDP) {
	return( -((4*EDM*EDM*Me*Me+2*EDM*Me*MDP*MDP+MDP*MDP*MDM*MDM)/(2*Me*(Ee-Me)+MDP*MDP)+(2*EDM*Me+MDM*MDM)*log(2*Me*(Ee-Me)+MDP*MDP))/(Me*Me*(EDM-MDM)*(EDM+MDM)) );
}
// DM - electron scattering total cross section sigma
double sigma (double EDM, double MDM, double MDP, double kappa, double alphaD) {
	double rsig;
	double coef;
	coef = Pi*kappa*kappa*alphaEM*alphaD;
	double EeMaxA, EeMinA;
	EeMaxA = EeTMax(EDM,MDM);
	EeMinA = EeTMin(EDM,MDM);
    cout << coef*F2(EeMaxA,EDM,MDM,MDP) << " " << coef*F2(EeMinA,EDM,MDM,MDP) << endl;
    rsig = coef*(F2(EeMaxA,EDM,MDM,MDP)-F2(EeMinA,EDM,MDM,MDP));
	return(rsig);
}

double sigma2 (double EDM, double MDM, double MDP, double kappa, double alphaD, double Emax, double Emin) {
	if(Emax<Emin)
		return 0;
    //cout << EDM << " " << MDM << " " << MDP << " " << kappa << " " << alphaD << " " << Emax << " " << Emin << endl;
	//cout << Pi*kappa*kappa*alphaEM*alphaD*(F2(Emax,EDM,MDM,MDP)) - Pi*kappa*kappa*alphaEM*alphaD*F2(Emin,EDM,MDM,MDP) << endl;
    return(Pi*kappa*kappa*alphaEM*alphaD*(F2(Emax,EDM,MDM,MDP)-F2(Emin,EDM,MDM,MDP)));
}
