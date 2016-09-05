#ifndef GUARD_DMNscattering_Baryonic_h
#define GUARD_DMNscattering_Baryonic_h

#include "Particle.h"
#include "Integrator.h"
#include "detector.h"

class DMNscattering_Baryonic{
public:
	static double dsigmadEdmP(double E, double Edm,  double mdm, double mV, double alD, double kappa);
	static double dsigmadEdmP_coherent(double E, double Edm,  double mdm, double mV, double alphaprime, double kappa, double A, double Z);
	static double dsigmadEdmN(double E, double Edm,  double mdm, double mV, double alD, double kappa);
	//Kinematic Limits
	static double Efmin(double En, double mdm, double m);
	static double Ef_to_N_Theta(double En, double Ef, double mdm, double m);
	//void TEST();
	static double RadiusFunction(double A);
	static double CoherentFormFactor(double q, double A);
private:
	static double A(double, double, double, const double);
	static double B(double, double, double, const double);
	static double C(double, double ,double, const double);
	static double GEPact(double);
	static double GEP03(double);
	static double GEN03(double);
	static double GMP03(double);
	static double GMN03(double);
	static double F1P(double);
	static double F2P(double);
	static double F1N(double);
	static double F2N(double);
	static double Fs1(double);	
	static double FS1(double);	
	static double FS2(double);	
	static double Fs2(double);	
	static double Fv1(double);	
	static double Fv2(double);	
	static double F1Pb(double,double,double);
	static double F2Pb(double,double,double);
	static double F1Nb(double,double,double);
	static double F2Nb(double,double,double);
	static double CoherentFormFactor(double q);
	static double gu(double, double);
	static double gd(double, double);
	static double tau(double, double);
	static double BesselJ1(double);
};

#endif
