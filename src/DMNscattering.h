#ifndef GUARD_DMNscattering_h
#define GUARD_DMNscattering_h

#include <cmath>
#include <memory>

double dsigmadEdmP(double E, double Edm, double mdm, double mV, double alphaprime, double kappa);
double dsigmadEdmP_coherent(double E, double Edm, double mdm, double mV, double alphaprime, double kappa, double A, double Z);
double dsigmadEdmN(double E, double Edm, double mdm, double mV, double alphaprime, double kappa);

//Kinematic limits
double Efmin(double E, double mdm, double m);
double Ef_to_N_Theta(double E, double Ef, double mdm, double m);

//Form factors
double CoherentFormFactor(double q, double A);

double RadiusFunction(double);

#endif
