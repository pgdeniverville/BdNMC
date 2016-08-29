#ifndef GUARD_DMscattering_h
#define GUARD_DMscattering_h

#include <memory>

#include "Particle.h"
#include "detector.h"

double ThetaEe (double,double,double);
double EeTheta (double,double,double);
double EeTMax (double,double);
double EeTMin (double,double);
double F1 (double,double,double,double);
double dsigmadEe (double,double,double,double,double,double);
double dsigmadEe_scaled (double,double,double,double,double,double);
double F2 (double,double,double,double);
double sigma (double,double,double,double,double);
double sigma2 (double,double,double,double,double,double,double);
//double quickscatter(double EDM, double MDM, double MDP, double kap, double alD);
//double quickscatter2(double EDM, double MDM, double MDP, double kap, double alD, double Emax, double Emin);
//void scatterevent (double MDP, double kap, double alD, Particle DM, Particle &electron);
//bool probscatterelectron(double &,double,double,double, const std::shared_ptr<detector>& det,Particle DM);
//bool probscatterelectron(double &,double,double,double, const std::shared_ptr<detector>& det,Particle DM, Particle &electron);

#endif
