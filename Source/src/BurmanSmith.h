#ifndef GUARD_BURMANSMITH_H
#define GUARD_BURMANSMITH_H

//pi^+ distribution from Burman and Smith.
//LA-11502-MS
//www.osti.gov/scitech/servlets/purl/6167579

#include <cmath>
#include <functional>
#include <iostream>
#include "Random.h"
#include "Distribution.h"

double Bspline(int, int, double);

class BurmanSmith : public Distribution {
    public:
        void sample_particle(Particle &);
        void sample_momentum(double &,double &, double &);
        BurmanSmith(double beam_kinetic_energy, int proton_number);
        double fpimax;
        double dsigma(double Tpi, double theta);
   private:
        double Tp, Ep, NormZ, TA, TAZ730, TAZ585, sigmaAZ585, sigmaAZ730;
        double sigmaA;
        double Tpimax;
        int Z;
        double TF(double);
        double Amp(double);
        double Tbar(double);
        double sigma(double);
        double ampz[5];
        double Z1coefficient;
};

#endif
