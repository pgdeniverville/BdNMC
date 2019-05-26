#ifndef GUARD_decay_h
#define GUARD_decay_h

#include "Particle.h"
#include "Random.h"
#include <memory>

/*
 * decay.h and decay.cpp supplies a number of useful functions for calculating
 * decay kinematics.
 */


void DecayDP (Particle &, Particle &);
void DecayDM (Particle &, Particle &, Particle &, Particle &);
void DecayDM_Off_Shell(Particle &daughter1, Particle &daughter2, Particle &mediato, Particle &parent, double theta);
void TwoBodyDecay (Particle &parent, Particle &daughter1, Particle &daughter2);
void TwoBodyDecay (Particle &parent, Particle &daughter1, Particle &daughter2, double thetad);
void TwoBodyDecay (Particle &parent, Particle &daughter1, Particle &daughter2, double thetad, double phid);
void Meson_Capture_Off_Shell(Particle &daughter1, Particle &daughter2, Particle &mediator, double theta);

namespace Three_Body_Decay_Space{
    void Three_Body_Decay(Particle &parent, Particle &daughter1, Particle &daughter2, Particle &daughter3, double &d_width_max, std::function<double(double,double,double,double,double,double)> &amplitude);
    
    double E2star(double m12s, double m1, double m2);

    //m0 is the mass of the parent particle.
    double E3star(double m12s, double m0, double m3);

    double p2star(double m12s, double m1, double m2);

    double p3star(double m12s, double m0, double m3);

    double m23min(double m12s, double m0, double m1, double m2, double m3);

    double m23max(double m12s, double m0, double m1, double m2, double m3);
    //!!!!!Note that this is using m12, not m12s=m12^2!!!!
    double p1star(double m12, double m1, double m2);

    double p3rest(double m12, double m0, double m3);

    //p2.p3 dot product evaluated in p1+p2 rest frame!
    double p2p3(double m12s, double CosTheta, double m0, double m1, double m2, double m3);

    double Cos_Theta_to_m23s(double m12s, double CosTheta, double m0, double m1, double m2, double m3);

    //Eq 47.22 in the PDG.
    //Function expects amp(m12^2, m23^2, m0, m1, m2, m3)
    double d_decay_width(std::function<double(double,double,double,double,double,double)> amp, double m12s, double m23s, double m0, double m1, double m2, double m3);

    //This includes the angular dependence
    //EQ 47.20 in the PDG
    //amp(m12s,cos_t,m0, m1,m2,m3), where cos_t is the angle between p1 and p3 in the p1-p2 rest frame.
    //Remember to multiply by 8 pi^2 when integrating to account for d\phi_1 dcos\theta_3 d\phi_3!
    double d_decay_width_2(std::function<double(double,double,double,double,double,double)> amp, double m12, double cos_t, double m0, double m1, double m2, double m3);

    double d_decay_width_3(std::function<double(double,double)> amp, double m12, double cos_t, double m0, double m1, double m2, double m3);

    double integrate_decay_width(std::function<double(double,double)> amp, double m0, double m1, double m2, double m3);

}

#endif

