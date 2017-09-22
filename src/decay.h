#ifndef GUARD_decay_h
#define GUARD_decay_h

#include "Particle.h"
#include "Random.h"

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
#endif

