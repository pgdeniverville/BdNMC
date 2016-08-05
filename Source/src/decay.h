#ifndef GUARD_decay_h
#define GUARD_decay_h

// class decay

#include "Particle.h"
#include "Random.h"

void DecayDP (Particle &, Particle &);
void DecayDM (Particle &, Particle &, Particle &, Particle &);
void DecayDM_Off_Shell(Particle &daughter1, Particle &daughter2, Particle &mediato, Particle &parent, double theta);
void TwoBodyDecay (Particle &parent, Particle &daughter1, Particle &daughter2);
void TwoBodyDecay (Particle &parent, Particle &daughter1, Particle &daughter2, double thetad);
void TwoBodyDecay (Particle &parent, Particle &daughter1, Particle &daughter2, double thetad, double phid);
void Meson_Capture_Off_Shell(Particle &daughter1, Particle &daughter2, Particle &mediator, double theta);
#endif

