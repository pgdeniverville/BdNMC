#include "Distribution.h"

void BeamDistribution::sample_particle(Particle& part){
	part.Set_Mass(part.m);
	part.ThreeMomentum(0,0,sqrt(pow(beam_energy,2)-pow(beam_mass,2)));
}