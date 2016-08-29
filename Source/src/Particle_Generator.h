#ifndef GUARD_PARTICLE_GENERATOR_H
#define GUARD_PARTICLE_GENERATOR_H

#include "Particle.h"
#include "Distribution.h"
#include <memory>

class Particle_Generator{
public:
	Particle_Generator(double m, std::shared_ptr<Distribution> DIST);
	void Generate_Particle(Particle &);
private:
	std::shared_ptr<Distribution> dist;
	double mass_particle;
};

#endif
