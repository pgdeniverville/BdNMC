#include "Particle_Generator.h"
#include <list>

void Particle_Generator::Generate_Particle(Particle &part){
	part.Set_Mass(mass_particle);
	dist->Sample_Particle(part);
}
Particle_Generator::Particle_Generator(double m, std::shared_ptr<Distribution> DIST){
	mass_particle=m; dist=DIST;
}

