#ifndef DISTRIBUTION_H
#define DISTRIBUTION_H

#include "Particle.h"

//Distribution can alter the origin_coords, momentum or both.
class Distribution {
	public:
		Distribution(){};
		//virtual void sample_momentum(double& p, double& theta, double& phi) = 0;
		virtual void sample_particle(Particle &part)=0;
	protected:
};

#endif
