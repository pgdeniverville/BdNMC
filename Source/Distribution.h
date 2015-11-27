#ifndef DISTRIBUTION_H
#define DISTRIBUTION_H

#include "Particle.h"
#include <list>
#include <memory>

//Distribution can alter the origin_coords, momentum or both.
class Distribution {
	public:
		Distribution(){
			dist_list = std::shared_ptr<std::list<Distribution> > (new std::list<Distribution>);	
		};
		//virtual void sample_momentum(double& p, double& theta, double& phi) = 0;
		void Sample_Particle(Particle &part)
		{
			sample_particle(part);
			for(std::list<Distribution>::iterator i = dist_list->begin(); i!=dist_list->end(); i++)
				i->sample_particle(part);
		}
		void Add_Dist(Distribution &dist){dist_list->push_back(dist);}//This might be dangerous.
	protected:
		void sample_particle(Particle &part){};
		std::shared_ptr<std::list<Distribution> > dist_list;
};

#endif
