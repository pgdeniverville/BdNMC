#ifndef DISTRIBUTION_H
#define DISTRIBUTION_H

#include "Particle.h"
#include <list>
#include <memory>

//Distribution can alter the origin_coords, momentum or both.
class Distribution {
	public:
		Distribution(){part_mass=0;}
		void Sample_Particle(Particle &part)
		{
			part.Set_Mass(part_mass);
			sample_particle(part);
			for(std::list<std::shared_ptr<Distribution> >::iterator i = dist_list.begin(); i!=dist_list.end(); i++){
				//part.report(std::cout);
				(*i)->sample_particle(part);
			}
		}
		void set_mass(double m){part_mass=m;}
		void Add_Dist(std::shared_ptr<Distribution> distptr){dist_list.push_back(distptr);}//This should be unique, need to figure out proper syntax.
	protected:
		virtual void sample_particle(Particle &part) = 0;
		std::list<std::shared_ptr<Distribution> > dist_list;
	private:
		double part_mass;
};

class DoNothingDist : public Distribution{
    public:
        void sample_particle(Particle &part){return;} 
};

#endif
