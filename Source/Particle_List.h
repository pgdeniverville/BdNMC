#ifndef PARTICLE_LIST_H
#define PARTICLE_LIST_H

#include <string>
#include <iostream>
#include <fstream>
#include <list>
#include <memory>
#include "Distribution.h"

struct _3momentum{
	double px,py,pz,E;
	_3momentum(double x, double y, double z, double e){px=x; py=y; pz=z; E=e;};
};

class Particle_List: public Distribution{

public:
	Particle_List(const std::string&);
	~Particle_List(){instream.close();}
	//Get_Particle returns the momentum and direction of a particle from a list.
	void sample_momentum(double& pmom, double& theta, double& phi);
	void Get_Particle(double& pmom, double& theta, double& phi);
private:
	//determines how many particles should be read in from file at once.
	std::ifstream instream;
	//int batch_size;
	//std::list<std::shared_ptr<_3momentum> > partlist;
	//std::list<std::shared_ptr<_3momentum> >::iterator iter;
	std::list<_3momentum> partlist;
	std::list<_3momentum>::iterator iter;
	int parse_line(std::string &line);
	int load_particle_batch();
	int part_count;
	//int load_particle_batch();
};



#endif
