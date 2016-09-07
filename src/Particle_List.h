#ifndef PARTICLE_LIST_H
#define PARTICLE_LIST_H

#include <string>
#include <iostream>
#include <fstream>
#include <list>
#include <memory>
#include "Distribution.h"

struct _part{
	double px,py,pz,E,x,y,z,t;
	_part(double Px, double Py, double Pz, double e, double X, double Y, double Z, double T){px=Px; py=Py; pz=Pz; E=e;x=X;y=Y;z=Z;t=T;}
};

class Particle_List: public Distribution{

public:
	Particle_List(const std::string&, bool set_pos = false);
	~Particle_List(){instream.close();}
	//Get_Particle returns the momentum and direction of a particle from a list.
	void sample_particle(Particle &part);
	void Get_Particle(double& pmom, double& theta, double& phi);
private:
	//determines how many particles should be read in from file at once.
	std::ifstream instream;
	//int batch_size;
	std::list<_part> partlist;
	std::list<_part>::iterator iter;
	int parse_line(std::string &line);
	int load_particle_batch();
	int part_count;
	bool POS;
	//int load_particle_batch();
};



#endif
