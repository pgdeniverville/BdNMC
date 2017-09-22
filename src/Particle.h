// class for lorentz transformation
#ifndef GUARD_Particle_h
#define GUARD_Particle_h

#include <cmath>
#include <memory>
#include <string>
#include <iostream>
#include <functional>

#include "Random.h"

class Particle{
	
public:
    Particle(double);
    Particle(){Particle(0);}
    Particle(const Particle &);

    Particle& operator=(const Particle& part);
    

    double   m, px, py, pz, E;//Eventually, these variables are going to be made private.
    std::string name;
    double   Theta();
    double   Phi(); 
    void     Set_Mass(double);
    void     FourMomentum(double, double, double, double);
    void     ThreeMomentum(double PX, double PY, double PZ);
    void	 ThreeMomentumPolar(double mom, double theta, double phi);
	void     Lorentz(Particle&);
    void     Rotate_x(double);
    void     Rotate_y(double);
    void     Rotate_z(double);
    void     report(std::ostream&) const;
    //void    Generate_4Vector(double s);
	void 	Set_Origin(double x, double y, double z);
//    void    Set_Position(double x, double y, double z);
	void 	Set_Creation_Time(double t);
    void    Set_Time(double t); 
    double Momentum();
    double Speed();
    double Kinetic_Energy(){return E-m;}
    void Generate_Position();//Generates an end_coords position between crossing[0] and crossing[1].
    void Generate_Position(double);//Generates an end_coords position
    //I should design custom 4 vector objects. Also 3 vector objects which 4-vectors hold.
    double origin_coords[4];//space-time location where the particle was created
	double end_coords[4];//Stores the space-time location where the particle scatters or is otherwise destroyed.
    double crossing[2];//This keeps track of when a particle crosses through a detector. This will eventually replace det_cross in detector objects!
    bool EVENT_SET;//I am going to have to be careful with this variable.
};

void Link_Particles(Particle &, Particle &);
void Link_Particles_Immediate(Particle &, Particle &);

#endif
