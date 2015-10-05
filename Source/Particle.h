// class for lorentz transformation
#ifndef GUARD_Particle_h
#define GUARD_Particle_h

#include <cmath>
#include <memory>
#include <string>
#include <iostream>
#include <functional>
#include "Distribution.h"
class Particle{
	
public:
    Particle(double);
    double   m, px, py, pz, E;//Eventually, these variables are going to be made private.
    std::string name;
    double   Theta();
    double   Phi(); 
    void     Set_Mass(double);
    void     FourMomentum(double, double, double, double);
    void     ThreeMomentum(double PX, double PY, double PZ);
    void     Lorentz(Particle);
    void     Rotate_x(double);
    void     Rotate_y(double);
    void     Rotate_z(double);
    void     report(std::ostream&);
    void    Generate_4Vector(double s);
    void    Set_Position(double x, double y, double z);
    void    Set_Time(double t); 
    double Momentum();
    double Speed();
    void Generate_Position(double);
    double coords[4];//Stores the space-time location of the particle.
    bool EVENT_SET;//I am going to have to be careful with this variable.
};
//Not sure where else to put this.

class Particle_Generator{
public:
	Particle_Generator(double m, std::shared_ptr<Distribution> DIST){mass_particle=m; dist=DIST;}
	void Generate_Particle(Particle &);
private:
	std::shared_ptr<Distribution> dist;
	double mass_particle;
};

#endif
