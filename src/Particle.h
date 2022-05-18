// class for lorentz transformation
#ifndef GUARD_Particle_h
#define GUARD_Particle_h

#include "Random.h"

#include <cmath>
#include <memory>
#include <string>
#include <iostream>
#include <functional>

class Particle{
	
public:
    Particle(double);
    Particle(){Particle(0);}
    Particle(const Particle &);

    Particle& operator=(const Particle& part);
    

    double m=0;
    double px=0, py=0, pz=0, E=0;//Eventually, these variables are going to be made private.
    double width=0;
    std::string name;
    double   Theta();
    double   Phi();
    double Width(){return width;} 
    void     Set_Mass(double);
    void     FourMomentum(double, double, double, double);
    void     ThreeMomentum(double PX, double PY, double PZ);
    void	 ThreeMomentumPolar(double mom, double theta, double phi);
	void     Lorentz(Particle&);
    void     Lorentz(double beta, double betax, double betay, double betaz);
    void     Rotate_x(double);
    void     Rotate_y(double);
    void     Rotate_z(double);
    void     Rotate(Particle&);
    void     report(std::ostream&) const;
    void     report() const {report(std::cout);}
    //void    Generate_4Vector(double s);
	void 	Set_Origin(double x, double y, double z);
//    void    Set_Position(double x, double y, double z);
	void 	Set_Creation_Time(double t);
    //If END_SET=true, this also updates dec_time.
    void    Set_Time(double t);
    void    Increment_Time(double t){END_SET=true; Set_Time(t+origin_coords[3]);} 
    double Momentum() const;
    double Speed();
    double Kinetic_Energy() const {return E-m;}
    void Generate_Position();//Generates an end_coords position between crossing[0] and crossing[1].
    void Generate_Position(double);//Generates an end_coords position
    //I should design custom 4 vector objects. Also 3 vector objects which 4-vectors hold.
    double origin_coords[4];//space-time location where the particle was created
	double end_coords[4];//Stores the space-time location where the particle scatters or is otherwise destroyed.
    double crossing[2];//This keeps track of when a particle crosses through a detector. This will eventually replace det_cross in detector objects!
    double dec_time=0;//tracks when a particle is destroyed by scattering or decay if END_SET=true. Multiply by momentum() to retrieve distance.
    bool END_SET=false;//Indicates that this particle decays or scatters at a distance of Momentum()*dec_time.
    bool EVENT_SET=false;//I am going to have to be careful with this variable.
    
    //Going to revamp output soon. my_id will be its position in the particle list. Every particle will keep track of its parent. Option for both parents if produced in scattering, but parent_2 mostly unused.
    int my_id=-1;
    int parent_1_id=-1;
    int parent_2_id=-1;
};

void Link_Particles(Particle &, Particle &);
void Link_Particles_Immediate(Particle &, Particle &);

#endif
