#include "Distribution.h"
#include "Particle.h"

//Adds an offset to the end position and time of particles. Ideally this would add to the starting position with a propagation Distribution in between, but I don't want to set that up, yet.
//Call this with position_offset
class Position_Offset : public Distribution{
	public:
		Position_Offset(double x, double y, double z, double t) : Distribution(), X(x), Y(y), Z(z), T(t){};
		void sample_particle(Particle &part);
	private:
	   double X,Y,Z,T;	

};
/*
class Position_Offset_Flat : public Distribution{
	public:
		Position_Offset_Flat(double xmin, double xmax, double ymin, double ymax, double zmin, double zmax, double tmin, double tmax): Distribution(), Xmin(xmin), Xmax(xmax), Ymin(xmin), Ymax(xmax), Zmin(xmin), Zmax(xmax), Tmin(xmin), Tmax(xmax);
		void sample_particle(Particle &part);
	private:
		double Xmin, Xmax, Ymin, Ymax, Zmin, Zmax, Tmin, Tmax;
}*/
