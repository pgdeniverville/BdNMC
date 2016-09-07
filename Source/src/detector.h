#ifndef GUARD_detector_h
#define GUARD_detector_h

// class detector
#include "Particle.h"
#include <string>
#include <vector>
class Material {
    public:
        Material(double nd, double np, double nn, double ne, double m, std::string name) {nDensity=nd; Proton_Number=np; Neutron_Number=nn; Electron_Number=ne; matname=name; mass=m;}
       /* double get_Density() {return nDensity;}
        double PN() {return Proton_Number;}
        double NN() {return Neutron_Number;}
        double EN() {return Electron_Number;}
        double matmass() {return mass;}
        std::string getname() {return name;}*/
   // private:
        double nDensity, Proton_Number, Neutron_Number, Electron_Number, mass;
        std::string matname;
};

class detector {
    //double xdet, ydet, zdet, Rdet;
public:
    virtual double Ldet (const Particle &) = 0;
//	virtual void intersect (int &, int &, Particle) = 0;
    detector(){ p_num_tot=0; n_num_tot=0; e_num_tot=0;}
    virtual ~detector(){};
    void add_material(double nd, double np, double nn, double ne, double m, std::string matname);
    double get_nDensity(int index) {return matvec[index].nDensity;}
    double PN(int index) {return matvec[index].Proton_Number;}//These should be an int most of the time.
    double NN(int index) {return matvec[index].Neutron_Number;}
    double EN(int index) {return matvec[index].Electron_Number;}
    double M(int index) {return matvec[index].mass;}
    std::string matname(int index) {return matvec[index].matname;}
    int mat_num() {return matvec.size();}
    
    double PNtot(){ return p_num_tot;}
    double NNtot(){ return n_num_tot;}
    double ENtot(){ return e_num_tot;}
    
    double cross_point[2];//holds the last entrance and exit points of a particle.
protected:
    double r[3];
    double b[3];
	double o[3];
    std::vector<Material> matvec;
    double p_num_tot, n_num_tot, e_num_tot;
};

class detector_sphere: public detector{
public:
    detector_sphere(double x, double y, double z, double R);
    ~detector_sphere(){}
    double Ldet (const Particle &);
//  void intersect (const int &, const int &, const Particle &);
private:
    double Rdet;
};

class detector_cylinder: public detector{
    public:
        detector_cylinder(double xdet, double ydet, double zdet, double detlength, double radius, 
                double detTheta, double detPhi);
        ~detector_cylinder(){}
        double Ldet (const Particle &);
        double Ldeto (const Particle &, const double offset[3]);
    private:
        double Rdet, Ldetector;
        double l[3];
};

#endif