#include "decay.h"
#include "constants.h"
#include <iostream>
#include <fstream>
#include "Kinematics.h"


//std::ofstream datalog_decay("decaylog.dat",std::ios::out);

using std::cout; using std::endl;

// two body decay of parent ->  daughter + X (e.g. Pi0 -> gamma + V)
// output momentum of daughter particle in the lab frame
void DecayDP(Particle &daughter, Particle &parent){
	double pdx, pdy, pdz, Ed, pd;
	double thetad, phid;		
	double lam;	
	thetad = acos(-1+2*Random::Flat(0,1));
	phid = Random::Flat(0,2.0*pi*Random::Flat(0,1));
	double Md, Mp;
	Md = daughter.m;
	Mp = parent.m;
	lam = lambda(1,Md*Md/Mp/Mp,0);
	Ed = Mp/2.0*(1.0+Md*Md/Mp/Mp);
	pd = Mp/2.0*sqrt(lam);
	pdx = pd*sin(thetad)*cos(phid);
	pdy = pd*sin(thetad)*sin(phid);
	pdz = pd*cos(thetad);
	daughter.FourMomentum(pdx,pdy,pdz,Ed);
	// boost darkphoton to lab frame
	daughter.Lorentz(parent);
	Link_Particles(parent,daughter);	
}
// cascade decay for parent -> mediator + X, mediator > daughter1 + daughter 2
//This assumes no meaningful propagation before the decay occurs
void DecayDM(Particle &daughter1, Particle &daughter2, Particle &mediator, Particle &parent){
	double pmx, pmy, pmz, Em, pm;
	double ym, zm, thetam, phim;	
	double pdx, pdy, pdz, Ed, pd;
	double yd, zd, thetad, phid;	
	double lam;	
	//double polarX1, polarX2;
	//Particle darkphoton(MV);
	ym = Random::Flat(0,1);
	zm = Random::Flat(0,1);
    thetam = acos(-1+2*ym);
	phim = 2.0*pi*zm;
	double Md, Mm, Mp;
	Md = daughter1.m;
	Mm = mediator.m;
	Mp = parent.m;
	lam = lambda(1,Mm*Mm/Mp/Mp,0);
	Em = Mp/2.0*(1.0+Mm*Mm/Mp/Mp);
	pm = Mp/2.0*sqrt(lam);
	pmx = pm*sin(thetam)*cos(phim);
	pmy = pm*sin(thetam)*sin(phim);
	pmz = pm*cos(thetam);
	mediator.FourMomentum(pmx,pmy,pmz,Em);
    // boost darkphoton to lab frame
	mediator.Lorentz(parent);	
    // generate dark matter momentum in dark photon rest frame
    yd = Random::Flat(0,1);
	zd = Random::Flat(0,1);
	thetad = acos(-1+2*yd);
	phid = 2.0*pi*zd;
	lam = lambda(1,Md*Md/Mm/Mm,Md*Md/Mm/Mm);
	Ed = Mm/2.0;
	pd = Mm/2.0*sqrt(lam);
	pdx = pd*sin(thetad)*cos(phid);
	pdy = pd*sin(thetad)*sin(phid);
	pdz = pd*cos(thetad);
	daughter1.FourMomentum(pdx,pdy,pdz,Ed);
	daughter2.FourMomentum(-pdx,-pdy,-pdz,Ed);
	// boost dark matter particles to lab frame
	daughter1.Lorentz(mediator);	
	daughter2.Lorentz(mediator);
	//cout << parent.name << " " << parent.end_coords[0] << " " << parent.end_coords[1] << " " << parent.end_coords[2] << endl;
	Link_Particles(parent,mediator);
	Link_Particles_Immediate(mediator,daughter1);
	Link_Particles_Immediate(mediator,daughter2);
	//cout << daughter1.name << " " << daughter1.end_coords[0] << " " << daughter1.end_coords[1] << " " << daughter1.end_coords[2] << endl;
	//cout << daughter2.name << " " << daughter2.end_coords[0] << " " << daughter2.end_coords[1] << " " << daughter2.end_coords[2] << endl;
}

//isotropic decay for parent -> daughter1 + daughter2
void TwoBodyDecay(Particle &parent, Particle &daughter1, Particle &daughter2){
    double thetad = acos(Random::Flat(-1,1));
    double phid = Random::Flat(0,1)*2*pi;
    double mom = TriangleFunc(parent.m, daughter1.m, daughter2.m);
    daughter1.ThreeMomentum(mom*sin(thetad)*cos(phid),mom*sin(thetad)*sin(phid),mom*cos(thetad));
    daughter2.ThreeMomentum(-mom*sin(thetad)*cos(phid),-mom*sin(thetad)*sin(phid),-mom*cos(thetad));
    daughter1.Lorentz(parent);
    daughter2.Lorentz(parent);
	Link_Particles(parent,daughter1);
	Link_Particles(parent,daughter2);
}

//Non-isotropic decay for parent -> daughter1 + daughter2.
void TwoBodyDecay(Particle &parent, Particle &daughter1, Particle &daughter2, double thetad){
    double phid = Random::Flat(0,1)*2*pi;
    double mom = TriangleFunc(parent.m, daughter1.m, daughter2.m);
    daughter1.ThreeMomentum(mom*sin(thetad)*cos(phid),mom*sin(thetad)*sin(phid),mom*cos(thetad));
    daughter2.ThreeMomentum(-mom*sin(thetad)*cos(phid),-mom*sin(thetad)*sin(phid),-mom*cos(thetad));
    daughter1.Lorentz(parent);
    daughter2.Lorentz(parent);
	Link_Particles(parent,daughter1);
	Link_Particles(parent,daughter2);
}

void TwoBodyDecay(Particle &parent, Particle &daughter1, Particle &daughter2, double thetad, double phid){
    double mom = TriangleFunc(parent.m, daughter1.m, daughter2.m);
    daughter1.ThreeMomentum(mom*sin(thetad)*cos(phid),mom*sin(thetad)*sin(phid),mom*cos(thetad));
    daughter2.ThreeMomentum(-mom*sin(thetad)*cos(phid),-mom*sin(thetad)*sin(phid),-mom*cos(thetad));
    daughter1.Lorentz(parent);
    daughter2.Lorentz(parent);
	Link_Particles(parent,daughter1);
	Link_Particles(parent,daughter2);
}

/*I should check that this is handled properly. Not quite sure that it is, though
it probably won't change very much if the angle is handled differently.
*/
void DecayDM_Off_Shell(Particle &daughter1, Particle &daughter2, Particle &mediator, Particle &parent, double theta){
    TwoBodyDecay(mediator, daughter1, daughter2, theta);
    mediator.ThreeMomentum(0,0,TriangleFunc(parent.m,mediator.m,0));
    daughter1.Lorentz(mediator);
    daughter2.Lorentz(mediator);
    double mtheta=acos(Random::Flat(-1,1));
    double mphi = Random::Flat(0,1)*2*pi;
    mediator.Rotate_z(mtheta);
    mediator.Rotate_y(mphi);
    daughter1.Rotate_z(mtheta);
    daughter1.Rotate_y(mphi);
    daughter2.Rotate_z(mtheta);
    daughter2.Rotate_y(mphi);
    mediator.Lorentz(parent);
    daughter1.Lorentz(parent);
    daughter2.Lorentz(parent);
	Link_Particles(parent,mediator);
	Link_Particles_Immediate(mediator,daughter1);
	Link_Particles_Immediate(mediator,daughter2);
}

void Meson_Capture_Off_Shell(Particle &daughter1, Particle &daughter2, Particle &mediator, double theta){
    double thetad = mediator.Theta();
    double phid = mediator.Phi();
    mediator.Rotate_y(-phid);
    mediator.Rotate_z(-thetad);
    TwoBodyDecay(mediator, daughter1, daughter2, theta); 
    mediator.Rotate_z(thetad);
    mediator.Rotate_y(phid);
    daughter1.Rotate_z(thetad);
    daughter1.Rotate_y(phid);
    daughter2.Rotate_z(thetad);
    daughter2.Rotate_y(phid);
	Link_Particles(mediator,daughter1);
	Link_Particles(mediator,daughter2);
}
