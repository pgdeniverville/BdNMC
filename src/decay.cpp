#include "decay.h"
#include "constants.h"
#include "Kinematics.h"
#include "Integrator.h"

#include <iostream>
#include <fstream>

//std::ofstream datalog_decay("decaylog.dat",std::ios::out);

using std::cout; using std::endl;

using namespace std::placeholders;

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

namespace Three_Body_Decay_Space{

    double E2star(double m12s, double m1, double m2){
        return (m12s-m1*m1+m2*m2)/(2*sqrt(m12s));
    }

    //m0 is the mass of the parent particle.
    double E3star(double m12s, double m0, double m3){
        return (m0*m0-m12s-m3*m3)/(2*sqrt(m12s));
    }

    double p2star(double m12s, double m1, double m2){
        return sqrt(pow(E2star(m12s,m1,m2),2)-m2*m2);
    }

    double p3star(double m12s, double m0, double m3){
        return sqrt(pow(E3star(m12s,m0,m3),2)-m3*m3);
    }

    double m23min(double m12s, double m0, double m1, double m2, double m3){
        return pow(E2star(m12s,m1,m2)+E3star(m12s,m0,m3),2)-pow(sqrt(E3star(m12s,m0,m3)-m3*m3)+sqrt(E2star(m12s,m1,m2)-m2*m2),2);
    }

    double m23max(double m12s, double m0, double m1, double m2, double m3){
        return pow(E2star(m12s,m1,m2)+E3star(m12s,m0,m3),2)-pow(sqrt(E3star(m12s,m0,m3)-m3*m3)-sqrt(E2star(m12s,m1,m2)-m2*m2),2);
    }
    //!!!!!Note that this is using m12, not m12s=m12^2!!!!
    double p1star(double m12, double m1, double m2){
        return sqrt((m12*m12-pow(m1+m2,2))*(m12*m12-pow(m1-m2,2)))/(2*m12);
    }

    double p3rest(double m12, double m0, double m3){
        return sqrt((m0*m0-pow(m12+m3,2))*(m0*m0-pow(m12-m3,2)))/(2*m0);
    }

    //p2.p3 dot product evaluated in p1+p2 rest frame!
    double p2p3(double m12s, double CosTheta, double m0, double m1, double m2, double m3){
        return E2star(m12s,m1,m2)*E3star(m12s,m0,m3)-p2star(m12s,m1,m2)*p3star(m12s,m0,m3)*CosTheta;
    }

    double Cos_Theta_to_m23s(double m12s, double CosTheta, double m0, double m1, double m2, double m3){
        return m2*m2+m3*m3+2*p2p3(m12s,CosTheta,m0,m1,m2,m3);
    }

    //Eq 47.22 in the PDG.
    //Function expects amp(m12s, s, m0, m1, m2, m3)
    double d_decay_width(std::function<double(double,double,double,double,double,double)> amp, double m12s, double m23s, double m0, double m1, double m2, double m3){
        return amp(m12s,m23s,m0,m1,m2,m3)/(32*pow(2*pi*m0,3));
    }

    //This includes the angular dependence
    //EQ 47.20 in the PDG
    //amp(m12s,s,m0, m1,m2,m3). 
    //cos_t is the angle between p1 and p3 in the p1-p2 rest frame.
    //Remember to multiply by 8 pi^2 when integrating to account for d\phi_1 dcos\theta_3 d\phi_3!
    double d_decay_width_2(std::function<double(double,double,double,double,double,double)> amp, double m12, double cos_t, double m0, double m1, double m2, double m3){
        if(m12<=m1+m2){
            return 0;
        }
        //cout << endl << "d_decay_width_2 test " << amp(m12*m12,Cos_Theta_to_m23s(m12*m12, cos_t, m0, m1, m2, m3),m0,m1,m2,m3) << " " << m12*m12 << " " << Cos_Theta_to_m23s(m12*m12, cos_t, m0, m1, m2, m3) << " " << p2p3(m12*m12,cos_t,m0,m1,m2,m3) << " " << p1star(m12,m1,m2) << " " << p3rest(m12,m0,m3) << " " << 1/(pow(2*pi,5)*16*m0*m0) << endl;
        return amp(m12*m12,Cos_Theta_to_m23s(m12*m12, cos_t, m0, m1, m2, m3),m0,m1,m2,m3)/(pow(2*pi,5)*16*m0*m0)*p1star(m12,m1,m2)*p3rest(m12,m0,m3);
    }

    double d_decay_width_3(std::function<double(double,double)> amp, double m12, double cos_t, double m0, double m1, double m2, double m3){
        if(m12<=m1+m2){
            return 0;
        }
        //cout << endl << "d_decay_width_2 test " << amp(m12*m12,Cos_Theta_to_m23s(m12*m12, cos_t, m0, m1, m2, m3),m0,m1,m2,m3) << " " << m12*m12 << " " << Cos_Theta_to_m23s(m12*m12, cos_t, m0, m1, m2, m3) << " " << p2p3(m12*m12,cos_t,m0,m1,m2,m3) << " " << p1star(m12,m1,m2) << " " << p3rest(m12,m0,m3) << " " << 1/(pow(2*pi,5)*16*m0*m0) << endl;
        return amp(m12*m12,Cos_Theta_to_m23s(m12*m12, cos_t, m0, m1, m2, m3))/(pow(2*pi,5)*16*m0*m0)*p1star(m12,m1,m2)*p3rest(m12,m0,m3);
    }

    //Need to test this, may just want to do a random sample instead.
    double integrate_decay_width(std::function<double(double,double)> amp, double m0, double m1, double m2, double m3, double r_accuracy_goal){
        if(m0<m1+m2+m3){
            return 0;
        }

        std::function<double(double, double)> d2amp = bind(d_decay_width_3,amp,_1,_2,m0,m1,m2,m3);

        //return 8*pi*pi*RandomIntegrate2(d2amp,m1+m2,m0-m3,-1,1,1e6);
        //return 8*pi*pi*SimpsonCubature(d2amp,m1+m2,m0-m3,100,-1,1,100);
        //return 8*pi*pi*SimpsonCubature_adapt(d2amp,m1+m2,m0-m3,100,-1,1,100);
        return 8*pi*pi*RandomIntegrate2_adapt(d2amp,m1+m2,m0-m3,-1,1,1e5,r_accuracy_goal);
    }

    void Three_Body_Decay(Particle &parent, Particle &daughter1, Particle &daughter2, Particle &daughter3, double &d_width_max, std::function<double(double, double, double, double, double, double)> &amp){
        
        double d_width = 0;
        double m12, cos_t;
        double m0=parent.m;
        double m1=daughter1.m;
        double m2=daughter2.m;
        double m3=daughter3.m;

        do{
            //Very specifically using m12
            m12=Random::Flat(m1+m2,m0-m3);
            cos_t=Random::Flat(-1,1);
            d_width=d_decay_width_2(amp,m12,cos_t,m0,m1,m2,m3);
        } while(d_width<Random::Flat()*d_width_max);

        if(d_width>d_width_max){
            d_width_max=d_width;
        }

        double phi_1 = Random::Flat(0,2*pi);
        double phi_3 = Random::Flat(0,2*pi);
        double cos_t_3 = Random::Flat(-1,1);

        double t_1 = acos(cos_t);
        //double t_2 = t_1+pi;

        double m12s = m12*m12;

        //transform to rest frame, need to check that this is the right sign.
        double m12_mom3=p3star(m12s,m0,m3);
        double bet=m12_mom3/sqrt(pow(m12_mom3,2)+m0*m0);
//      double gam=beta_gamma(bet);
        double m12_mom = TriangleFunc(m12,m1,m2);


  //      cout << "Testing before everything...\n";

        daughter1.ThreeMomentumPolar(m12_mom,t_1,phi_1);
        daughter2.ThreeMomentumPolar(m12_mom,t_1+pi,phi_1);
        daughter3.ThreeMomentumPolar(m12_mom3,0,0);
/*
        daughter1.report(cout);
        daughter2.report(cout);
        daughter3.report(cout);
*/
        //Transform from theta^* to theta_rest
        //t_1 = atan(sin(t_1)/cos(t_1)/gam);
        //t_2 = atan(sin(t_2)/cos(t_2)/gam);
        double t_3 = acos(cos_t_3);

        //Need to check that these rotations are in the right direction.
        //daughter1.ThreeMomentumPolar(sqrt(pow(E1_rest,2)-m1*m1), t_1, phi_1);
        //daughter2.ThreeMomentumPolar(sqrt(pow(E2_rest,2)-m2*m2), t_2, phi_1);

        //daughter1.report(cout);
        //daughter2.report(cout);
        daughter1.Lorentz(bet,0,0,-bet);
        daughter2.Lorentz(bet,0,0,-bet);
        daughter3.Lorentz(bet,0,0,-bet);
/*
        cout << "Before rotation...\n";

        daughter1.report(cout);
        daughter2.report(cout);
        daughter3.report(cout);        
*/

        daughter1.Rotate_y(t_3);
        daughter1.Rotate_z(phi_3);
        daughter2.Rotate_y(t_3);
        daughter2.Rotate_z(phi_3);
      // daughter3.Rotate_y(t_3);
      //  daughter3.Rotate_z(phi_3);
/*
        cout << "4-momentum report in meson rest frame\n";
        parent.report(cout);
        daughter1.report(cout);
        daughter2.report(cout);
        daughter3.report(cout);
*/

        //Creating Daughter3!
        daughter3.ThreeMomentumPolar(p3rest(m12, m0, m3),t_3,phi_3);
/*
        cout << "Compare alternate c3...\n";

        daughter3.report(cout);
*/
        //boost to lab frame!
        daughter1.Lorentz(parent);
        daughter2.Lorentz(parent);
        daughter3.Lorentz(parent);
/*
        cout << "4-momentum report in lab frame\n";
        parent.report(cout);
        daughter1.report(cout);
        daughter2.report(cout);
        daughter3.report(cout);

        cout << endl << endl;
*/
        //link the particles so they have the right initial coordinates.
        Link_Particles(parent, daughter1);
        Link_Particles(parent, daughter2);
        Link_Particles(parent, daughter3);

        //At this point, each of the daughter particles should be given a chance to decay by its own DMGenerator. Handled in Three_Body_Decay_Gen.cpp!
    }
}