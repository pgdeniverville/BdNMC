#include "Particle.h"
#include "Kinematics.h"
#include <iostream>
#include "constants.h"
#include <iomanip>

using std::cout; using std::endl;

//This moves a daughter particle's origin to the parent particle's end.
void Link_Particles(Particle &parent, Particle &child){
	child.Set_Origin(parent.end_coords[0],parent.end_coords[1],parent.end_coords[2]);
	child.Set_Creation_Time(parent.end_coords[3]);
    child.parent_1_id=parent.my_id;
    child.EVENT_SET = parent.EVENT_SET;
}

void Link_Particles_Immediate(Particle &parent, Particle &child){
	//parent.Set_Position(parent.origin_coords[0],parent.origin_coords[1],parent.origin_coords[2]);
	parent.Set_Time(parent.origin_coords[3]);
	child.Set_Origin(parent.end_coords[0],parent.end_coords[1],parent.end_coords[2]);
	child.Set_Creation_Time(parent.end_coords[3]);
    child.parent_1_id=parent.my_id;
    child.EVENT_SET = parent.EVENT_SET;
}


//using std::cout; using std::endl;
//Should update to use Set_Origin and Set_Time
Particle::Particle(double mass){
    Set_Mass(mass);
    ThreeMomentum(0,0,0);
    for(int i=0; i<4; i++){
        origin_coords[i]=0.0;
		end_coords[i]=0.0;
	}
    crossing[0]=0; crossing[1]=0;
    EVENT_SET = false;
    END_SET = false;
}

//Should one of these just call the other one?
Particle::Particle(const Particle &part){
	Set_Mass(part.m);
	ThreeMomentum(part.px,part.py,part.pz);
	for(int i=0; i<4; i++){
        origin_coords[i]=part.origin_coords[i];
		end_coords[i]=part.end_coords[i];
	}
    crossing[0]=part.crossing[0];
    crossing[1]=part.crossing[1];
	EVENT_SET = part.EVENT_SET;
    END_SET = part.END_SET;
    dec_time = part.dec_time;
	name = part.name;
    width = part.width;
}

Particle& Particle::operator=(const Particle& part){
	Set_Mass(part.m);
	ThreeMomentum(part.px,part.py,part.pz);
	for(int i=0; i<4; i++){
        origin_coords[i]=part.origin_coords[i];
		end_coords[i]=part.end_coords[i];
	}
	EVENT_SET = part.EVENT_SET;
    END_SET = part.END_SET;
    dec_time = part.dec_time;
	width = part.width;
    name = part.name;
    crossing[0]=part.crossing[0];
    crossing[1]=part.crossing[1];
	return *this;
}

double Particle::Theta(){
    return calculatetheta(px, py, pz, E); 
}  

double Particle::Phi(){
    return phi(px, py, pz, E);
}

void Particle::Rotate(Particle &part){
    Rotate_y(part.Theta());
    Rotate_z(part.Phi());
}

void Particle::FourMomentum(double PX, double PY, double PZ, double P0){

	px = PX;
	py = PY;
	pz = PZ;
	E = P0;
}

void Particle::ThreeMomentum(double PX, double PY, double PZ){
    
    px = PX;
    py = PY;
    pz = PZ;
    E = sqrt(m*m+px*px+py*py+pz*pz);

}

void Particle::ThreeMomentumPolar(double mom, double theta, double phi){
	ThreeMomentum(mom*cos(phi)*sin(theta),mom*sin(phi)*sin(theta),mom*cos(theta));
}

void Particle::Set_Mass(double m_new){
    m = m_new;
    E = sqrt(m_new*m_new+px*px+py*py+pz*pz);
}
//Currently, I am not enforcing any joining of the origin coordinates to the end coordinates.
//This should be corrected. Probably by not allowing end point to be set at all, only the time?
//Need to update end point to reflect new origin point.
void Particle::Set_Origin(double x, double y, double z){
	origin_coords[0]=x;
	origin_coords[1]=y;
	origin_coords[2]=z;
}
//This should be handled by Set_Time.
//It's been removed, will probably break a bunch of things, but that's okay.
/*
void Particle::Set_Position(double x, double y, double z){
    end_coords[0]=x;
    end_coords[1]=y;
    end_coords[2]=z;
}
*/
void Particle::Set_Creation_Time(double t){
	origin_coords[3]=t;
}

void Particle::Set_Time(double t){
    if(END_SET){
        dec_time = t*Speed()*speed_of_light/Momentum();
    }

//    cout << "Testing Set_Time" << endl;
//    cout << "Time to be set=" << t << endl;
//    cout << "origin coords = " << origin_coords[0] << " " << origin_coords[1] << " " << origin_coords[2] << " "  <<  origin_coords[3] << endl;

    end_coords[0]=(t-origin_coords[3])*Speed()*speed_of_light*px/Momentum()+origin_coords[0];
    end_coords[1]=(t-origin_coords[3])*Speed()*speed_of_light*py/Momentum()+origin_coords[1];
    end_coords[2]=(t-origin_coords[3])*Speed()*speed_of_light*pz/Momentum()+origin_coords[2];
    end_coords[3]=t;
}

//Optimize the lorentz transformations at a later date, I can squeeze more speed out of these.
void Particle::Lorentz(double beta, double betax, double betay, double betaz){
    
    if(beta==0.0)
        return;
    double gamma = 1/sqrt(1-beta*beta);

    double nx = betax/beta;
    double ny = betay/beta;
    double nz = betaz/beta;

    double Lam11 = gamma;
    double Lam12 = gamma*betax;
    double Lam13 = gamma*betay;
    double Lam14 = gamma*betaz;
    double Lam21 = Lam12;
    double Lam22 = 1+(gamma-1)*nx*nx;
    double Lam23 = (gamma-1)*nx*ny;
    double Lam24 = (gamma-1)*nx*nz;
    double Lam31 = Lam13;
    double Lam32 = Lam23;
    double Lam33 = 1+(gamma-1)*ny*ny;
    double Lam34 = (gamma-1)*ny*nz;
    double Lam41 = Lam14;
    double Lam42 = Lam24;
    double Lam43 = Lam34;
    double Lam44 = 1+(gamma-1)*nz*nz;

    double E2   = Lam11*E+Lam12*px+Lam13*py+Lam14*pz;
    double p2x  = Lam21*E+Lam22*px+Lam23*py+Lam24*pz;  
    double p2y  = Lam31*E+Lam32*px+Lam33*py+Lam34*pz;  
    double p2z  = Lam41*E+Lam42*px+Lam43*py+Lam44*pz;
        
    FourMomentum(p2x, p2y, p2z, E2);
}
//Need to test that this works, but should be fine!

void Particle::Lorentz(Particle& parent){
    double gamma = parent.E/parent.m;
    double beta = sqrt(1.0-1.0/gamma/gamma);
    if(beta==0)
        return;
    Lorentz(beta, parent.px/parent.E, parent.py/parent.E, parent.pz/parent.E);
}

void Particle::Rotate_x(double psi){
    double zr = pz*cos(psi)+py*sin(psi);
    double yr = py*cos(psi)-pz*sin(psi);
    pz = zr; py = yr;
}

void Particle::Rotate_y(double theta){
    double zr = pz*cos(theta)-px*sin(theta);
    double xr = px*cos(theta)+pz*sin(theta);
    px=xr; pz=zr;
}

void Particle::Rotate_z(double phi){
    double yr = py*cos(phi)+px*sin(phi);
    double xr = px*cos(phi)-py*sin(phi);
    px=xr; py=yr;
}

double Particle::Momentum() const  {
    return sqrt(px*px+py*py+pz*pz);
}
//Speed is in units of c
double Particle::Speed(){
    if(Momentum()==0.0)
        return 0.0;
    return Momentum()/E;
}

void Particle::Generate_Position(){
    Generate_Position(Random::Flat(crossing[0],crossing[1]));
}

void Particle::Generate_Position(double rngpoint){
    if(Momentum()==0.0){
        Set_Time(origin_coords[3]);
    }
    else{
        Set_Time(Momentum()*rngpoint/Speed()/speed_of_light+origin_coords[3]);
    }
}

void Particle::report(std::ostream& ostr) const{
    ostr << std::setprecision(10) << name << " " << E << " " << px << " " << py << " " << pz << " " << m << " " << origin_coords[0] << " " << origin_coords[1] << " " << origin_coords[2] << " " << origin_coords[3] << " " << end_coords[0] << " " << end_coords[1] << " " << end_coords[2] << " " << end_coords[3] << " " << Kinetic_Energy() << " " << Momentum()*crossing[0] << " " << Momentum()*crossing[1] << std::endl;
} 
