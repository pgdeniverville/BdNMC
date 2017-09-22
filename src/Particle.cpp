#include "Particle.h"
#include "Kinematics.h"
#include <iostream>
#include "constants.h"

//This moves a daughter particle's origin to the parent particle's end.
void Link_Particles(Particle &parent, Particle &child){
	child.Set_Origin(parent.end_coords[0],parent.end_coords[1],parent.end_coords[2]);
	child.Set_Creation_Time(parent.end_coords[3]);
}

void Link_Particles_Immediate(Particle &parent, Particle &child){
	//parent.Set_Position(parent.origin_coords[0],parent.origin_coords[1],parent.origin_coords[2]);
	parent.Set_Time(parent.origin_coords[3]);
	child.Set_Origin(parent.end_coords[0],parent.end_coords[1],parent.end_coords[2]);
	child.Set_Creation_Time(parent.end_coords[3]);
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
}

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
	name = part.name;
}

Particle& Particle::operator=(const Particle& part){
	Set_Mass(part.m);
	ThreeMomentum(part.px,part.py,part.pz);
	for(int i=0; i<4; i++){
        origin_coords[i]=part.origin_coords[i];
		end_coords[i]=part.end_coords[i];
	}
	EVENT_SET = part.EVENT_SET;
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
    end_coords[0]=(t-origin_coords[3])*Speed()*speed_of_light*px/Momentum()+origin_coords[0];
    end_coords[1]=(t-origin_coords[3])*Speed()*speed_of_light*py/Momentum()+origin_coords[1];
    end_coords[2]=(t-origin_coords[3])*Speed()*speed_of_light*pz/Momentum()+origin_coords[2];
    end_coords[3]=t;
}

void Particle::Lorentz(Particle& parent){
	
	double M1 = parent.m;	
	double E1 = parent.E;
	double p1x = parent.px;
	double p1y = parent.py;
	double p1z = parent.pz;
	double E2 = E;	
	double p2x = px;
	double p2y = py;
	double p2z = pz;
	double E3, p3x, p3y, p3z;		
	double gamma, v;	
	double vx, vy, vz;	
	double Lam11, Lam12, Lam13, Lam14;	
	double Lam21, Lam22, Lam23, Lam24;	
	double Lam31, Lam32, Lam33, Lam34;
	double Lam41, Lam42, Lam43, Lam44;	
    gamma = E1/M1;
    v = sqrt(1.0-1.0/gamma/gamma);
	if(v==0.0)
        return;
    vx = p1x/E1;
	vy = p1y/E1;
	vz = p1z/E1;
	Lam11 = gamma;
	Lam12 = gamma*vx;
	Lam13 = gamma*vy;
	Lam14 = gamma*vz;
	Lam21 = gamma*vx;
	Lam22 = 1+(gamma-1)*vx*vx/v/v;
	Lam23 = (gamma-1)*vx*vy/v/v;
	Lam24 = (gamma-1)*vx*vz/v/v;
	Lam31 = gamma*vy;
	Lam32 = (gamma-1)*vy*vx/v/v;
	Lam33 = 1+(gamma-1)*vy*vy/v/v;
	Lam34 = (gamma-1)*vy*vz/v/v;
	Lam41 = gamma*vz;
	Lam42 = (gamma-1)*vz*vx/v/v;
	Lam43 = (gamma-1)*vz*vy/v/v;
	Lam44 = 1+(gamma-1)*vz*vz/v/v;
	E3   = Lam11*E2+Lam12*p2x+Lam13*p2y+Lam14*p2z;	
	p3x  = Lam21*E2+Lam22*p2x+Lam23*p2y+Lam24*p2z;	
	p3y  = Lam31*E2+Lam32*p2x+Lam33*p2y+Lam34*p2z;	
	p3z  = Lam41*E2+Lam42*p2x+Lam43*p2y+Lam44*p2z;
		
	FourMomentum(p3x, p3y, p3z, E3);
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

double Particle::Momentum(){
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
    ostr << name << " " << E << " " << px << " " << py << " " << pz << " " << m << " " << origin_coords[0] << " " << origin_coords[1] << " " << origin_coords[2] << " " << origin_coords[3] << " " << end_coords[0] << " " << end_coords[1] << " " << end_coords[2] << " " << end_coords[3] << std::endl;
} 
