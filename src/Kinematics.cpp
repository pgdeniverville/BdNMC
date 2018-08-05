#include "Kinematics.h"
#include <math.h>
#include <iostream>

namespace elastic_scattering{
    //This namespace deals with 2->2 elastic scattering
    
    //Returns recoil angle of particle 2 in the lab frame.
    //E2f: Recoil energy of particle 2
    //E1i: Incident energy of particle 1
    //m1,2: masses of particle 1 and 2
    double Theta_from_E2f(double E2f, double E1i, double m1, double m2){
        return acos(sqrt((E2f-m2)/(E2f+m2))*(E1i+m2)/sqrt(E1i*E1i-m1*m1));
    }
    //This is the inverse of Theta_from_E2f.
    double E2f_from_Theta(double theta2, double E1i, double m1, double m2){
        return m2*((E1i - m1)*(E1i + m1)*pow(cos(theta2),2) + pow((E1i + m2),2))/(pow(E1i + m2,2) + (-pow(E1i,2) + pow(m1,2))*pow(cos(theta2),2));
    }
    //Maximum and minimum energy of outgoing recoil particle
    //as function of masses and incoming energy E1i
    double E2fMax(double E1i, double m1, double m2){
        return E2f_from_Theta(0.0, E1i, m1, m2);
    }
    double E2fMin(double E1i, double m1, double m2){
        return E2f_from_Theta(M_PI/2.0, E1i, m1, m2);
    }
}
/*
namespace annihilation{
	//Center of mass frame annihilation, 1+2->3+4.

	//s as calculated with lab frame variables
	double shat_lab(double E1, double m1, double m2){
		return m1*m1+m2*m2+2*E1*m2;
	}

	//s as calculated with center of mass frame Energy
	double shat_cm(double Ecm){
		return 4*Ecm*Ecm;
	}
}
*/
namespace annihilation_to_pair{
	//lab frame 1+2->3+4 where m3=m4, but m1!=m2.
	double shat(double E1, double m1, double m2){
		return m1*m1+m2*m2+2*E1*m2;
	}

	double t_hat(double E1, double m1, double m2, double E3, double m3){
		return m3*m3-m2*m2+2*m2*(E3-E1);
	}

	double p1cm(double E1, double m1, double m2){
		return sqrt(E1*E1-m1*m1)*m2/sqrt(shat(E1, m1, m2));
	}

	double p3cm(double E1, double m1, double m2, double m3){
		return sqrt(shat(E1, m1, m2)/4.0-m3*m3);
	}

	double t0(double E1, double m1, double m2, double m3){
		return pow(m1*m1-m2*m2,2)/4.0/shat(E1,m1,m2)-pow(p1cm(E1,m1,m2)-p3cm(E1,m1,m2,m3),2);
	}

	double t1(double E1, double m1, double m2, double m3){
		return pow(m1*m1-m2*m2,2)/4.0/shat(E1,m1,m2)-pow(p1cm(E1,m1,m2)+p3cm(E1,m1,m2,m3),2);
	}

	double E4_from_t(double m2, double m4, double t){
		return (m2*m2+m4*m4-t)/(2*m2);
	}

	double E3_from_t(double E1, double m1, double m2, double m3, double t){
		return m2+E1-E4_from_t(m2, m3, t);
	}

/*	double Theta_from_E3(double E1,double m1,double m2,double E3,double m3){
		return acos((2*E1*E3+2*(E1-E3)*m2-m1*m1-m2*m2)/2/sqrt(E1*E1-m1*m1)/sqrt(E3*E3-m3*m3));
	}*/

	double Theta_from_E3(double E1, double m1, double m2, double E3, double m3){
		double bot = (2*sqrt(E1*E1-m1*m1)*sqrt(E3*E3-m3*m3));
		if(bot==0){
			return M_PI;
		}
		return acos((t_hat(E1,m1,m2,E3,m3)-m1*m1-m3*m3+2*E1*E3)/bot);
	}

	double E3Min(double E1, double m1, double m2, double m3){
		return E3_from_t(E1, m1, m2, m3, t1(E1,m1,m2,m3));
	}

	double E3Max(double E1, double m1, double m2, double m3){
		return E3_from_t(E1, m1, m2, m3, t0(E1,m1,m2,m3));
	}
}

//
double pAbs (double px, double py, double pz, double p0) {
	double rpAbs;
	rpAbs = sqrt(px*px+py*py+pz*pz);
	return(rpAbs);
}
//
double pt (double px, double py, double pz, double p0) {
	double rpt;
	rpt = sqrt(px*px+py*py);
	return(rpt);
}
//
double calculatetheta (double px, double py, double pz, double p0) {
	double rtheta;
	rtheta = acos(pz/pAbs(px,py,pz,p0));
	return(rtheta);
}
//
double phi (double px, double py, double pz, double p0) {
    if(px==0 && py==0){
        return 0;
    }
	if (px > 0 && py > 0) 
	{
		return atan(fabs(py/px));
	}
	else if (px < 0 && py > 0) 
	{
		return M_PI-atan(fabs(py/px));
	}	
	else if (px < 0 && py < 0) 
	{
		return M_PI+atan(fabs(py/px));
	}	
	else
	{
		return 2*M_PI - atan(fabs(py/px));
	}	
}
//
double invariantmass (double px, double py, double pz, double p0) {
	double rinvariantmass;
	rinvariantmass = sqrt(p0*p0-px*px-py*py-pz*pz);
	return(rinvariantmass);
}
//
double lambda (double a, double b, double c) {
	double rlam;
	rlam = a*a+b*b+c*c-2*a*b-2*a*c-2*b*c;
	return(fabs(rlam));
}
//
double TriangleFunc(double m1, double m2, double m3){
    return 1/(2*m1)*sqrt(pow(m1,4)+pow(m2,4)+pow(m3,4)-2*pow(m1*m2,2)-2*pow(m1*m3,2)-2*pow(m2*m3,2));
}

double TriangleFunc2(double m1, double m2, double m3){
	return pow(m1,4)+pow(m2,4)+pow(m3,4)-2*pow(m1*m2,2)-2*pow(m1*m3,2)-2*pow(m2*m3,2);
}

double beta_gamma(double beta){
    return 1.0/sqrt(1-pow(beta,2));
}

double Angle_Spread(double x1, double y1, double z1, double x2, double y2, double z2){
    return acos((x1*x2+y1*y2+z1*z2)/sqrt(x1*x1+y1*y1+z1*z1)/sqrt(x2*x2+y2*y2+z2*z2));
}

//Kinematics for DM+X->DM+X scattering, with X at rest.

// X Angle as a function of X final energy, DM initial energy, and masses.
//double ThetaEe (double Ee, double EDM, double MDM, MX) {
//	return(acos(sqrt((Ee-MX)/(Ee+MX))*(EDM+MX)/sqrt(EDM*EDM-MDM*MDM)));
//}
