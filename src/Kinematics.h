#ifndef GUARD_Kinematics_h
#define GUARD_Kinematics_h

// class Kinematics
double pAbs (double,double,double,double);
double pt (double,double,double,double);
double calculatetheta (double,double,double,double);
double phi (double,double,double,double);
double invariantmass (double,double,double,double);
double lambda (double,double,double);
double TriangleFunc(double m1, double m2, double m3);
double TriangleFunc2(double m1, double m2, double m3);
//Turns a beta=speed/c into gamma
double beta_gamma(double beta);
double boost_calc(double mass, double Enegry);
//Calculates the angle between vectors 1 and 2.
double Angle_Spread(double x1, double y1, double z1, double x2, double y2, double z2);

namespace one_to_two_decay{
    double two_body_momentum(double m0, double m1, double m2);
}

namespace elastic_scattering{
    double Theta_from_E2f(double E2f, double E1i, double m1, double m2);
    double E2f_from_Theta(double theta2, double E1i, double m1, double m2);
    double E2fMax(double E1i, double m1, double m2);
    double E2fMin(double E1i, double m1, double m2);
}

namespace annihilation_to_pair{
	//lab frame 1+2->3+4 where m3=m4, but m1!=m2.
	double shat(double E1, double m1, double m2);
	double p1cm(double E1, double m1, double m2);
	double p3cm(double E1, double m1, double m2, double m3);
	double t0(double E1, double m1, double m2, double m3);
	double t1(double E1, double m1, double m2, double m3);
	double E4_from_t(double m2, double m4, double t);
	double E3_from_t(double E1, double m1, double m2, double m3, double t);
	double E3Min(double E1, double m1, double m2, double m3);
	double E3Max(double E1, double m1, double m2, double m3);
	double Theta_from_E3(double E1,double m1,double m2,double E3,double m3);
}

namespace two_to_two_scattering{
    //1+2->3+4.
    //2 is stationary in lab frame.
    double s_lab(double E1lab, double m1, double m2);
    //t=(p2-p4)^2
    double t_lab(double E4, double m2, double m4);
    //E1 is E1lab. This is equation 47.37
    double p1cm(double E1, double m1, double m2);
    double E3cm(double E1, double m1, double m2, double m3, double m4);
    double p3cm(double E1, double m1, double m2, double m3, double m4);
    double t0(double E1lab, double m1, double m2, double m3, double m4);
    double t1(double E1lab, double m1, double m2, double m3, double m4);
    double E4_from_t(double m2, double m4, double t);
    //E1 is in lab frame
    double E3_from_t(double E1, double m1, double m2, double m3, double m4, double t);
    //E1 is in lab frame
    double Theta_from_E4(double E1, double E4, double m1, double m2 , double m3, double m4);
    //Need to check that these are correct, might need to swap t0 and t1. t1 gives E3min, so that should be E4max.
    double E4min(double E1lab, double m1, double m2, double m3, double m4);
    double E4max(double E1lab, double m1, double m2, double m3, double m4);
}
/*
namespace annihilation{
	double E3_from_Theta(double E1, double m1, double theta, double m3);
	double p1cm(double E1, double m1);
	double p3cm(double E1, double m1, double m3);
	double t0(double E1, double m1, double m3);
	double t1(double E1, double m1, double m3);
	double E3_from_t(double E1, double m1, double m3, double t);
	double E3Min(double E1, double m1, double m3);
	double E3Max(double E1, double m1, double m3);
}
*/
#endif
