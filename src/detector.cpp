#include "detector.h"
#include <cmath>
#include <vector>
#include <iostream>

using std::vector;
using std::cout;
using std::endl;

//inner product for 3-vectors
double ip(double vec1[],double vec2[]){
    return vec1[0]*vec2[0]+vec1[1]*vec2[1]+vec1[2]*vec2[2];
}

void detector::add_material(double nd, double np, double nn, double ne, double m, std::string matname){
    p_num_tot+=np*nd;
    n_num_tot+=nn*nd;
    e_num_tot+=ne*nd;
    matvec.push_back(Material(nd, np, nn, ne, m, matname));
}


/*********************
 * detector_sphere
 * ******************/

detector_sphere::detector_sphere (double x, double y, double z, double radius){
    r[0] = x; r[1] = y; r[2] = z; Rdet = radius;
}

// distance DM travels through detector
double detector_sphere::Ldet (const Particle &DM) {
	double Ldetenter, Ldetexit;
    double A, B, C;	
//cout << DM.name << " " << DM.E << " " << DM.px << " " << DM.py << " " << DM.pz << " " << DM.m << " " << DM.origin_coords[0] << " " << DM.origin_coords[1] << " " << DM.origin_coords[2] << " " << DM.origin_coords[3] << " " << DM.end_coords[0] << " " << DM.end_coords[1] << " " << DM.end_coords[2] << " " << DM.end_coords[3] << std::endl;
    b[0]=DM.px;b[1]=DM.py;b[2]=DM.pz;
	o[0]=r[0]-DM.origin_coords[0];o[1]=r[1]-DM.origin_coords[1];o[2]=r[2]-DM.origin_coords[2];
	/*
	for(int i=0; i<3; i++)
		cout << " r" << i << "=" << r[i];
	
	for(int i=0; i<3; i++)
		cout << " DM.origin_coords" << i << "=" << DM.origin_coords[i];


	for(int i=0; i<3; i++)
		cout << " o" << i << "=" << o[i];
	cout << endl;	
	*/
    A = ip(b,b);
   // B = -2*(ip(r,b)-ip(o,b));
   // C = (ip(r,r)-2*ip(o,r)+ip(o,o))-pow(Rdet,2);
    B = -2*ip(o,b);
    C = ip(o,o)-pow(Rdet,2);

	/*
	cout << "A = " << A << " B = " << B << " C = " << C << endl;
    */
	
	if((pow(B,2)-4*A*C)<0||A==0){
        return 0.0;
    }
	
	//DM.report(std::cout);

    Ldetenter = (-B-sqrt(pow(B,2)-4*A*C))/(2*A); 
	Ldetexit = (-B+sqrt(pow(B,2)-4*A*C))/(2*A);
	
	//std::cout << Ldetenter << " " << Ldetexit << std::endl;
	if(Ldetenter<0)
		Ldetenter=0;

	if(Ldetexit<0)
		Ldetexit=0;

    if(Ldetenter<Ldetexit){
        cross_point[0] = Ldetenter;
        cross_point[1] = Ldetexit;
        return (Ldetexit-Ldetenter)*sqrt(A);
    }
    else
        return 0;
}

/*******************
 *detector_cylinder*
 *******************/


//May want to add an extra constructor that instead asks for the coordinates of the two faces.
detector_cylinder::detector_cylinder (double x, double y, double z, double detlength, double radius, double detTheta, double detPhi){
    //r points to the center of the detector
    r[0]=x;
    r[1]=y;
    r[2]=z;
    Rdet=radius;
    Ldetector=detlength;

    //l points from the center of the detector to the center of the one of the circular detector faces.
    l[0] = detlength*cos(detPhi)*sin(detTheta)/2;
    l[1] = detlength*sin(detPhi)*sin(detTheta)/2;
    l[2] = detlength*cos(detTheta)/2;
    
}
//This is a temporary setup.
/*
double detector_cylinder::Ldeto (const Particle &DM, const double offset[3]){
    b[0]=DM.px;b[1]=DM.py;b[2]=DM.pz;
    double ro[3];
    for(int i=0; i<3; i++){
        ro[i] = r[i] - offset[i];
    } 
    vector<double> crossings;
    double B1,B2;
    

    //Checking crossing point of circular faces.
    if(ip(b,l)!=0){//if b.l==0, beam is parallel to the circular faces.
        

        B1 = (ip(ro,l)+ip(l,l))/ip(b,l);
        double Rhold[3];//Need to hold this value for use in calculation.
        if(B1>0){
            for(int iter=0; iter <3; iter++){
                Rhold[iter] = B1*b[iter]-l[iter]-ro[iter];
            }
            if(sqrt(ip(Rhold,Rhold))<=Rdet){
                crossings.push_back(B1);
            }
        }
            

        B2 = (ip(ro,l)-ip(l,l))/ip(b,l);//For B2 l->-l.
        if(B2>0){ 
            
            for(int iter=0; iter <3; iter++){
                Rhold[iter] = B2*b[iter]+l[iter]-ro[iter];
            }
            
            if(sqrt(ip(Rhold,Rhold))<=Rdet){
                crossings.push_back(B2);
            }
        }
    }
    

    double B3, B4, A3, A4;

    if(crossings.size()<2){
        double X = -pow(ip(b,l),2) + ip(b,b)*ip(l,l);
        double Y = -2*ip(ro,b)*ip(l,l)+2*ip(b,l)*ip(ro,l);
        double Z = -pow(ip(ro,l),2)-pow(Rdet,2)*ip(l,l)+ip(ro,ro)*ip(l,l);
        double sqr_arg=Y*Y-4*X*Z; //argument in the quadratic equation square root.
        if(sqr_arg > 0){
            B3 = (-Y+sqrt(sqr_arg))/(2*X);
            B4 = (-Y-sqrt(sqr_arg))/(2*X);

            
            if(B3 > 0 && (A3 = (B3*ip(b,l)-ip(ro,l))/ip(l,l) ) > -1 && A3<1){
                crossings.push_back(B3);
            }
             
            if(B4 > 0 && (A4 = (B4*ip(b,l)-ip(ro,l))/ip(l,l) )>-1 && A4<1){
                crossings.push_back(B4);
            }
        }
    }

    
    if(crossings.size()==0)
        return 0.0;
    else if(crossings.size()==2){
        cross_point[0]=crossings[0];
        cross_point[1]=crossings[1];
        return (crossings[1]>crossings[0] ? (crossings[1]-crossings[0])*sqrt(ip(b,b)) :\
                (crossings[0]-crossings[1])*sqrt(ip(b,b)));
    }
    else if(crossings.size()==1)
        //return crossings[0]*sqrt(ip(b,b));
        //needs to be written better for inside of detector case.
        return 0.0;
    else
        throw crossings.size();
}
*/

//I should write a 3 vector class, and overload multiplication and addition.
//This can also be improved by allowing the intersection points to be negative, but
//setting them to zero is so. That means the intersection happened in the past, and we either
//originating in the detector, or going in the opposite direction. Both cases are handled
//properly by the Ldet calculation.
double detector_cylinder::Ldet (const Particle &DM){
    b[0]=DM.px;b[1]=DM.py;b[2]=DM.pz;
   
	for(int i=0; i<3; i++){
        o[i] = r[i] - DM.origin_coords[i];
    }

    vector<double> crossings;
    double B1,B2;

    //Checking crossing point of circular faces.
    if(ip(b,l)!=0){//if b.l==0, beam is parallel to the circular faces.
        

        B1 = (ip(o,l)+ip(l,l))/ip(b,l);
        double Rhold[3];//Need to hold this value for use in calculation.
        if(B1>0){
            for(int iter=0; iter <3; iter++){
                Rhold[iter] = B1*b[iter]-l[iter]-o[iter];
            }
            if(sqrt(ip(Rhold,Rhold))<=Rdet){
                crossings.push_back(B1);
            }
        }
            

        B2 = (ip(o,l)-ip(l,l))/ip(b,l);//For B2 l->-l.
        if(B2>0){ 
            
            for(int iter=0; iter <3; iter++){
                Rhold[iter] = B2*b[iter]+l[iter]-o[iter];
            }
            
            if(sqrt(ip(Rhold,Rhold))<=Rdet){
                crossings.push_back(B2);
            }
        }
    }
    

    double B3, B4, A3, A4;

    if(crossings.size()<2){
        double X = -pow(ip(b,l),2) + ip(b,b)*ip(l,l);
        double Y = -2*ip(o,b)*ip(l,l)+2*ip(b,l)*ip(o,l);
        double Z = -pow(ip(o,l),2)-pow(Rdet,2)*ip(l,l)+ip(o,o)*ip(l,l);
        double sqr_arg=Y*Y-4*X*Z; //argument in the quadratic equation square root.
        if(sqr_arg > 0){
            B3 = (-Y+sqrt(sqr_arg))/(2*X);
            B4 = (-Y-sqrt(sqr_arg))/(2*X);

            
            if(B3 > 0 && (A3 = (B3*ip(b,l)-ip(o,l))/ip(l,l) ) > -1 && A3<1){
                crossings.push_back(B3);
            }
             
            if(B4 > 0 && (A4 = (B4*ip(b,l)-ip(o,l))/ip(l,l) )>-1 && A4<1){
                crossings.push_back(B4);
            }
        }
    }

    
    if(crossings.size()==0)
        return 0.0;
    else if(crossings.size()==2){
        cross_point[0]=crossings[0];
        cross_point[1]=crossings[1];
        return (crossings[1]>crossings[0] ? (crossings[1]-crossings[0])*sqrt(ip(b,b)) :\
                (crossings[0]-crossings[1])*sqrt(ip(b,b)));
    }
    else if(crossings.size()==1)
        //return crossings[0]*sqrt(ip(b,b));
        //needs to be written better for inside of detector case.
        return 0.0;
    else
        throw crossings.size();
}
