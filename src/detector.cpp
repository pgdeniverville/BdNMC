#include "detector.h"
#include <cmath>
#include <vector>
#include <iostream>

using std::vector;
using std::cout;
using std::endl;

//inner product for 3-vectors
double ip(const double vec1[], const double vec2[]){
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
    B = -2*ip(o,b);
    C = ip(o,o)-pow(Rdet,2);
	
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

    if(Ldetexit>Ldetenter){
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

//I should write a 3 vector class, and overload multiplication and addition.
//This can also be improved by allowing the intersection points to be negative. 
//That means the intersection happened in the past, and we are either
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

/*******************
 *detector_cuboid*
 *******************/

//In place rotation of 3-vector arr.
void XYX_rotate(const double t1, const double t2, const double t3, double *arr){
    double x=arr[0];
    double y=arr[1];
    double z=arr[2];
    //cout << arr[0] << " " << arr[1] << " " << arr[2] << endl;
    arr[0] = x*cos(t2)+z*cos(t3)*sin(t2)+y*sin(t2)*sin(t3);
    arr[1] = x*sin(t1)*sin(t2)+z*(-cos(t2)*cos(t3)*sin(t1)-cos(t1)*sin(t3))+y*(cos(t1)*cos(t3)-cos(t2)*sin(t1)*sin(t3));
    arr[2] = -x*cos(t1)*sin(t2)+y*(cos(t3)*sin(t1)+cos(t1)*cos(t2)*sin(t3))+z*(cos(t1)*cos(t2)*cos(t3)-sin(t1)*sin(t3));
    //cout << arr[0] << " " << arr[1] << " " << arr[2] << endl;
}

detector_cuboid::detector_cuboid(double x, double y, double z, double detwidth, double detheight, double detlength, double detPhi, double detTheta, double detPsi){
    r[0]=x; r[1]=y; r[2]=z;
    //Length oriented along z-axis
    face[0][0]=detwidth/2.0;
    face[1][1]=detheight/2.0;
    face[2][2]=detlength/2.0;
    //dim[0]=pow(detwidth/2.0,2);
    //dim[1]=pow(detheight/2.0,2);
    //dim[2]=pow(detlength/2.0,2);
    for(int i = 0; i<3; i++){
        XYX_rotate(detPhi,detTheta,detPsi,face[i]);
        for(int j=0; j<3; j++)
            face[i+3][j]=-face[i][j];
    }
    for(int i=0; i<6; i++){
        //cout << face[i][0] << " " <<face[i][1] << " " <<face[i][2] << endl;
    }

}

double detector_cuboid::Ldet (const Particle &DM){
    b[0]=DM.px;b[1]=DM.py;b[2]=DM.pz;
     
    double o[3];
	for(int i=0; i<3; i++){
        o[i] = r[i] - DM.origin_coords[i];
    }
    double entry=0;
    double exit=0;
    //cout << "DETECTOR_START\n";
    //DM.report(cout);
    bool consistency=false;
    for(int i=0; i<6; i++){
        double bot = ip(b,face[i]);
        //cout << "bot=" << bot << endl;
        if(bot==0){
            continue;
        }
        double ft = ip(o,face[i])+ip(face[i],face[i]);
        //cout << "ft" << i << "= " << ft << endl;
        //cout << face[i][0] << " " << face[i][1] << " " << face[i][2] << endl;
        //cout << ip(o,face[i]) << endl;
        //cout << ip(face[i],face[i]) << endl;
        if(ft==0){
            if(consistency)//Two faces are parallel to R=r+n.
            {
                std::cerr << "Only one face of cuboid detector can be parallel to o.";
                throw -1;
            }
            else
                consistency=true;
            continue;
        }
        double detcross=ft/bot;
        //cout << "detcross=" << detcross << endl;
        if(ft<0&&(detcross>entry||entry==0.0)){
            entry=detcross;
            //cout << "entry" << endl;
        }
        else if(detcross>0){
            if(exit==0||detcross<exit){
                exit=detcross;
                //cout << "exit" << endl;
            }
        }
        
    }
    if(entry>=exit||entry<0)
        return 0.0;
    else{
        //DM.report(cout);
        //cout << "CrossReport " << entry << " " << exit << " " << (exit-entry)*sqrt(ip(b,b)) << endl;
        cross_point[0]=entry;
        cross_point[1]=exit;
        return (exit-entry)*sqrt(ip(b,b));
    }
}


