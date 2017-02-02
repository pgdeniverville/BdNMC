/*
 * NOTES
 *
 *	These minimizaiton algorithms are adapted from  NUMERICAL RECIPES IN FORTRAN 77: THE 
 *	ART OF SCIENTIFIC COMPUTING (ISBN 0-521-43064-X) Copyright (C) 1986-1992 by Cambridge 
 *	University Press.
 *	Programs Copyright (C) 1986-1992 by Numerical Recipes Software. 
 *  
 *
 *  Pointless to use tolerance better than Sqrt(Machine_Precision).
 */

#include "minimization.h"
#include <math.h>
#include <iostream>

using std::cout; using std::endl;

const double GOLD = 1.618034;
const double GLIMIT = 100;
const double TINY = 1e-20;

double mnbrak(double& ax, double& bx, std::function<double(double)> func){
	double r,q,u,ulim,fu;
	double fa = func(ax)+0.0;
	double fb = func(bx)+0.0;
	if(fb>fa){
		double dum = ax;
		ax = bx;
		bx = dum;
		dum = fb;
		fb = fa;
		fa = dum;
	}
	double cx = bx+GOLD*(bx-ax);
	double fc=func(cx)+0.0;
    if(fb==0&&fc==0)
        return cx;
	while(fb>=fc){
		r = (bx - ax)*(fb - fc);
		q = (bx - cx)*(fb - fa);
		u = bx - ((bx - cx)*q - (bx - ax)*r)/(2.0*copysign(fmax(fabs(q - r), TINY),q-r));
		ulim = bx + GLIMIT*(cx - bx);
		if((bx-u)*(u-cx)>0){
			fu = func(u);
			if(fu<fc){
				ax=bx; fa=fb; bx=u; fb=fu;
				return cx;
			}
			else if(fu>fb){
				cx=u;
				fc=fu;
				return cx;
			}
			u=cx+GOLD*(cx-bx);
			fu=func(u);
		}
		else if((cx-u)*(u-ulim)>0){
			fu=func(u);
			if(fu<fc){
				bx=cx; cx=u; u=cx+GOLD*(cx-bx); fb=fc; fc=fu; fu=func(u);
			}
		}
		else if((u-ulim)*(ulim-cx)>=0){
			u=ulim;
			fu=func(u);
		}
		else{
			u=cx+GOLD*(cx-bx);
			fu=func(u);
		}
		ax=bx; bx=cx; cx=u; fa=fb; fb=fc; fc=fu;
	}
	return cx;
}

const double R=0.61803399;
const double C=1.0-R;

/*
 * ax<bx<cx, and f(ax),f(cx)>f(bx)
 * 
 */
double golden(double ax, double bx, double cx, std::function<double(double)> f, double frac_tol, double abs_tol, double& xmin){
	double f1,f2,x0,x1,x2,x3;
	x0=ax;
	x3=cx;
	if(fabs(cx-bx)>fabs(bx-ax)){
		x1=bx;
		x2=bx+C*(cx-bx);
	}
	else{
		x2=bx;
		x1=bx-C*(bx-ax);
	}
	f1=f(x1);
	f2=f(x2);
	while((fabs(x3-x0)>frac_tol*(fabs(x1)+fabs(x2)))&&(fabs(x3-x0)>abs_tol)){
		if(f2<f1){
			x0=x1;
			x1=x2;
			x2=R*x1+C*x3;
			f1=f2;
			f2=f(x2);
		}
		else{
			x3=x2;
			x2=x1;
			x1=R*x2+C*x0;
			f2=f1;
			f1=f(x1);
		}
	}
	if(f1<f2){
		xmin=x1;
		return(f1);
	}
	else{
		xmin=x2;
		return(f2);
	}
	return 0;
}

//Evaluates to def outside of [xmin, xmax]. Flips the sign of func.
double lim_func_wrapper(double x, double def, std::function<double(double)> func, double xmin, double xmax){
	if(x<xmin||x>xmax)
		return def;
	else
		return -1.0*func(x);
}
