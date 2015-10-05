#include "minimization.h"
#include "DMNscattering.h"
#include <iostream>
#include <functional>

using std::cout;
using std::endl;
using namespace std::placeholders;

int main(){
	std::function<double(double)> fp = bind(dsigmadEdmN,1,_1,0.03,0.1,0.1,1e-3);
	cout << lim_func_wrapper(0.7, 0, fp, Efmin(1,0.03,0.938),1) << endl;
	std::function<double(double)> fplim = bind(lim_func_wrapper,_1,0,fp,Efmin(1,0.03,0.939),1);
	cout << fplim(0.7) << endl;
	cout << fplim(1.1) << endl;
	double a=0.7;
	double b=0.8;
	double c=mnbrak(a,b,fplim);
	cout << a << " " << b << " " << c << endl;
	double xmin=0;
	double fmax = -1.0*golden(a,b,c,fplim,1e-10,1e-20,xmin);
	cout << xmin << " " << fmax << endl;
}
