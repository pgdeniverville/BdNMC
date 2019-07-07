#include "Integrator.h"
#include "minimization.h"
#include "Random.h"

#include <cmath>
#include <stdexcept>
#include <iostream>
#include "constants.h"

using std::cout; using std::endl;
using std::cerr;
using namespace std::placeholders;
const int DDMAX = 20;

using std::vector;

double ddweight(double k, double h){
    return 1.0/2.0*h*pi*cosh(k*h)/pow(cosh(1.0/2.0*pi*sinh(k*h)),2);
}

long double ddabscissa(int k, double h){
    return tanh(1.0/2.0*pi*sinh(h*k));
}
//rescales an abscissa on integration region [min,max] to [-1,1]
long double res_int(double x, double min, double max){
    return 1.0/2.0*(min+max-min*x+max*x);
}

double func_check(std::function<double(double)> f, double abscissa, double min, double max){
	double fp,fm;
	fp=f(res_int(abscissa,min,max));
	//cout << "fp=" << fp << endl;
	if(std::isnan(fp))
		fp=0;
	fm=f(res_int(-abscissa,min,max));
	//cout << "fm = " << fm<< endl;
	if(std::isnan(fm))
		fm=0;
	return fm+fp;
}

const int _N_MAX = 1e9;

double RandomIntegrate2_adapt(std::function<double(double,double)> f, double xmin, double xmax, double ymin, double ymax, int n_start, double r_accuracy_goal){
    int n = n_start;
    double estimate = RandomIntegrate2(f, xmin, xmax, ymin, ymax, n);
    double estimate_2;
    while(n<_N_MAX){
        n*=2;
        estimate_2 = RandomIntegrate2(f, xmin, xmax, ymin, ymax, n);
        //cout << n << " " << estimate << " " << estimate_2 << endl;
        if((abs(estimate-estimate_2))/(0.5*(estimate+estimate_2))<r_accuracy_goal){
            return (estimate+2*estimate_2)/3.0;
        }
        else{
            estimate = (estimate+2*estimate_2)/3.0;
            n*=1.5;
        }
    }
    cerr << "RandomIntegrate2_adapt reached max steps without achieving accuracy goal\n";
    cerr << "Estimate = " << estimate << " next step = " << estimate_2 << endl;
    return estimate_2;
}

//Unfortunately this retraces its steps.
double SimpsonCubature_adapt(std::function<double(double,double)> f, double xmin, double xmax, int xsteps_start, double ymin, double ymax, int ysteps_start, double r_accuracy_goal){
    int xsteps = xsteps_start;
    int ysteps = ysteps_start;
    double estimate = SimpsonCubature(f, xmin, xmax, xsteps, ymin, ymax, ysteps);
    double estimate_2;
    while(xsteps < _N_MAX && ysteps < _N_MAX){
        xsteps*=2;
        ysteps*=2;
        estimate_2=SimpsonCubature(f, xmin, xmax, xsteps, ymin, ymax, ysteps);
        cout << estimate << " " << estimate_2 << endl;
        if((abs(estimate-estimate_2))/(0.5*(estimate+estimate_2))<r_accuracy_goal){
            return estimate_2;
        }
        else{
            estimate = estimate_2;
        }
    }
    cerr << "SimpsonCubature_adapt reached max steps without achieving accuracy goal\n";
    cerr << "Estimate = " << estimate << " next step = " << estimate_2 << endl;
    return estimate_2;
}

double RandomIntegrate2(std::function<double(double,double)> f, double xmin, double xmax, double ymin, double ymax, int n){
    double tot=0;
    for(int i = 0; i < n; i++){
        tot+=f(Random::Flat(xmin,xmax),Random::Flat(ymin,ymax));
    }
    return tot/(double)n*(xmax-xmin)*(ymax-ymin);
}

//From "Numerical Algorithms with C"
//Format is SimpsonCubature(function, xmin, xmax, x_steps/2, ymin, ymax, y_steps/2)
double SimpsonCubature(std::function<double(double, double)> f, double a, double b, int P, double c, double d, int Q){
	double hx = (b-a)/(2.0*P);
	double hy = (d-c)/(2.0*Q);
	double hold = f(a,c)+f(b,c)+f(a,d)+f(b,d);
	double hold2=0;
	for(int p = 0; p<=P-1; p++)
		hold2+=f(a+(2*p+1)*hx,c)+f(a+(2*p+1)*hx,d);
	hold+=4*hold2;
	hold2=0;
	for(int q = 0; q<=Q-1; q++)
		hold2+=f(a,c+(2*q+1)*hy)+f(b,c+(2*q+1)*hy);
	hold+=4*hold2;
	hold2=0;
	for(int p = 1; p<=P-1; p++)
		hold2+=f(a+2*p*hx,c)+f(a+2*p*hx,d);
	hold+=2*hold2;
	hold2=0;
	for(int q = 1; q<=Q-1; q++)
		hold2+=f(a,c+2*q*hy)+f(b,c+2*q*hy);
	hold+=2*hold2;
	hold2=0;
	for(int p = 1; p<=P-1; p++)
		for(int q = 1; q<=Q-1; q++)
			hold2+=f(a+2*p*hx,c+2*q*hy);
	hold+=4*hold2;
	hold2=0;
	for(int p = 0; p<=P-1; p++)
		for(int q = 0; q<=Q-1; q++)
			hold2+=f(a+(2*p+1)*hx,c+(2*q+1)*hy);
	hold+=16*hold2;
	hold2=0;
	for(int p = 1; p<=P-1; p++)
		for(int q = 0; q<=Q-1; q++)
			hold2+=f(a+2*p*hx,c+(2*q+1)*hy);
	hold+=8*hold2;
	hold2=0;
	for(int p = 0; p<=P-1; p++)
		for(int q = 1; q<=Q-1; q++)
			hold2+=f(a+(2*p+1)*hx,c+2*q*hy);
	hold+=8*hold2;
	hold2=0;
	return hold*hx*hy/9.0;
}


//h is the width between abscissa.
double DoubleExponential(std::function<double(double)> f, double min, double max, int N, double h){
    double scale = 2.0/(max - min);
    if(min>=max)
        return 0;
    double sum = 0; double abscissa = 0;
    //cout << "scale=" << scale << endl;
	for(int k=0; k<=N; k++){
        if((abscissa=ddabscissa(k,h))>=1)
            break;
		//cout << " k=" << k << " ddweight=" << ddweight(k,h) << " f+=" << f(res_int(abscissa,min,max)) << " res_int+=" << res_int(+abscissa,min,max) << " f-=" << f(res_int(-abscissa,min,max)) <<  " res_int-=" << res_int(-abscissa,min,max) << endl;
        //sum+=ddweight(k,h)*(f(res_int(abscissa,min,max))+f(res_int(-abscissa,min,max)));
        sum+=ddweight(k,h)*func_check(f,abscissa,min,max);
		//cout << "cumulative sum " << sum << endl;
    }
	//cout << "double exp sum = " << sum << " double exp scale " << scale << endl;
    return sum/scale;
}

//Double exponential routine used for extending range in n.
double DoubleExponential_Nout(std::function<double(double)> f, double min, double max, int N, int Nmin, double h){
    double scale = 2.0/(max - min);
    if(min>=max)
        return 0;
    double sum = 0; double abscissa = 0;
    for(int k=Nmin+1; k<=N; k++){
        if((abscissa=ddabscissa(k,h))>=1)
            break;
        //sum+=ddweight(k,h)*(f(res_int(abscissa,min,max))+f(res_int(-abscissa,min,max)));
        sum+=ddweight(k,h)*func_check(f,abscissa,min,max);
    }
    return sum/scale;
}

//Only sums over odd k. Used to increase resolution.
double DoubleExponential_hdub(std::function<double(double)> f, double min, double max, int N, double h){
    double scale = 2.0/(max - min);
    if(min>=max)
        return 0;
    double sum = 0; double abscissa = 0;
    //cout << "hdub\n";
	for(int k=1; k<=N-1; k+=2){
        if((abscissa=ddabscissa(k,h))>=1)
            break;
        //sum+=ddweight(k,h)*(f(res_int(abscissa,min,max))+f(res_int(-abscissa,min,max)));
        //cout << "k=" << k << " ddweight = " << ddweight(k,h) << endl;
		sum+=ddweight(k,h)*func_check(f,abscissa,min,max);
		//cout << "cumulative sum= " << sum << endl;
    }
	//cout << "hdub sum = " << sum << " hdub scale = " << scale  << endl;
    return sum/scale;
}
//slow and inefficient. Still good enough. Could be improved by adding a means of subdividing the integration region.
double DoubleExponential_adapt(std::function<double(double)> f, double min, double max, int N, double h, double precision){
    //cout << "initial estimate" << endl;
	double dd0 = DoubleExponential(f, min, max, N, h);
    double dd1,dd2;
    int attempts;
	if(min>=max)
		return 0;
    //cout << "N=" << N << " h=" << h << " dd0=" << dd0 << endl;
    for(attempts=0; attempts<DDMAX; attempts++){
	   	//cout << "attempts=" << attempts << endl;	
        if(fabs(dd0-(dd1 = dd0+DoubleExponential_Nout(f, min, max, N*2.0, N, h)))>fabs(dd0*precision)){
            //cout << fabs(dd0-(dd1 = dd0+DoubleExponential_Nout(f, min, max, N*2.0, N, h))) << endl;
            dd0=dd1;
            N*=2.0;
            //cout << "N=" << N << " h=" << h << " dd1=" << dd0 << endl;
        }
        else if(fabs(dd0-(dd2 = 0.5*dd0+DoubleExponential_hdub(f, min, max, 2*N, h/2.0)))>fabs(dd0*precision)){
            dd0=dd2;
            h*=0.5; N*=2.0;
            //cout << "N=" << N << " h=" << h << " dd2=" << dd0 << endl;
			//cout << "double exp comparison\n";
			//double tmp = DoubleExponential(f, min, max, N, h);
			//cout << "double exp out = " << tmp << endl;
        }
        else
            break;
    }
    if(attempts==DDMAX){
        std::cerr << "Integrator DoubleExponential_adapt reached max number of iterations=" << DDMAX << endl;
        std::cerr << "Estimate=" << dd0 << " " << dd1 << " " << dd2 << endl; 
	}
    return dd0;
}

//Not good enough anymore
double SimpsonsRule(std::function<double(double)> f, double min, double max, int steps){
    if(min>=max)
        return 0.0;
    if(steps%2!=0){
        std::cerr << "Invalid Argument: steps must be a multiple of 2\n";
    }
    double dx = (double)(max-min)/steps;
    double sum1 = 0;
    double sum2 = 0;

    for(int i=1;i<=steps;i+=2){
        sum1 += 4.0 * f(min+i*dx);
    }
    for(int i=2;i<=steps-1;i+=2){
        sum2 += 2.0 * f(min+i*dx);
    }
    return (sum1+sum2+f(min)+f(max)) * dx/3.0;
}

//This works on CDFs only.
bool Interpolation1D::Find_Y(double y, double &x){
    if(!IS_CDF){
        cerr << "Is not a CDF\n";
        return false;
    }
    if(Interpolate(xmin)>Interpolate(xmax)){
        if(Interpolate(xmin)<y||Interpolate(xmax)>y){
            cerr << "xmin= " << Interpolate(xmin) << "<" << y << " or xmax=" << Interpolate(xmax) << ">" << y << endl; 
            return false;
        }
    }
    else{
        if(Interpolate(xmin)>y||Interpolate(xmax)<y){
            cerr << "ymin= " << Interpolate(xmin) << ">" << y << " or xmax=" << Interpolate(xmax) << "<" << y << endl;
            return false;
        }
    }
    std::function<double(double)> f = std::bind(&Interpolation1D::Dif,this,_1,y);
    double a = xmin;
    double c = xmax;
    double b = (xmax+xmin)/2.0;
    golden(a,b,c,f,1e-3,std::abs(xmin-xmax)/100.0,x);
    return true;
}

/*
 *  Linear_Interpolation takes a set of data Yvals whose elements map to xvals between Xmin
 *  and Xmax such that Yvals[0] = f(Xmin), Yvals[i] = f(Xmin + i*xres), where xres is defined
 *  in the function body and f(x) is some function approximating the data set Yvals.
 */
Linear_Interpolation::Linear_Interpolation(vector<double> Yvals, double Xmin, double Xmax){
    yvals.resize(Yvals.size());
    std::copy(Yvals.begin(),Yvals.end(),yvals.begin());
    xmin = Xmin;
    xmax = Xmax;
    xres = (Xmax-Xmin)/(Yvals.size()-1);
    yvals.push_back(0);
}

/*
 * Interpolate returns the value of f(xval) for xmin<=xval<=xmax. Values outside this range will
 * likely lead to a segmenation fault. Should perhaps add error handling, not certain of the 
 * cost at run time to check the domain at every evaluation.
 *
 * It has to be done somewhere, added error checking. Returns 0.0 rather
 * than trying to extrapolate.
 */
double Linear_Interpolation::Interpolate(double xval){
    if(xval<xmin)
        return 0.0;
    else if(xval>xmax)
        return 0.0;
    double index = (xval-xmin)/xres;
    int lowindex = (int)floor(index);
    return (lowindex+1-index)*yvals[lowindex]+(index-lowindex)*yvals[lowindex+1];
}

bool Linear_Interpolation::QCDF(){
    double sign = 1;
    if(yvals[0]>yvals[1]){
        sign = -1;
    }

    for(vector<double>::iterator it = yvals.begin()+1; it!=yvals.end(); it++){
        //This is DICEY
        if(sign*(*it)<(sign)*(*(it-1))){
            return false;
        }
    }
    return true;
}

void Linear_Interpolation::Convert_to_CDF(){
    for(vector<double>::iterator it = yvals.begin()+1; it!=yvals.end(); it++){
        *it += *(it-1);
    }
    IS_CDF=true;
}

//Convert f(x) = y to f(y) = x. This will only work if the interpolation is monotonically increasing (or decreasing).
void Linear_Interpolation::Invert(){
    double ymin,ymax;
    if(yvals[0]>yvals[1]){
        ymax=Interpolate(xmin);
        ymin=Interpolate(xmax);
    }
    else{
        ymin=Interpolate(xmin);
        ymax=Interpolate(xmax);
    }

    if(!QCDF()){
        std::cerr << "Cannot invert this function, not monotonically increasing or decreasing." << endl;
        //Is there some form of range error?
        throw -1;
    }

    vector<double> xvals;
    double n = (xmax-xmin)/xres; 
    double yres=(ymax-ymin)/n;
    double x;
    for(double i=ymin;i<=ymax;i+=yres){
        if(Find_Y(i,x)){
            xvals.push_back(x);
        }
        else{
            cerr << "Invert() could not find expected y_val." << endl;
            throw std::domain_error("yval not in range");
        }
    }
    xmin = ymin;
    xmax = ymax;
    xres = yres;
    yvals = xvals;
}

Linear_Interpolation::Linear_Interpolation(const Linear_Interpolation &LI){
    yvals = LI.yvals;
    xmin = LI.xmin;
    xmax = LI.xmax;
    xres = LI.xres;
}

Linear_Interpolation& Linear_Interpolation::operator=(const Linear_Interpolation& LI){
    yvals = LI.yvals;
    xmin = LI.xmin;
    xmax = LI.xmax;
    xres = LI.xres;
    return *this;
}

void Linear_Interpolation::Report(){
    cout << "Interpolation report: Linear_Interpolation" << endl;
    cout << "x_min=" << xmin << " x_max=" << xmax << endl;
    for(unsigned it=0; it!= yvals.size(); it++){
        cout << it*xres+xmin << " " << yvals[it] << endl;
    }
}

/*
 *  Linear_Interpolation2 takes a set of data Yvals whose elements map to a set of x values Xvals, such that Y[i] = f(Xvals[i]). 
 *  It is assumed that Xvals[i]<Xvals[i+1].
 */
Linear_Interpolation2::Linear_Interpolation2(vector<double> Xvals, vector<double> Yvals){
    if(Yvals.size()!=Xvals.size()){
        cerr << Yvals.size() << "!=" << Xvals.size() << endl;
        throw std::domain_error("Linear_Interpolation2 requires the length of Yvals to equal the length of Xvals.");
    }
    yvals.resize(Yvals.size());
    std::copy(Yvals.begin(),Yvals.end(),yvals.begin());
    xvals.resize(Xvals.size());
    std::copy(Xvals.begin(),Xvals.end(),xvals.begin());
    xmin = Xvals.front();
    xmax = Xvals.back();
}

//Finds the index i corresponding to the largest x[i]<=xval.
double Find_Index(const double xval, const vector<double>& xvals){
    unsigned L = 0;
    unsigned R = xvals.size()-1;
    unsigned m;
    while(R-L>1){ 
        m=(unsigned)std::floor((L+R)/2);
        //cout << "Find_Index L=" << L << " R=" << R << " m" << m << " xvals[m]=" << xvals[m] << endl;
        if(xvals[m]>xval){
            R=m-1;
        }
        else if(xvals[m]<xval){
            L=m+1;
        }
        else{
            return m;
        }
    }
    return L;
}

/*
 * Interpolate returns the value of f(xval) for xmin<=xval<=xmax. Values outside this range will
 * evaluate to zero.
 * It has to be done somewhere, added error checking. Returns 0.0 rather
 * than trying to extrapolate.
 */
double Linear_Interpolation2::Interpolate(const double xval){
    if(xval<xmin)
        return 0.0;
    else if(xval>xmax)
        return 0.0;
/*
    if(xmin==xval){
        return yvals[0];
    }
    else if(xmax == xval){
        return yvals[-1];
    }*/
    //cout << "Begin Interpolation for xval=" << xval << endl;
    double lowindex = Find_Index(xval, xvals);
    //cout << "Found lowindex=" << lowindex << " " << xvals[lowindex] << " " << xvals[lowindex+1] << endl;
    //cout << "Returning " << yvals[lowindex] << " " << yvals[lowindex+1] << " " << ((xvals[lowindex+1]-xval)*yvals[lowindex]+(xval-xvals[lowindex])*yvals[lowindex+1])/(xvals[lowindex+1]-xvals[lowindex]) << endl;
    //throw -1;
    if(xvals[lowindex]==xval){
        return yvals[lowindex];
    }
    return ((xvals[lowindex+1]-xval)*yvals[lowindex]+(xval-xvals[lowindex])*yvals[lowindex+1])/(xvals[lowindex+1]-xvals[lowindex]);
}

bool Linear_Interpolation2::QCDF(){
    double sign = 1;
    if(yvals[0]>yvals[1]){
        sign = -1;
    }

    for(vector<double>::iterator it = yvals.begin()+1; it!=yvals.end(); it++){
        //This is DICEY
        if(sign*(*it)<(sign)*(*(it-1))){
            return false;
        }
    }
    return true;
}

void Linear_Interpolation2::Convert_to_CDF(){
    for(vector<double>::iterator it = yvals.begin()+1; it!=yvals.end(); it++){
        *it += *(it-1);
    }
    IS_CDF=true;
}

//Convert f(x) = y to f(y) = x. This will only work if the interpolation is monotonically increasing (or decreasing).
void Linear_Interpolation2::Invert(){
    if(!QCDF()){
        std::cerr << "Cannot invert this function, not monotonically increasing or decreasing." << endl;
        //Is there some form of range error?
        throw -1;
    }
//    double ymin,ymax;
    vector<double> tmpvals;
    if(yvals[0]>yvals[1]){
        tmpvals.resize(yvals.size());
        std::reverse_copy(yvals.begin(), yvals.end(), tmpvals.begin());
    }
    else{
        tmpvals = yvals;
    }
    yvals = xvals;
    xvals = tmpvals;
    xmin = xvals.front();
    xmax = xvals.back();
}

Linear_Interpolation2::Linear_Interpolation2(const Linear_Interpolation2 &LI){
    yvals = LI.yvals;
    xvals = LI.xvals;
    xmin = LI.xmin;
    xmax = LI.xmax;
}

Linear_Interpolation2& Linear_Interpolation2::operator=(const Linear_Interpolation2& LI){
    yvals = LI.yvals;
    xvals = LI.xvals;
    xmin = LI.xmin;
    xmax = LI.xmax;
    return *this;
}

void Linear_Interpolation2::Report(){
    cout << "Interpolation report: Linear_Interpolation" << endl;
    cout << "x_min=" << xmin << " x_max=" << xmax << endl;
    for(unsigned it=0; it!= yvals.size(); it++){
        cout << xvals[it] << " " << yvals[it] << endl;
    }
}
