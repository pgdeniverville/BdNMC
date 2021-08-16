#ifndef GUARD_Integrator_h
#define GUARD_Integrator_h

#include <stdexcept>
#include <functional>
#include <vector>
#include <algorithm>
#include <memory>
#include <cmath>

double SimpsonsRule(std::function<double(double)>, double min, double max, int steps);
double DoubleExponential(std::function<double(double)>, double min, double max, int N, double stepsize);
double DoubleExponential_adapt(std::function<double(double)> f, double min, double max, int N, double h, double precision);


double SimpsonCubature_adapt(std::function<double(double,double)> f, double xmin, double xmax, int xsteps_start, double ymin, double ymax, int ysteps_start, double r_accuracy_goal=0.01);
//From "Numerical Algorithms with C"
//Format is SimpsonCubature(function, xmin, xmax, x_steps/2, ymin, ymax, y_steps/2)
double SimpsonCubature(std::function<double(double, double)> f, double a, double b, int P, double c, double d, int Q);

double RandomIntegrate2_adapt(std::function<double(double,double)> f, double xmin, double xmax, double ymin, double ymax, int n_start, double r_accuracy_goal=0.01);
double RandomIntegrate2(std::function<double(double,double)> f, double xmin, double xmax, double ymin, double ymax, int n);

class Interpolation1D{
    public:
        Interpolation1D(){IS_CDF=false;}
        virtual ~Interpolation1D(){};
        virtual double Interpolate(double Xval) = 0;
        double operator()(const double x){return Interpolate(x);}
        double Max(){return xmax;}
        double Min(){return xmin;}
        //This is used in search algorithms. 
        double Dif(double x, double val){return std::abs(this->Interpolate(x)-val);}
        //This only runs for IS_CDF=true
        bool Find_Y(double Y, double &x);
        void Report();
        //Is this a strictly increasing or decreasing function?
        bool QCDF(){return IS_CDF;}

        //bool CDF(Interpolation1D& interp){return false;}
    protected:
        double xmin, xmax;
        bool IS_CDF;
};

class Linear_Interpolation : public Interpolation1D{
    public:
        Linear_Interpolation(std::vector<double>, double Xmin, double Xmax);
		Linear_Interpolation(){xmin=0; xmax=0;}
        Linear_Interpolation(const Linear_Interpolation &);
        Linear_Interpolation& operator=(const Linear_Interpolation&);
        //double operator()(double x){return Interpolate(x);}
        ~Linear_Interpolation(){};
        double Interpolate(const double Xval);
        //Converts the interpolated function in a cumulative function.
        //Should only be used on functions f(x)>0 for all x in domain [xmin, xmax].
        void Report();
        bool QCDF();
        void Convert_to_CDF();
        void Invert();
    private:
        std::vector<double> yvals;
        double xres;
};

//This is a more general version of Linear_Interpolation.
class Linear_Interpolation2 : public Interpolation1D{
    public:
        Linear_Interpolation2(std::vector<double> xvals, std::vector<double> yvals);
        Linear_Interpolation2(){xmin=0; xmax=0;}
        Linear_Interpolation2(const Linear_Interpolation2 &);
        Linear_Interpolation2& operator=(const Linear_Interpolation2&);
        //double operator()(double x){return Interpolate(x);}
        ~Linear_Interpolation2(){};
        double Interpolate(const double Xval);
        //Converts the interpolated function in a cumulative function.
        //Should only be used on functions f(x)>0 for all x in domain [xmin, xmax].
        void Report();
        bool QCDF();
        void Convert_to_CDF();
        void Invert();
    private:
        std::vector<double> yvals;
        std::vector<double> xvals;
};

#endif
