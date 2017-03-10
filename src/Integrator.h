#ifndef GUARD_Integrator_h
#define GUARD_Integrator_h

#include <stdexcept>
#include <functional>
#include <vector>
#include <algorithm>
double SimpsonsRule(std::function<double(double)>, double min, double max, int steps);
double DoubleExponential(std::function<double(double)>, double min, double max, int N, double stepsize);
double DoubleExponential_adapt(std::function<double(double)> f, double min, double max, int N, double h, double precision);

//From "Numerical Algorithms with C"
//Format is SimpsonCubature(function, xmin, xmax, x_steps/2, ymin, ymax, y_steps/2)
double SimpsonCubature(std::function<double(double, double)> f, double a, double b, int P, double c, double d, int Q);

class Linear_Interpolation{
    public:
        Linear_Interpolation(std::vector<double>, double Xmin, double Xmax);
		Linear_Interpolation(){xmin=0; xmax=0;}
        Linear_Interpolation(const Linear_Interpolation &);
        Linear_Interpolation& operator=(const Linear_Interpolation&);
        double operator()(double x){return Interpolate(x);}
        ~Linear_Interpolation(){};
        double Interpolate(double Xval);
    private:
        std::vector<double> yvals;
        double xmin,xmax,xres;
};

#endif
