#ifndef GUARD_Integrator_h
#define GUARD_Integrator_h

#include <stdexcept>
#include <functional>
#include <vector>
#include <algorithm>
//I should add an adaptive Simpson's Rule.
double SimpsonsRule(std::function<double(double)>, double min, double max, int steps);
double DoubleExponential(std::function<double(double)>, double min, double max, int N, double stepsize);
double DoubleExponential_adapt(std::function<double(double)> f, double min, double max, int N, double h, double precision);

class Linear_Interpolation{
    public:
        Linear_Interpolation(std::vector<double>, double Xmin, double Xmax);
		~Linear_Interpolation(){};
        double Interpolate(double Xval);
    private:
        std::vector<double> yvals;
        double xmin, xres;
};

#endif
