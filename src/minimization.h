#ifndef GUARD_Minimization_h
#define GUARD_Minimization_h
#include <functional>

double mnbrak(double&, double&, std::function<double(double)> func);
double golden(double, double, double, std::function<double(double)> func, double frac_tol, double abs_tol, double& xmin);
double lim_func_wrapper(double x, double def, std::function<double(double)> func, double xmin, double xmax);

#endif
