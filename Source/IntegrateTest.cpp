#include <functional>
#include <cmath>
#include "Integrator.h"
#include <iostream>

int main(){
	auto func = [](double x) {return (1+pow(x,3))/(pow(0.5-x,2)+1e-10*x);};
	
	std::cout << DoubleExponential_adapt(func, -1, 1, 100, 0.01, 1e-6) << std::endl;
	return 0;
}
