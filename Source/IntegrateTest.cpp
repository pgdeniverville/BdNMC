#include <functional>
#include <cmath>
#include "Integrator.h"
#include <iostream>

int main(){
	//auto func = [](double x) {return (1+pow(x,3))/(pow(0.5-x,2)+1e-10*x);};
	//auto func = [](double x, double y) {return (x+y+1.0)*(y+x+1.0);};
	auto func = [](double x, double y) {return (x+y+1)*(x+y+1);};
	//std::cout << DoubleExponential_adapt(func, -1, 1, 100, 0.01, 1e-6) << std::endl;
	std::cout << SimpsonCubature(func,0,2,10,1,4,20) << std::endl;
	return 0;
}
