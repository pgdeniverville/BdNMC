#ifndef GUARD_Random_h
#define GUARD_Random_h

/*  Random.h 
 Random Class  (Aug 2008)
 Returns random numbers from a special distributions of popular algorithms
 */

#include <math.h>
#include <stdlib.h>
#include <time.h>
#include <random>
#include <functional>
#include <chrono>

class Random{
public:
	Random();
	Random(unsigned);
    static double Flat       (double =0.0, double =1.0);
	static double Gauss      (double =0.0, double =1.0);
	static double Exponential(double, double =0.0, double =1.0e100);
	static double BreitWigner(double =0.0, double =1.0);
	static int    Integer    (int =0, int =100);
};
    
#endif
