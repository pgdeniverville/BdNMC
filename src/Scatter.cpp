#include "Scatter.h"
#include "minimization.h"
#include "constants.h"

using namespace std::placeholders;

double const tol_frac=1e-10;
double const tol_abs=1e-20;

void Prepare_Cross_Section(std::function<double(double,double)> f_init, std::function<double(double)> ER_Min, std::function<double(double)> ER_Max, std::shared_ptr<Linear_Interpolation> cross, std::shared_ptr<Linear_Interpolation> maxima, double E_in_min, double E_in_max, double E_in_res){
    double a,b,c,xmin;
    std::vector<double> vec_cross;
    std::vector<double> vec_maxima;
    for(double iter=E_in_min+E_in_res; iter<=E_in_max; iter+=E_in_res){
        double E_r_min=ER_Min(iter);
        double E_r_max=ER_Max(iter);
        if(E_r_min>=E_r_max){
            vec_cross.push_back(0);
            vec_maxima.push_back(0);
            continue;
        }
        std::function<double(double)> f = std::bind(f_init,_1,iter);
        std::function<double(double)> flim = std::bind(lim_func_wrapper,_1,0.0,f,E_r_min,E_r_max);

        vec_cross.push_back(DoubleExponential_adapt(f,E_r_min,E_r_max,100,0.1,1e-4));
        a=E_r_min; b=(E_r_max+E_r_min)/2.0;

        c=mnbrak(a,b,flim);
        xmin=0;
        vec_maxima.push_back(-1.0*golden(a,b,c,flim,tol_frac,tol_abs,xmin));
    }

    cross=std::shared_ptr<Linear_Interpolation>(new Linear_Interpolation(vec_cross,E_in_min,E_in_max));
	maxima=std::shared_ptr<Linear_Interpolation>(new Linear_Interpolation(vec_maxima,E_in_min,E_in_max));
}
