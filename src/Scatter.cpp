#include "Scatter.h"
#include "minimization.h"
#include "constants.h"

using namespace std::placeholders;

using std::cout; using std::endl;

double const tol_frac=1e-10;
double const tol_abs=1e-20;

//f_init(E1,E4) for 1+2->3+4 scattering.
void Prepare_Cross_Section(std::function<double(double,double)> f_init, std::function<double(double)> ER_Min, std::function<double(double)> ER_Max, std::shared_ptr<Linear_Interpolation>& cross, std::shared_ptr<Linear_Interpolation>& maxima, double E_in_min, double E_in_max, double E_in_res){
    double a,b,c,xmin;
    std::vector<double> vec_cross;
    std::vector<double> vec_maxima;
//    cout << "E_in_res=" << E_in_res << endl;
    for(double iter=E_in_min+E_in_res; iter<=E_in_max; iter+=E_in_res){

//        cout << "E=" << iter << endl;
        double E_r_min=ER_Min(iter);
        double E_r_max=ER_Max(iter);
//        cout << iter << " " << E_r_min << " " << E_r_max << endl;
        
        if(E_r_min>=E_r_max){
//            cout << "E_r_min > E_r_max? cross = 0\n";
            vec_cross.push_back(0);
            vec_maxima.push_back(0);
            continue;
        }
        std::function<double(double)> f = std::bind(f_init,iter,_1);
        std::function<double(double)> flim = std::bind(lim_func_wrapper,_1,0.0,f,E_r_min,E_r_max);
       
//        cout << "f(" << iter << "," << (E_r_min+E_r_max)/2.0 << ") = " << f((E_r_min+E_r_max)/2.0) << endl;

        vec_cross.push_back(DoubleExponential_adapt(f,E_r_min,E_r_max,200,0.1,1e-4));

//        cout << "cross = " << vec_cross.back() << endl;

        a=E_r_min; b=(E_r_max+E_r_min)/2.0;

        c=mnbrak(a,b,flim);
        xmin=0;
        vec_maxima.push_back(-1.0*golden(a,b,c,flim,tol_frac,tol_abs,xmin));
    }
    cross=std::shared_ptr<Linear_Interpolation>(new Linear_Interpolation(vec_cross,E_in_min,E_in_max));
	maxima=std::shared_ptr<Linear_Interpolation>(new Linear_Interpolation(vec_maxima,E_in_min,E_in_max));
}
