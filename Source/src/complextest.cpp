#include <iostream>
#include <complex>
#include <cmath>

using std::complex;

double mrh = 0.77, mome = 0.77, mrhop = 1.25, mrhopp = 1.45, momep = 1.25, momepp = 1.45, gammarho= 0.15, gammaome = 0.0085, gammarhop=0.3, gammaomep= 0.3, gammarhopp = 0.5, gammaomepp = 0.5;

double f1ra = 0.6165340033101271; double f1rb = 0.22320420111672623; double f1rc = -0.33973820442685326;
double f1wa = 1.0117544786579074; double f1wb = -0.8816565944110686; double f1wc = 0.3699021157531611;

complex<double> F1r(double q2){
 	return f1ra*mrh*mrh/(mrh*mrh-complex<double>(q2,mrh*gammarho))+f1rb*mrhop*mrhop/(mrhop*mrhop-complex<double>(q2,mrhop*gammarhop))\
 	+f1rc*mrhopp*mrhopp/(mrhopp*mrhopp-complex<double>(q2,mrhopp*gammarhopp));
 } 

complex<double> F1w(double q2){
 	return f1wa*mome*mome/(mome*mome-complex<double>(q2,mome*gammaome))+f1wb*momep*momep/(momep*momep-complex<double>(q2,momep*gammaomep))\
 	+f1wc*momepp*momepp/(momepp*momepp-complex<double>(q2,momepp*gammaomepp));
 }

complex<double> F_1_proton(double q2){
	return F1r(q2)+F1w(q2);
}

int main(){
	for(double i = 0; i<=6.0; i+=0.01)
		std::cout << i << " " << std::real(F_1_proton(i)) << " " << imag(F_1_proton(i)) << std::endl;
}