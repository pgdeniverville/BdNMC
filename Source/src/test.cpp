//#include "sanfordwang.h"
//#include "BMPT_dist.h"
//#include "BurmanSmith.h"
#include "Random.h"
#include "Proton_Brem_Distribution.h"

#include <iostream>
#include <functional>

using std::cout;
using std::endl;
using namespace std::placeholders;

int main(){
	Random();
<<<<<<< HEAD
	Proton_Brem_Distribution test(8,1e-3,0.1,0.2,0.7,0.3,0.0);
	//sanfordwang sw("pi0_sanfordwang");
	//BMPT bmpt(400,184);
	//BurmanSmith bs(0.8,6);
	//double angle[] = {0.05, 0.15, 0.3};
/*	for(double i = 0; i<=0.8; i+=0.01){
=======
	//sanfordwang sw("pi0_sanfordwang");
	//BMPT bmpt(400,184);
	BurmanSmith bs(0.8,6);
	double angle[] = {0.05, 0.15, 0.3};
	for(double i = 0; i<=0.8; i+=0.01){
>>>>>>> inelastic-pion-dev
//	cout << i << " " << 0.5*(sw.swpip(i,angle[0])+sw.swpim(i,angle[0])) << " " << 0.5*(sw.swpip(i,angle[1])+sw.swpim(i,angle[1]))<< " " << 0.5*(sw.swpip(i,angle[2])+sw.swpim(i,angle[2])) << endl;
		cout << i;
	   	for(int j = 0; j<sizeof(angle)/sizeof(double); j++){
			//cout << " " << 0.5*(bmpt.Invariant_Cross_Section(i,angle[j])+bmpt.Invariant_Cross_Section_pi_minus(i,angle[j]));
			cout << " " << 0.5*(bs.dsigma(i,angle[j]*180.0/3.14159));
	   	}
		cout << endl;

	}
<<<<<<< HEAD
*/	//	 << " " << 0.5*(bmpt.Invariant_Cross_Section(i,angle[1])+bmpt.Invariant_Cross_Section_pi_minus(i,angle[1]))<< " " << 0.5*(bmpt.Invariant_Cross_Section(i,angle[2])+bmpt.Invariant_Cross_Section_pi_minus(i,angle[2])) << endl;
=======
	//	 << " " << 0.5*(bmpt.Invariant_Cross_Section(i,angle[1])+bmpt.Invariant_Cross_Section_pi_minus(i,angle[1]))<< " " << 0.5*(bmpt.Invariant_Cross_Section(i,angle[2])+bmpt.Invariant_Cross_Section_pi_minus(i,angle[2])) << endl;
>>>>>>> inelastic-pion-dev
//	}

/*
	for(int i = 0; i<1e6; i++){
		Particle pion(0.134);
		//sw.sample_particle(pion);
		bs.sample_particle(pion);
		cout << pion.Momentum() << " "  <<  pion.Theta() << endl;
	}
*/
	return 0;
}
