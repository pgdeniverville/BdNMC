#include <iostream>
#include <string>
#include <memory>
#include "Kinematics.h"
#include "Distribution.h"
#include "Particle_List.h"
#include "sanfordwang.h"
#include "BMPT_dist.h"
#include "DMgenerator.h"
using std::cout;
using std::endl;

int main(){
/*	std::shared_ptr<Distribution> s = std::shared_ptr<Distribution>(new sanfordwang("pi0_sanfordwang"));
	std::shared_ptr<Distribution> s2 = std::shared_ptr<Distribution>(new BMPT(450,96));
	double theta, mom, phi;
	int count = 0;
	for(int i=0; i<1e6; i++){
		s->sample_momentum(mom, theta, phi);
		cout << cos(phi)*sin(theta)*mom << " " << sin(phi)*sin(theta)*mom << " " << mom*cos(theta) << " " << sqrt(mom*mom+0.134*0.134) << endl;
	}*/
/*	for(int i=0; i<1e3; i++){
		s2->sample_momentum(mom,theta,phi);
		if(theta<0.01&&mom>5)
			count++;
	}
	cout <<	(double)count/1e3 << endl;
*/
	parton_V_gen parV(1,0.1,1e-3,0.1,"parton_production");
	cout << parV.Channel_Name() << endl;
	cout << parV.BranchingRatio() << endl;
	pion_decay_gen pid(0.1,0.01,1e-3,0.1);
	cout << pid.Channel_Name() << endl;
	return 0;
}
