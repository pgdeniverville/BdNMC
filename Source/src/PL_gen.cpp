#include <iostream>
#include <string>
#include <memory>
#include "Kinematics.h"
#include "Distribution.h"
#include "Particle_List.h"
#include "sanfordwang.h"
#include "BMPT_dist.h"
#include "Parameter.h"
using std::cout;
using std::endl;

int main(int argc, char* argv[]){
    if(argc==2)
        parameter_file = argv[1];
    else
        parameter_file = "../parameter.dat";
	
	//-------------------
    // Read in parameters
    //-------------------
    std::ifstream parstream(parameter_file);
    Parameter * par;//should make this a unique_ptr
    try{
    par = new Parameter(parstream);
    } catch(int e_int){
        if(e_int==-1){
            cerr << "Could not open " << parameter_file << ". Terminating run." << endl; 
            parstream.close();      
            return -1;
        }
		if(e_int==-2){
			cerr << "No particle_list_file supplied for particle_list distribution. Terminating run." << endl; 
            parstream.close();      
            return -1;
		}
        else{
            cerr << "Something unexpected went wrong. Terminating run." << endl;
            parstream.close();      
            return -1;
        }
    }
    catch(exception &e){
       cerr << "Exception in " << e.what() << " thrown. Terminating run." << endl;
        parstream.close();      
       return -1; 
    }
    
    if(par->Integrity()!=0){
        cerr << "Missing parameters, terminating run. Please check previous error messages." << endl;
        parstream.close();
        return -1;
    }

    //Initializing Random Number Generator
    if(par->Seed()>=0)
		Random(par->Seed());
	else
		Random();


	
	std::shared_ptr<Distribution> s = std::shared_ptr<Distribution>(new sanfordwang("pi0_sanfordwang"));
	std::shared_ptr<Distribution> s2 = std::shared_ptr<Distribution>(new BMPT(450,96));
	double theta, mom, phi;
	int count = 0;
/*	for(int i=0; i<1e6; i++){
		s->sample_momentum(mom, theta, phi);
		cout << cos(phi)*sin(theta)*mom << " " << sin(phi)*sin(theta)*mom << " " << mom*cos(theta) << " " << sqrt(mom*mom+0.134*0.134) << endl;
	}*/
	for(int i=0; i<1e3; i++){
		s2->sample_momentum(mom,theta,phi);
		if(theta<0.01&&mom>5)
			count++;
	}
	cout <<	(double)count/1e3 << endl;
	return 0;
}
