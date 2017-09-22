
#include "Scatter.h"
#include "minimization.h"
#include "DMNscattering.h"
#include "Kinematics.h"
#include <cmath>
#include <iostream>
#include "Random.h"
#include <algorithm>
#include "constants.h"
#include <fstream>
#include <memory>

const double mp = MASS_PROTON;
const double mn = MASS_NEUTRON;
const double convmcm = 100.0;
const double tol_abs=1e-20;
const double tol_frac=1e-10;

using std::list;
using std::cout;
using std::endl;
using namespace std::placeholders;


Inelastic_Nucleon_Scatter::Inelastic_Nucleon_Scatter(double MDM, double MV, double alphaprime, double kappa, const std::string &chan_name, const std::string &filename){
	//std::cout << Edmmax << " " << Edmmin << endl;
	channel_name = chan_name;
	set_Model_Parameters(MDM, MV, alphaprime, kappa);
	load_cross_sections(filename);
}

void Inelastic_Nucleon_Scatter::load_cross_sections(const std::string &filename){
    std::ifstream instream(filename);
	if(!instream.is_open()){
		std::cerr << "Inelastic_Nucleon_Scatter cannot open " << filename << std::endl;
		throw -1;
	}
	
	double temp_eps, temp_alpha_D,scale;
	instream >> temp_eps;
	instream >> temp_alpha_D;

	//scaling to switch from the parameters used to generate the cross sections to the 
	//parameters used by the simulation.
	if(channel_name=="Inelastic_Nucleon_Scattering_Baryonic"){
		scale = pow(alD/temp_alpha_D,2);
	}
	else{
		scale = pow(kap/temp_eps,2)*alD/temp_alpha_D;
	}

	std::vector<double> s_dist;//Protons
	std::vector<double> s_dist_2;//Neutrons
    double in;
    instream >> E_min;
    instream >> in;
    s_dist.push_back(in*scale*femtobarn*pow(convmcm,2));
    instream>>in;
    s_dist_2.push_back(in*scale*femtobarn*pow(convmcm,2));
    E_max=E_min;
    while(instream >> E_max){
        instream >> in;
        //cout << E_max << " " << in << endl;
        s_dist.push_back(in*scale*femtobarn*pow(convmcm,2));
        //cout << E_max << " " << in*scale*femtobarn*pow(convmcm,2) << endl;
        instream >> in;
       	s_dist_2.push_back(in*scale*femtobarn*pow(convmcm,2));
    }
    instream.close();

    scatter_dist = std::unique_ptr<Linear_Interpolation>(new Linear_Interpolation(s_dist, E_min, E_max));
    scatter_dist_n = std::unique_ptr<Linear_Interpolation>(new Linear_Interpolation(s_dist_2, E_min, E_max));
    //for(double i=3.0; i<=300; i+=1.0)
    //	cout << i << " " << scatter_dist->Interpolate(i) << endl;
}


bool Inelastic_Nucleon_Scatter::probscatter(std::shared_ptr<detector>& det, Particle &DM){ 
    using namespace std::placeholders;
    if(DM.E<E_min || DM.E>E_max){
    	return false;
    }
    double LXDet = convmcm*(det->Ldet(DM));
    double XDMp = scatter_dist->Interpolate(DM.E)*det->PNtot();
    double XDMn = scatter_dist_n->Interpolate(DM.E)*(det->NNtot());
    double prob=LXDet*(XDMp+XDMn);
    //cout << prob << " " << LXDet << " " << DM.E << " " << scatter_dist->Interpolate(DM.E) << " " << (det->PNtot()+det->NNtot()) << endl;
	//std::cout << DM.E << " " << XDMp << " " << XDMn << " " << prob << " " << endl;
	if(prob > pMax*Random::Flat(0,1)){
        if(prob > pMax){
        	pMax = prob;
        }
        return true;
    }
    else
        return false;
}

//This only checks for scattering occurence, no cuts are possible other than those imposed
//when the inelastic cross section data was generated.
bool Inelastic_Nucleon_Scatter::probscatter(std::shared_ptr<detector>& det, list<Particle>& partlist, list<Particle>::iterator& DMit){
	Particle nucleon(0);
	if(probscatter(det, *DMit)){
		DMit->Generate_Position();
        Link_Particles(*DMit, nucleon);
		partlist.insert(std::next(DMit),nucleon);
		return true;
	}
	return false;
}
