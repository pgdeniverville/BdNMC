#include "Scatter.h"
#include "minimization.h"
#include "DMNscattering_Baryonic.h"
#include "Kinematics.h"
#include <cmath>
#include <iostream>
#include "Random.h"
#include <algorithm>
#include "constants.h"
#include <list>

using std::list;
using std::cout; using std::endl;

const double mp = MASS_PROTON;
const double mn = MASS_NEUTRON;
const double convmcm = 100.0;
const double tol_abs=1e-20;
const double tol_frac=1e-10;

using namespace std::placeholders;

Nucleon_Scatter_Baryonic::Nucleon_Scatter_Baryonic(double Emini, double Emaxi, double Eresi, double MDM, double MV, double alphaprime, double kappa, double NEmax, double NEmin){
	Edmmin=Emini; Edmmax=Emaxi; Edmres=Eresi; Escatmin=NEmin; Escatmax=NEmax;
	set_Model_Parameters(MDM, MV, alphaprime, kappa);
}

void Nucleon_Scatter_Baryonic::generate_cross_sections(){
	std::vector<double> vec_proton;
	std::vector<double> vec_neutron;
	std::vector<double> vec_proton_maxima;
	std::vector<double> vec_neutron_maxima;
	double a,b,c,xmin,sigmatot;
    for(double iter=Edmmin; iter<=Edmmax; iter+=Edmres){
		//cout << "Energy=" << iter << endl;
		std::function<double(double)> fp = std::bind(DMNscattering_Baryonic::dsigmadEdmP,iter,_1,mdm,MDP,alD,kap);
		std::function<double(double)> fplim = std::bind(lim_func_wrapper,_1,0.0,fp,scatmin(iter,mdm,mp),scatmax(iter));	
		//for(double ii=scatmin(iter,mdm,mp);ii<scatmax(iter);ii+=0.01){
	    //	cout << ii << " " << fp(ii) << " " << endl;
		//}	  
		sigmatot=DoubleExponential_adapt(fp,scatmin(iter,mdm,mp),scatmax(iter),100,0.1,1e-3);
		//Build an array of integrated proton scattering cross sections.
		vec_proton.push_back(sigmatot);
		//cout << fp(iter) << endl;
		if(sigmatot==0){
			vec_proton_maxima.push_back(0);
		}
		else{
			a=scatmin(iter,mdm,mp);
			b=scatmax(iter);
			c=mnbrak(a,b,fplim);
			xmin=0;
			vec_proton_maxima.push_back(-1.0*golden(a,b,c,fplim,tol_frac,tol_abs,xmin));
		}	
		std::function<double(double)> fn = std::bind(DMNscattering_Baryonic::dsigmadEdmN,iter,_1,mdm,MDP,alD,kap);
        std::function<double(double)> fnlim = std::bind(lim_func_wrapper,_1,0.0,fn,scatmin(iter,mdm,mn),scatmax(iter));	
		sigmatot=DoubleExponential_adapt(fn,scatmin(iter,mdm,mn),scatmax(iter),100,0.1,1e-3);
		vec_neutron.push_back(sigmatot);
        if(sigmatot==0){
			vec_neutron_maxima.push_back(0);
		}
		else{
			a=scatmin(iter,mdm,mn);
			b=scatmax(iter);
			c=mnbrak(a,b,fnlim);
			xmin=0;
			vec_neutron_maxima.push_back(-1.0*golden(a,b,c,fnlim,tol_frac,tol_abs,xmin));
		}
		// cout << iter << " " << vec_proton.back() << " " << vec_neutron.back() << endl;
		//cout << iter <<  " " << vec_proton.back() << endl;
		//cout << "scatmin and max = " << scatmin(iter,mdm,mp) << " " << scatmax(iter) << endl;
    }
	proton_cross = std::unique_ptr<Linear_Interpolation>(new Linear_Interpolation(vec_proton,Edmmin,Edmmax));
    neutron_cross = std::unique_ptr<Linear_Interpolation>(new Linear_Interpolation(vec_neutron,Edmmin,Edmmax));
	proton_cross_maxima = std::unique_ptr<Linear_Interpolation>(new Linear_Interpolation(vec_proton_maxima,Edmmin,Edmmax));
	neutron_cross_maxima = std::unique_ptr<Linear_Interpolation>(new Linear_Interpolation(vec_neutron_maxima,Edmmin,Edmmax));
}

bool Nucleon_Scatter_Baryonic::probscatter(std::shared_ptr<detector>& det, list<Particle>& partlist, list<Particle>::iterator& DMit){
	Particle nucleon(0);
	if(probscatter(det, *DMit, nucleon)&&(min_angle<=0||nucleon.Theta()>min_angle)&&(max_angle>2*pi||nucleon.Theta()<max_angle)){
		Generate_Position(det, *DMit, nucleon);
		partlist.insert(std::next(DMit),nucleon);
		return true;
	}
	return false;
}

bool Nucleon_Scatter_Baryonic::probscatter(std::shared_ptr<detector>& det, Particle &DM){ 
    using namespace std::placeholders;

    double LXDet = convmcm*(det->Ldet(DM));
    double XDMp = proton_cross->Interpolate(DM.E)*(det->PNtot());
    double XDMn = neutron_cross->Interpolate(DM.E)*(det->NNtot());
    double prob=LXDet*convGeV2cm2*(XDMp+XDMn);
    if(prob > pMax*Random::Flat(0,1)){
        if(prob > pMax){
        	pMax = prob;
        }
 
		return true;
    }
    else{
        return false;
    }
}

bool Nucleon_Scatter_Baryonic::probscatter(std::shared_ptr<detector>& det, Particle &DM, Particle &Nucleon){ 
    using namespace std::placeholders;

    double LXDet = convmcm*(det->Ldet(DM));
    double XDMp = proton_cross->Interpolate(DM.E)*(det->PNtot());
    double XDMn = neutron_cross->Interpolate(DM.E)*(det->NNtot());
    double prob=LXDet*convGeV2cm2*(XDMp+XDMn);
    //DM.report(cout);
    //cout << "prob = " << prob << " pMax = " << pMax << endl;
	if(prob > pMax*Random::Flat(0,1)){
        if(prob > pMax){
        	pMax = prob;
        }
 
		if(XDMp/(XDMn+XDMp) >= Random::Flat(0,1)){
            std::function<double(double)> Xsec = std::bind(DMNscattering_Baryonic::dsigmadEdmP,DM.E,_1,DM.m,MDP,alD,kap);
            Nucleon.Set_Mass(mp);
            Nucleon.name = "proton";
            scatterevent(DM, Nucleon, Xsec, *proton_cross_maxima);
        }
        else{
            std::function<double(double)> Xsec = std::bind(DMNscattering_Baryonic::dsigmadEdmN,DM.E,_1,DM.m,MDP,alD,kap);
            Nucleon.Set_Mass(mn);
            Nucleon.name = "neutron";
            scatterevent(DM, Nucleon, Xsec, *neutron_cross_maxima);
        }
        return true;
    }
    else
        return false;
}

void Nucleon_Scatter_Baryonic::scatterevent (Particle &DM, Particle &Nucleon, std::function<double(double)> Xsec, Linear_Interpolation& Xmax){
    double EDMMax = scatmax(DM.E);
    double EDMMin = scatmin(DM.E, DM.m, Nucleon.m); 
    double dsigmax = std::max(Xsec(EDMMax),Xmax.Interpolate(DM.E));
    //DM.report(cout);
	//Nucleon.report(cout);
	//cout << "max = " << EDMMax << " min = " << EDMMin << " dsigmax = " << dsigmax << endl;
	double xe,thetaN,phiN,pN;
    while(true){
        xe = Random::Flat(0,1)*(EDMMax-EDMMin)+EDMMin;
        if(Xsec(xe)/dsigmax > Random::Flat(0,1)){
            thetaN = DMNscattering_Baryonic::Ef_to_N_Theta(DM.E,xe,DM.m,Nucleon.m);
            phiN = Random::Flat(0,1)*2*pi;
            pN = sqrt(pow(DM.E+Nucleon.m-xe,2)-pow(Nucleon.m,2));
            Nucleon.ThreeMomentum(pN*sin(thetaN)*cos(phiN),pN*sin(thetaN)*sin(phiN),cos(thetaN)*pN);
            Nucleon.Rotate_y(DM.Theta());
            Nucleon.Rotate_z(DM.Phi());
            break;
        }
    }
}

double Nucleon_Scatter_Baryonic::scatmax(double DME){
	return std::min(DME,DME-Escatmin);
}

double Nucleon_Scatter_Baryonic::scatmin(double DME, double DMM, double NMM){
	return std::max(DMNscattering_Baryonic::Efmin(DME, DMM, NMM), DME-Escatmax);
}
