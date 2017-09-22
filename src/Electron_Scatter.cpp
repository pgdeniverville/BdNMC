#include "Scatter.h"
#include "DMscattering.h"
#include <cmath>
#include <stdlib.h>
#include <iostream>
#include "Random.h"
#include "Kinematics.h"
#include <algorithm>
#include "constants.h"
#include <list>


using std::list;
using std::cout;
using std::endl;
const double Me = MASS_ELECTRON;
const double Pi = pi;
const double convmcm = 100.0;

void scatterevent (double MDP, double kap, double alD, Particle DM, Particle &electron);

Electron_Scatter::Electron_Scatter(double MDM, double MV, double alphaprime, double kappa, double eEmax, double eEmin){
	Escatmax = eEmax+Me;
	Escatmin = eEmin+Me;
	set_Model_Parameters(MDM, MV, alphaprime, kappa);
}

bool Electron_Scatter::probscatter (std::shared_ptr<detector>& det, Particle &DM){
	
    double LXdet = det->Ldet(DM)*convmcm;
	double prob=sigma2(DM.E,DM.m,MDP,kap,alD,scatmax(DM.E, DM.m),scatmin(DM.E,DM.m))*(det->ENtot())*LXdet*convGeV2cm2;
	if(prob > Random::Flat(0,1)*pMax){
		if (prob > pMax)
			pMax = prob;
        return true;
    }
    return false;
}

bool Electron_Scatter::probscatter (std::shared_ptr<detector>& det, Particle &DM, Particle &electron){
	double prob;
	
    double LXdet = det->Ldet(DM)*convmcm;
	if((prob = sigma2(DM.E,DM.m,MDP,kap,alD,scatmax(DM.E,DM.m),scatmin(DM.E,DM.m))*(det->ENtot())*LXdet*convGeV2cm2) > Random::Flat(0,1)*pMax){
		if (prob > pMax)
			pMax = prob;
        scatterevent(DM, electron);
        return true;
    }
    return false;
}

bool Electron_Scatter::probscatter(std::shared_ptr<detector>& det, list<Particle>& partlist, list<Particle>::iterator& DMit){
	Particle electron(0);
	if(probscatter(det, *DMit, electron)&&(min_angle<=0||electron.Theta()>min_angle)&&(max_angle>2*pi||electron.Theta()<max_angle)){
		DMit->Generate_Position();
        Link_Particles(*DMit, electron);
		partlist.insert(std::next(DMit),electron);
        Particle DMout(DMit->m);
        DMout.name = "Recoil_DM";
        DMout.ThreeMomentum(DMit->px-electron.px,DMit->py-electron.py,DMit->pz-electron.pz);
        Link_Particles(*DMit, DMout);
        partlist.insert(std::next(DMit),DMout);
		return true;
	}
	return false;
}


double Electron_Scatter::scatmax(double DME, double DMM){
	return std::min(Escatmax,EeTMax(DME,DMM));
}

double Electron_Scatter::scatmin(double DME, double DMM){
	return std::max(Escatmin,EeTMin(DME,DMM));
}

void Electron_Scatter::scatterevent (Particle &DM, Particle &electron) {
	double EeMin = scatmin(DM.E,DM.m);
    double EeMax = scatmax(DM.E,DM.m);
	double dsigmax = dsigmadEe(EeMin,DM.E,DM.m,MDP,kap,alD);
    double xe, thetae, phie, pe;

    while(true){
        xe = Random::Flat(0,1)*(EeMax-EeMin)+EeMin;
        if(dsigmadEe(xe,DM.E,DM.m,MDP,kap,alD)>dsigmax*Random::Flat(0,1)){
		electron.Set_Mass(Me);
		electron.name = "electron";
            	thetae = ThetaEe(xe,DM.E,DM.m);
            	phie = Random::Flat(0,1)*2*Pi;
            	pe = sqrt(xe*xe-pow(electron.m,2));
		electron.ThreeMomentum(pe*sin(thetae)*cos(phie),pe*sin(thetae)*sin(phie),pe*cos(thetae));
		electron.Rotate_y(DM.Theta()); 
		electron.Rotate_z(DM.Phi());
            
		break;
        }
    }    
}
/*
double quickscatter (double EDM, double MDM, double MDP, double kap, double alD){
   return sigma(EDM, MDM, MDP, kap, alD)*convGeV2cm2*convmcm;
}

double quickscatter2(double EDM, double MDM, double MDP, double kappa, double alphaD, double Emax, double Emin)
{
    return sigma2(EDM, MDM, MDP, kappa, alphaD, Emax, Emin)*convGeV2cm2*convmcm;
}
*/
