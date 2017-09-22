#include "Scatter.h"
#include "minimization.h"
#include "DMNscattering.h"
#include "DMNscattering_Baryonic.h"
#include "Kinematics.h"
#include <cmath>
#include <iostream>
#include "Random.h"
#include <algorithm>
#include "constants.h"


const double mp = MASS_PROTON;
const double mn = MASS_NEUTRON;
const double convmcm = 100.0;
const double tol_abs=1e-20;
const double tol_frac=1e-10;

using std::list;
using std::cout;
using std::endl;
using std::vector;
using namespace std::placeholders;

//NEmin and NEmax are kinetic energies!
Nucleon_Scatter::Nucleon_Scatter(double Emini, double Emaxi, double Eresi, double MDM, double MV, double alphaprime, double kappa, double NEmax, double NEmin, std::string model, const bool cohere, std::shared_ptr<detector> det){
	Edmmin=std::max(Emini,MDM); Edmmax=Emaxi; Edmres=Eresi;
	Escatmin=NEmin;
	Escatmax=NEmax;

    //If things start returning zero, I made a scoping error right here.
    if(model == "Baryonic_V"){
        dsig_prot = &DMNscattering_Baryonic::dsigmadEdmP;
        dsig_neut = &DMNscattering_Baryonic::dsigmadEdmN;
        dsig_cohp = &DMNscattering_Baryonic::dsigmadEdmP_coherent;
    }
    else{
        dsig_prot = &dsigmadEdmP;
        dsig_neut = &dsigmadEdmN;
        dsig_cohp = &dsigmadEdmP_coherent;
    }

    //cout << "Escatmin= " << Escatmin << " Escatmax= " << Escatmax << endl; 
	min_angle=-1;//Set min_angle to less than zero so it always passes cut.
	max_angle=2*pi+1;//Set max_angle to larger than 2*pi so it always passes cut.
    coherent=cohere;
    if(coherent){
        if(det==NULL){
            std::cerr << "Nucleon_Scatter requires a detector ptr to operate in coherent mode.\n";
            throw -31;
        }
        set_Model_Parameters(MDM, MV, alphaprime, kappa, det);
    }
    else{
	    set_Model_Parameters(MDM, MV, alphaprime, kappa);
    }
}

void Nucleon_Scatter::cross_gen_handler(std::function<double(double)> fp, std::function<double(double)> fplim, std::vector<double> &cross_vec, std::vector<double> &cross_vec_maxima, double iter, double m1, const double m2){ 
    
    double sigmatot=DoubleExponential_adapt(fp,scatmin(iter,m1,m2),scatmax(iter),100,0.1,1e-3);
	
    //Build an array of integrated scattering cross sections.
	cross_vec.push_back(sigmatot);
	
    if(sigmatot==0){
        //If sigma is zero it is unlikely to have a maximum.
	    cross_vec_maxima.push_back(0);
	}
	else{
        double a=scatmin(iter,m1,m2);
        double b=scatmax(iter);
        double c=mnbrak(a,b,fplim);
        double xmin=0;
        cross_vec_maxima.push_back(-1.0*golden(a,b,c,fplim,tol_frac,tol_abs,xmin));
    }
}


void Nucleon_Scatter::generate_cross_sections(){
	std::vector<double> vec_proton_maxima;
	std::vector<double> vec_neutron_maxima;
	std::vector<double> vec_proton;
	std::vector<double> vec_neutron;

//	double emax; double emin;
    for(double iter=Edmmin; iter<=Edmmax; iter+=Edmres){
		
        //Proton
        std::function<double(double)> fp = std::bind(dsig_prot,iter,_1,mdm,MDP,alD,kap);
		std::function<double(double)> fplim = std::bind(lim_func_wrapper,_1,0.0,fp,scatmin(iter,mdm,mp),scatmax(iter));	
	    cross_gen_handler(fp, fplim, vec_proton, vec_proton_maxima,iter,mdm,mp);
        //Neutron
		std::function<double(double)> fn = std::bind(dsig_neut,iter,_1,mdm,MDP,alD,kap);
        std::function<double(double)> fnlim = std::bind(lim_func_wrapper,_1,0.0,fn,scatmin(iter,mdm,mn),scatmax(iter));	
	    cross_gen_handler(fn, fnlim, vec_neutron, vec_neutron_maxima,iter,mdm,mn);
		//cout << iter << " " << vec_proton.back() << endl; 
		//cout << "scatmin and max = " << scatmin(iter,mdm,mp) << " " << scatmax(iter) << endl;
    }
    //throw -1;
	proton_cross = Linear_Interpolation(vec_proton,Edmmin,Edmmax);
    neutron_cross = Linear_Interpolation(vec_neutron,Edmmin,Edmmax);
	proton_cross_maxima = Linear_Interpolation(vec_proton_maxima,Edmmin,Edmmax);
	neutron_cross_maxima = Linear_Interpolation(vec_neutron_maxima,Edmmin,Edmmax);
}

void Nucleon_Scatter::generate_coherent_cross_sections(std::shared_ptr<detector>& det){
    cout << "Generating Coherent Cross Sections\n";
    for(unsigned i =0; i!=det->mat_num(); i++){
        std::vector<double> vec_atom_maxima;
        std::vector<double> vec_atom;
        double mass=det->M(i);
        for(double iter=Edmmin; iter<=Edmmax; iter+=Edmres){
		    std::function<double(double)> fa = std::bind(dsig_cohp,iter,_1,mdm,MDP,alD,kap,det->PN(i)+det->NN(i),det->PN(i));
		    std::function<double(double)> falim = std::bind(lim_func_wrapper,_1,0.0,fa,scatmin(iter,mdm,mass),scatmax(iter));	
	        cross_gen_handler(fa, falim, vec_atom, vec_atom_maxima,iter,mdm,mass);
        }
        atom_maxima.push_back(Linear_Interpolation(vec_atom_maxima,Edmmin,Edmmax));
        atom_cross.push_back(Linear_Interpolation(vec_atom,Edmmin,Edmmax));
    }
}

bool Nucleon_Scatter::probscatter(std::shared_ptr<detector>& det, Particle &DM){ 
    //cout << "Calling probscatter properly" << endl;
    
    double LXDet = convmcm*(det->Ldet(DM));
    double XDM=0;
    
    if(coherent){
       //cout << "Coherent\n";
        for(unsigned i=0; i<det->mat_num(); i++){
           //cout << "Iterating\n";
            XDM +=atom_cross[i].Interpolate(DM.E)
                *(det->get_nDensity(i));
            //cout << "atom_cross:" << atom_cross[i].Interpolate(DM.E) << " DM.E:" << DM.E << endl;
        }
    }//coherent
    else{
         XDM = proton_cross.Interpolate(DM.E)*(det->PNtot())
             +neutron_cross.Interpolate(DM.E)*(det->NNtot());
    }//incoherent

    double prob=LXDet*convGeV2cm2*XDM;
    //cout << prob << endl;
	if(prob > pMax*Random::Flat(0,1)){
        if(prob > pMax){
        	pMax = prob;
        }
 
		return true;
    }
    else
        return false;
}


bool Nucleon_Scatter::probscatter(std::shared_ptr<detector>& det, list<Particle>& partlist, list<Particle>::iterator& DMit){
	Particle nucleon(0);
	if(probscatter(det, *DMit, nucleon)&&(min_angle<=0||nucleon.Theta()>min_angle)&&(max_angle>2*pi||nucleon.Theta()<max_angle)){
		DMit->Generate_Position();
        Link_Particles(*DMit, nucleon);
		partlist.insert(std::next(DMit),nucleon);
        Particle DMout(DMit->m);
        DMout.name = "Recoil_DM";
        DMout.ThreeMomentum(DMit->px-nucleon.px,DMit->py-nucleon.py,DMit->pz-nucleon.pz);
        Link_Particles(*DMit, DMout);
        partlist.insert(std::next(DMit),DMout);
		return true;
	}
	return false;
}
//Supplied particle may no longer be a Nucleon!
bool Nucleon_Scatter::probscatter(std::shared_ptr<detector>& det, Particle &DM, Particle &Nucleon){ 
	using namespace std::placeholders;
    double LXDet = convmcm*(det->Ldet(DM));
    if(coherent){    
        vector<double> XDMv;
        double XDM =0;
        //cout << "DM.E=" << DM.E << endl; 
        for(unsigned int i=0; i<atom_cross.size(); i++){
            XDMv.push_back(atom_cross[i].Interpolate(DM.E)*(det->get_nDensity(i)));
            XDM+=XDMv.back();
            //cout << "XDM=" << XDM << endl;
        }
        
        double prob=LXDet*convGeV2cm2*XDM;
        //cout << "prob=" << prob << endl; 
        if(prob>pMax*Random::Flat(0,1)){ 
            if(prob>pMax){
                pMax=prob;
            }
            
            unsigned int i;
            double u = Random::Flat(0,1);
            for(i=0;i<atom_cross.size();i++){
                if(u<XDMv[i]/XDM){
                    break;
                }
                u-=XDMv[i]/XDM;
            }
            
            std::function<double(double)> Xsec = 
            std::bind(dsig_cohp,DM.E,_1,DM.m,MDP,alD,
                kap,det->PN(i)+det->NN(i),
                det->PN(i));
            
            Nucleon.Set_Mass(det->M(i));
            Nucleon.name = det->matname(i);

            if(scatmax(DM.E)<scatmin(DM.E, DM.m,Nucleon.m))
                return false;
            scatterevent(DM,Nucleon,Xsec,atom_maxima[i]);
            return true;
        }
            
   }

    else{
        double XDMp = proton_cross.Interpolate(DM.E)*(det->PNtot());
        double XDMn = neutron_cross.Interpolate(DM.E)*(det->NNtot());
        double prob=LXDet*convGeV2cm2*(XDMp+XDMn);
        if(prob > pMax*Random::Flat(0,1)){
            if(prob > pMax)
                pMax = prob;
            if(XDMp/(XDMn+XDMp) >= Random::Flat(0,1)){
                std::function<double(double)> Xsec = std::bind(dsig_prot,DM.E,_1,DM.m,MDP,alD,kap);
                Nucleon.Set_Mass(mp);
                Nucleon.name = "proton";
                if(scatmax(DM.E) < scatmin(DM.E, DM.m, Nucleon.m))
                    return false;
                scatterevent(DM, Nucleon, Xsec, proton_cross_maxima);
            }
            else{
                std::function<double(double)> Xsec = std::bind(dsig_neut,DM.E,_1,DM.m,MDP,alD,kap);
                Nucleon.Set_Mass(mn);
                Nucleon.name = "neutron";
                if(scatmax(DM.E) < scatmin(DM.E, DM.m, Nucleon.m))
                    return false;
                scatterevent(DM, Nucleon, Xsec, neutron_cross_maxima);
            }
            return true;
        }
    }
    
    return false;
}

void Nucleon_Scatter::scatterevent (Particle &DM, Particle &Nucleon, std::function<double(double)> Xsec, Linear_Interpolation& Xmax){
    double EDMMax = scatmax(DM.E);
	double EDMMin = scatmin(DM.E, DM.m, Nucleon.m); 
	double dsigmax = std::max(Xsec(EDMMax),Xmax.Interpolate(DM.E));
	double xe,thetaN,phiN,pN;
    while(true){
        xe = Random::Flat(0,1)*(EDMMax-EDMMin)+EDMMin;
        if(Xsec(xe)/dsigmax > Random::Flat(0,1)){
            thetaN = Ef_to_N_Theta(DM.E,xe,DM.m,Nucleon.m);
            phiN = Random::Flat(0,1)*2*pi;
            pN = sqrt(pow(DM.E+Nucleon.m-xe,2)-pow(Nucleon.m,2));
            Nucleon.ThreeMomentum(pN*sin(thetaN)*cos(phiN),pN*sin(thetaN)*sin(phiN),cos(thetaN)*pN);
            Nucleon.Rotate_y(DM.Theta());
            Nucleon.Rotate_z(DM.Phi());
            break;
        }
    }
}

double Nucleon_Scatter::scatmax(double DME){
	return std::min(DME,DME-Escatmin);
}

double Nucleon_Scatter::scatmin(double DME, double DMM, double NMM){
	return std::max(Efmin(DME, DMM, NMM), DME-Escatmax);
}
