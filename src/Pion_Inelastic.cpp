#include "Scatter.h"

#include "minimization.h"
#include "DMNscattering.h"
#include "Kinematics.h"
#include "Random.h"
#include "constants.h"
#include "decay.h"

#include <cmath>
#include <iostream>
#include <algorithm>
#include <fstream>

const double mp = MASS_PROTON;
const double mn = MASS_NEUTRON;
const double Mdelta = MASS_DELTA;
const double convmcm = 100.0;
const double tol_abs=1e-20;
const double tol_frac=1e-10;


using std::cout;
using std::endl;
using namespace std::placeholders;
using std::string; using std::vector;
using std::shared_ptr;
using std::list;

//This may eventually get moved into Parameter.
const string form_factor_filename = "data/delta_production_form_factor.dat";

double Pion_Inelastic::Ermax(double E, double mx, double MN){
    return ((E + MN)*sqrt(pow(pow(Mdelta,2) + 2*E*MN + pow(MN,2),2)/(2*E*MN + pow(MN,2) + pow(mx,2)))*sqrt(pow(2*E*MN + pow(MN,2) + pow(mx,2),2)) + 
     sqrt((E - mx)*(E + mx))*sqrt((2*E*MN + pow(MN,2) + pow(mx,2))*(2*E*MN + pow(MN,2) + Mdelta*(-Mdelta + 2*mx))*(2*E*MN + pow(MN,2) - Mdelta*(Mdelta + 2*mx))))/
   (2.*(E + MN)*sqrt((2*E*MN + pow(MN,2) + pow(mx,2))/pow(E + MN,2))*sqrt(pow(2*E*MN + pow(MN,2) + pow(mx,2),2)));
}

double Pion_Inelastic::Ermin(double E, double mx, double MN){
    return ((E + MN)*sqrt(pow(pow(Mdelta,2) + 2*E*MN + pow(MN,2),2)/(2*E*MN + pow(MN,2) + pow(mx,2)))*sqrt(pow(2*E*MN + pow(MN,2) + pow(mx,2),2)) - 
     sqrt((E - mx)*(E + mx))*sqrt((2*E*MN + pow(MN,2) + pow(mx,2))*(2*E*MN + pow(MN,2) + Mdelta*(-Mdelta + 2*mx))*(2*E*MN + pow(MN,2) - Mdelta*(Mdelta + 2*mx))))/
   (2.*(E + MN)*sqrt((2*E*MN + pow(MN,2) + pow(mx,2))/pow(E + MN,2))*sqrt(pow(2*E*MN + pow(MN,2) + pow(mx,2),2)));
}

double Pion_Inelastic::Er_to_theta(double Ex, double EDelta, double mx, double MN){
    //cout << "Ex=" << Ex << " EDelta=" << EDelta << " MN=" << MN << " Bot=" << (2*sqrt(EDelta*EDelta-Mdelta*Mdelta)*sqrt(Ex*Ex-mx*mx)) << " theta=" << acos((2*EDelta*MN+2*EDelta*Ex-2*Ex*mn-Mdelta*Mdelta-mn*mn)/(2*sqrt(EDelta*EDelta-Mdelta*Mdelta)*sqrt(Ex*Ex-mx*mx))) << endl; 
    return acos((2*EDelta*MN+2*EDelta*Ex-2*Ex*mn-Mdelta*Mdelta-mn*mn)/(2*sqrt(EDelta*EDelta-Mdelta*Mdelta)*sqrt(Ex*Ex-mx*mx)));
}
/*
En: Incident DM energy
ER: Recoiling nucleon energy.
mN: Mass of nucleon. Looks like same form factor for both? Weird.
*/
double Pion_Inelastic::dsigma_dER_N(double En, double ER, double mx, double mA, double alphaprime, double kappa, double mN){
return -(alphaEM*alphaprime*pow(Mdelta + mN,2)*
      (2*pow(Mdelta,2)*(pow(Mdelta,2) - 2*ER*mN + pow(mN,2))*
         (pow(Mdelta,2) - 2*ER*mN + pow(mN,2) - 4*pow(mx,2)) + 
        pow(Mdelta,4)*(-pow(Mdelta,2) + 2*ER*mN - pow(mN,2) + 4*pow(mx,2)) + 
        pow(mN,4)*(-pow(Mdelta,2) + 2*ER*mN - pow(mN,2) + 4*pow(mx,2)) - 
        2*pow(mN,2)*(2*pow(Mdelta,2) - 2*ER*mN + pow(mN,2))*
         (-pow(Mdelta,2) + 2*ER*mN - pow(mN,2) + 4*pow(mx,2)) + 
        (pow(Mdelta,2) - 2*ER*mN + pow(mN,2))*
         ((pow(Mdelta,2) - 2*ER*mN + pow(mN,2))*
            (-pow(Mdelta,2) + 2*ER*mN - pow(mN,2) + 4*pow(mx,2)) + 
           pow(-2*ER*mN - 2*pow(mx,2) + 2*(2*En*mN + pow(mN,2) + pow(mx,2)),2)))
       *pi*pow(kappa,2)*pow(GM(-pow(Mdelta,2) + 2*ER*mN - pow(mN,2)),2))/
   (2.*mN*pow(pow(mA,2) - pow(Mdelta,2) + 2*ER*mN - pow(mN,2),2)*
     (-pow(Mdelta,2) + 2*ER*mN - pow(mN,2) + pow(Mdelta + mN,2))*
     TriangleFunc2(sqrt(2*En*mN + pow(mN,2) + pow(mx,2)),mx,mN))*final_branch;
}

double Pion_Inelastic::GM(double q2){
    if(q2>MAX_Q2){
        if(!MAX_Q2_WARNING){
            MAX_Q2_WARNING=true;
            std::cerr << "First request of out of bounds q^2 for pion inelastic form factor.\n";
            std::cerr << "q^2 tried =" << q2 << "\nMAX_Q2 = " <<  MAX_Q2 << endl;
        }
        return 0.0;
    }
    return form_factor->Interpolate(q2);
}

//Minimum total energy for a DM particle to produce a Delta
double Edm_kinetic_min(double mx, double mN){
    return (pow(Mdelta,2)+2*Mdelta*mx-mN*mN)/(2*mN);
}

/*
NEmin, NEmax are cuts on the energy of the recoiling nucleon.
Emini, Emaxi, Eresi are the min and max energy of the incoming dm 
particle, and the energy resolution that should be used to build
the scattering interpolation function.
final_state=0 is default and implies a pi0 final state
final_state=1 means a photon final state.
*/
Pion_Inelastic::Pion_Inelastic(double Emini, double Emaxi, double Eresi, double MDM, double MV, double alphaprime, double kappa, double NEmax, double NEmin, int final_state){
    if(final_state==1){
        final_mass=MASS_PHOTON;
        final_branch=Delta_to_gamma;
        final_name="photon";
    }
    //Inelegant, but straightforward.
    //This allows for the possibility of charged pion decays.
    if(final_state==2){
        final_branch=Delta_to_pi0+Delta_to_pion_charged;
        final_mass=0;
        final_name="pion";//The code will look for this final_name to trigger charged pion handling
    }
    else{
        final_mass=mpi0;
        final_branch=Delta_to_pi0;
        final_name="pi0";
    }
    MAX_Q2_WARNING = false;
    Edmmin=std::max(Emini,std::max(Edm_kinetic_min(MDM, mn),Edm_kinetic_min(MDM, mp))); Edmmax=Emaxi; Edmres=Eresi;
    Escatmin=NEmin;
	Escatmax=NEmax;
	//std::cout << Edmmax << " " << Edmmin << endl;
    form_factor = shared_ptr<Linear_Interpolation>();
    load_form_factor(form_factor_filename, form_factor);
    set_Model_Parameters(MDM, MV, alphaprime, kappa);
    
/*    for(double i = Ermin(1,0.01,mp); i < Ermax(1,0.01,mp); i+=0.001){
        cout << i; 
        double mxvals[]={0.01,0.1,0.2};
        for(int j = 0; j<3; j++)
            cout << " " << convGeV2cm2*dsigma_dER_N(1, i, mxvals[j], 0.5, 0.1, 1e-3, mp);
        cout << endl;
    }*/
}

//No error checking for q2_max! 
void Pion_Inelastic::load_form_factor(const string& filename, shared_ptr<Linear_Interpolation>& ff){
    std::ifstream instream(filename);
    if(!instream.is_open()){
        std::cerr << "Pion_Inelastic cannot open " << filename << " needed for magnetic dipole moment form factor." << std::endl;
        throw -1;
    }
    double q2_min,q2_max;
    std::vector<double> q2_dist;
    double in;
    instream >> q2_min;
    instream >> in;
    q2_dist.push_back(in);
    while(instream >> q2_max){
        instream >> in;
        q2_dist.push_back(in);
    }

    MAX_Q2 = q2_max;

    ff = shared_ptr<Linear_Interpolation>( new Linear_Interpolation(q2_dist, q2_min, q2_max));
}

void Pion_Inelastic::prep_ab(double &a, double &b, double Edm, double mx, double mN){
    a=Ermin(Edm,mx,mN);
    b=(Ermax(Edm,mx,mN)+Ermin(Edm,mx,mN))/2;
}

//Note that I should add the branching ratio of Delta->N\pi, but it seems to be 1.
//This is odd, a Delta should be able to go to a charged pion, right?
//Turns out it's 2/3!
void Pion_Inelastic::generate_cross_sections(){
	std::vector<double> vec_proton_maxima;
	std::vector<double> vec_neutron_maxima;
	std::vector<double> vec_proton;
	std::vector<double> vec_neutron;
//	double emax; double emin;
	
	double a,b,c,xmin;
    for(double iter=Edmmin+Edmres; iter<=Edmmax; iter+=Edmres){
        //cout << "Edm=" << iter << " Ermin=" <<  Ermin(iter,mdm,mp) << " Ermax=" << Ermax(iter,mdm,mp) << endl;
        if(Ermax(iter,mdm,mp)<=Ermin(iter,mdm,mp)){
            vec_proton.push_back(0);
            vec_proton_maxima.push_back(0);
        }
        else{
            std::function<double(double)> fp = std::bind(&Pion_Inelastic::dsigma_dER_N,this,iter,_1,mdm,MDP,alD,kap,mp);
            //lim_func_wrapper comes from minimization.h. It evaluates to zero outside
            //specified boundaries, and flips the sign of its function for use by
            //minimization.
			std::function<double(double)> fplim = std::bind(lim_func_wrapper,_1,0.0,fp,Ermin(iter,mdm,mp),Ermax(iter,mdm,mp));
    		vec_proton.push_back(DoubleExponential_adapt(fp,Ermin(iter,mdm,mp),Ermax(iter,mdm,mp),100,0.1,1e-4));
            prep_ab(a,b,iter,mdm,mp);
    		//cout << "a=" << a << "b=" << b << "Ermax=" << Ermax(iter,mdm,mp) << endl; 
            c=mnbrak(a,b,fplim);
    		xmin=0;
            vec_proton_maxima.push_back(-1.0*golden(a,b,c,fplim,tol_frac,tol_abs,xmin));
            //cout << "xsec=" << vec_proton.back() << " maxima=" << vec_proton_maxima.back() << endl;
        }
		
        if(Ermax(iter,mdm,mn)<=Ermin(iter,mdm,mn)){
            vec_neutron.push_back(0);
            vec_neutron_maxima.push_back(0);
        }
        else{
    		std::function<double(double)> fn = std::bind(&Pion_Inelastic::dsigma_dER_N,this,iter,_1,mdm,MDP,alD,kap,mn);
            std::function<double(double)> fnlim = std::bind(lim_func_wrapper,_1,0.0,fn,Ermin(iter,mdm,mn),Ermax(iter,mdm,mn));   
            vec_neutron.push_back(DoubleExponential_adapt(fn,Ermin(iter,mdm,mn),Ermax(iter,mdm,mn),100,0.1,1e-4));
            prep_ab(a,b,iter,mdm,mn);
            c=mnbrak(a,b,fnlim);
    		xmin=0;
    		vec_neutron_maxima.push_back(-1.0*golden(a,b,c,fnlim,tol_frac,tol_abs,xmin));
        }
    }

    //cout << "End of loop\n";
	proton_cross = std::unique_ptr<Linear_Interpolation>(new Linear_Interpolation(vec_proton,Edmmin,Edmmax));
    neutron_cross = std::unique_ptr<Linear_Interpolation>(new Linear_Interpolation(vec_neutron,Edmmin,Edmmax));
	proton_cross_maxima = std::unique_ptr<Linear_Interpolation>(new Linear_Interpolation(vec_proton_maxima,Edmmin,Edmmax));
	neutron_cross_maxima = std::unique_ptr<Linear_Interpolation>(new Linear_Interpolation(vec_neutron_maxima,Edmmin,Edmmax));


   //for(double i=Edmmin; i<=Edmmax; i+=0.01)
   //     cout << i << " " << proton_cross->Interpolate(i) << " " << neutron_cross->Interpolate(i) << endl;
}


/*
This function takes a list of particles involved in the scattering and the dark matter from the main code, runs 
probscatter with the final state scattered particle (pion) and the dark matter. It then checks that the pion satisfies
any angular cuts, and inserts it into the particle list for later recording.

This might not be a pion anymore if final_state==1, but it probably will be.
*/
bool Pion_Inelastic::probscatter(std::shared_ptr<detector>& det, list<Particle>& partlist, list<Particle>::iterator& DMit){
    Particle pion(final_mass);
    pion.name = final_name;
    Particle Delta(Mdelta);
    Delta.name = "Delta";
    Particle Nucleon(0);
    if(probscatter(det, *DMit, pion, Delta, Nucleon)&&(min_angle<=0||pion.Theta()>min_angle)&&(max_angle>2*pi||pion.Theta()<max_angle)&&(pion.Kinetic_Energy()>Escatmin)&&(pion.Kinetic_Energy()<Escatmax)){
        DMit->Generate_Position();
        Link_Particles(*DMit, Delta);
        Link_Particles_Immediate(Delta, pion);
        Link_Particles_Immediate(Delta, Nucleon);
        Particle DMout(DMit->m);
        DMout.name = "Recoil_DM";
        DMout.ThreeMomentum(DMit->px-Delta.px,DMit->py-Delta.py,DMit->pz-Delta.pz);
        Link_Particles(*DMit, DMout);
        //Insert in reverse display order.
        partlist.insert(std::next(DMit),pion);
        partlist.insert(std::next(DMit),Nucleon);
        partlist.insert(std::next(DMit),Delta);
        partlist.insert(std::next(DMit),DMout);
        return true;
    }
    return false;
}

/*
Attempts to scatter a DM particle off detector nucleons, without simulating end state particles.
This is mainly used for burn-in.
I need to add testing for decay to pi^0.
*/
bool Pion_Inelastic::probscatter(std::shared_ptr<detector>& det, Particle &DM){ 
    if(DM.E<Edmmin)
        return false;   
    double LXDet = convmcm*(det->Ldet(DM));
    double XDMp = proton_cross->Interpolate(DM.E)*(det->PNtot());
    double XDMn = neutron_cross->Interpolate(DM.E)*(det->NNtot());
    double prob=LXDet*convGeV2cm2*(XDMp+XDMn);
	//std::cout << DM.E << " " << XDMp << " " << XDMn << " " << prob << " " << endl;
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

/*
Scatters a DM particle off a detector nucleon. Stores the outgoing pion in &pion.
*/

bool Pion_Inelastic::probscatter(std::shared_ptr<detector>& det, Particle &DM, Particle &pion, Particle &Delta, Particle &Nucleon){ 
	if(DM.E<Edmmin)
        return false;
    double LXDet = convmcm*(det->Ldet(DM));
    double XDMp = proton_cross->Interpolate(DM.E)*(det->PNtot());
    double XDMn = neutron_cross->Interpolate(DM.E)*(det->NNtot());
    double prob=LXDet*convGeV2cm2*(XDMp+XDMn);

    if(prob > pMax*Random::Flat(0,1)){
        if(prob > pMax){
        	pMax = prob;
        }
		if(XDMp/(XDMn+XDMp) >= Random::Flat(0,1)){
			std::function<double(double)> Xsec = std::bind(&Pion_Inelastic::dsigma_dER_N,this,DM.E,_1,DM.m,MDP,alD,kap,mp);
			//Not sure if this check is still necessary.
            if(Ermin(DM.E,DM.m,mp) >= Ermax(DM.E, DM.m, mp)){
				return false;
            }
            
            scatterevent(DM, Delta, Xsec, *proton_cross_maxima,mp);
            //Temporary, may not be a proton!
            Particle Proton(mp);
            Proton.name="Proton";
            //This could be generalized, not sure the best way to do it yet.
            if(final_name=="pion"){
                if(Random::Flat()>1.0/3.0){
                    pion.name="pi0";
                    pion.m=mpi0;
                }
                else{
                    pion.name="pi+";
                    pion.m=MASS_PION_CHARGED;
                    Proton.name="Neutron";
                    Proton.m = mn;
                }
            }
            TwoBodyDecay(Delta, Proton, pion);
            Nucleon = Proton;
        }
        else{
			std::function<double(double)> Xsec = std::bind(&Pion_Inelastic::dsigma_dER_N,this,DM.E,_1,DM.m,MDP,alD,kap,mn);
            if(Ermin(DM.E,DM.m,mn) >= Ermax(DM.E, DM.m, mn))
				return false;
			scatterevent(DM, Delta, Xsec, *neutron_cross_maxima,mn);
            Particle Neutron(mn);
            Neutron.name = "Neutron";
            //This could be generalized, not sure the best way to do it yet.
            if(final_name=="pion"){
                if(Random::Flat()>1.0/3.0){
                    pion.name="pi0";
                    pion.m=mpi0;
                }
                else{
                    pion.name="pi-";
                    pion.m=MASS_PION_CHARGED;
                    Neutron.name="Proton";
                    Neutron.m = mp;
                }
            }
            TwoBodyDecay(Delta, Neutron, pion);
            Nucleon = Neutron;
        }
        return true;
    }
    else
        return false;
}

void Pion_Inelastic::scatterevent (Particle &DM, Particle &Delta, std::function<double(double)> Xsec, Linear_Interpolation& Xmax, double mN){
    double Delta_E_max = Ermax(DM.E, DM.m, mN);
    double Delta_E_min = Ermin(DM.E, DM.m, mN);
	double dsigmax = std::max(Xsec(Delta_E_min),Xmax.Interpolate(DM.E));
    double xe,thetaN,phiN,pN;
    while(true){
        xe = Random::Flat(Delta_E_min,Delta_E_max);
        //cout << "xe=" << xe << " Ermin" << Ermin(DM.E, DM.m, mN) << " Ermax " << Ermax(DM.E, DM.m, mN) << endl;
        if(Xsec(xe)/dsigmax > Random::Flat(0,1)){
            thetaN = Er_to_theta(DM.E,xe,DM.m,mN);
            //cout << "theta=" << thetaN << endl;
            phiN = Random::Flat(0,1)*2*pi;
            //pN = sqrt(pow(DM.E+mN-xe,2)-pow(Delta.m,2));
            pN = sqrt(xe*xe-pow(Delta.m,2));
            //cout << "pN=" << pN << endl;
            Delta.ThreeMomentum(pN*sin(thetaN)*cos(phiN),pN*sin(thetaN)*sin(phiN),cos(thetaN)*pN);
            Delta.Rotate_y(DM.Theta());
            Delta.Rotate_z(DM.Phi());
            break;
        }
    }
}
