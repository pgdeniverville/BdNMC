#include "Drell_Yan_Gen.h"
#include "Kinematics.h"
#include "constants.h"

#include <ios>
#include <iostream>
#include <exception>
#include <stdexcept>
#include <vector>

using std::string; using std::cerr;
using std::endl; using std::cout;
using std::vector
;
using namespace annihilation_to_pair;

const string default_file_path = "data/ct10n.00.pds";
//const string default_file_path = "data/CT14LN.pds";

//Only going over u/u-bar and d/d-bar
const double p_qf[] = {1.0/3.0,-2.0/3.0,2.0/3.0,-1.0/3.0};
const double n_qf[] = {-2.0/3.0,1.0/3.0,-1.0/3.0,2.0/3.0};
const int indices[] = {-2,-1,1,2};
const int qbar_indices[] = {2,1,-1,-2};

//const double QMIN = 1.3;
//const double QMIN = 1.3;

const double _XMIN_=1e-8;
const double _XMAX_=1.0;

double Drell_Yan_Gen::sig_test(double x, double y, int mode){
	double Q;
	double tot=0;
	if(mode==0){
		Q=sqrt(Q2_eval(BeamEnergy,MASS_PROTON,x,y));
		if(Q<QMIN){
			return 0;
		}
	//	cout << Q << " " << x << " " << y << endl;
		for(int i = 0; i<4; i++){
			tot+=proton_pdf.parton(indices[i],x,Q)*proton_pdf.parton(qbar_indices[i],y,Q)*sig_tot(x,y,MASS_PROTON);
		}
	}
	else{
		Q=sqrt(Q2_eval(BeamEnergy,MASS_NEUTRON,x,y));
		if(Q<QMIN){
			return 0;
		}
		for(int i = 0; i<4; i++){
			tot+=neutron_pdf.parton(indices[i],x,Q)*neutron_pdf.parton(qbar_indices[i],y,Q)*sig_tot(x,y,MASS_NEUTRON);
		}
	}
	return tot;
}

double Drell_Yan_Gen::Monte_Carlo_test(int mode){
	double tot1=0;
	double tot2=0;
	int n=10000;
	double x,y;
	do{
		tot1=(tot1+2*tot2)/3.0;
		//tot1=tot2;
		tot2=0;
		for(int i=0; i<n; i++){
			x=Random::Flat();
			y=Random::Flat();
			tot2+=sig_test(x,y,mode);
		}
		tot2/=n;
//		cout << "n=" << n << " tot1=" << tot1 << " tot2=" << tot2 << endl;
		n*=2;		
	} while(2*std::abs(tot1-tot2)/(tot1+tot2)>TOLERANCE);

	return (tot1+2*tot2)/3.0;
}

double Drell_Yan_Gen::Monte_Carlo_test_2(int mode){
	double tot1=0;
	double tot2=0;
	int n=10000;
	double x,y,Q,MASS;
	if(mode==0){
		MASS=MASS_PROTON;
	}
	else{
		MASS=MASS_NEUTRON;
	}
	do{
		tot1=(tot1+2*tot2)/3.0;
		//tot1=tot2;
		tot2=0;
		for(int i=0; i<n; i++){
			x=Random::Flat();
			y=Random::Flat();
			Q=sqrt(Q2_eval(BeamEnergy,MASS,x,y));
			if(Q<QMIN)
				continue;
			tot2+=(E3Max(BeamEnergy*x,MASS_PROTON*x,y*MASS,mx)-E3Min(BeamEnergy*x,MASS_PROTON*x,y*MASS,mx))*dsig(0,x,y,Q,mode);
		}
		tot2/=n;
//		cout << "n=" << n << " tot1=" << tot1 << " tot2=" << tot2 << endl;
		n*=2;		
	} while(2*std::abs(tot1-tot2)/(tot1+tot2)>TOLERANCE);
	return (tot1+2*tot2)/3.0;
}

double Drell_Yan_Gen::dsigma_dt(double x, double y, double mass_target, double s_hat){
	double ma = 1;
	double MASS=mass_target;
	double EA=BeamEnergy;
	double gq=1;
	double gchi=1;
	return pow(gq*gchi,2)*s_hat*(s_hat-pow(MASS_PROTON*x-MASS*y,2))/192.0/pi/pow(MASS,2)/pow(x*y,2)/(pow(EA,2)-pow(MASS_PROTON,2))/pow(ma*ma-s_hat,2);;
}

double Drell_Yan_Gen::sig_tot(double x, double y, double mass_target){
	double s_hat = shat(BeamEnergy*x,MASS_PROTON*x,mass_target*y);
	//double E1cm_hat = (s_hat+pow(x*MASS_PROTON,2)-pow(y*mass_target,2))/(2*sqrt(s_hat));
	//double p1cm_hat = sqrt(pow(E1cm_hat,2)-pow(x*MASS_PROTON,2));
	//double E3cm = sqrt(s_hat)/2;
	//double p3cm = sqrt(pow(E3cm,2)-pow(mx,2));
	double _t0_=t0(BeamEnergy*x,MASS_PROTON*x,mass_target*y,mx);
	double _t1_=t1(BeamEnergy*x,MASS_PROTON*x,mass_target*y,mx);
	if(_t0_<_t1_){
		return 0;
	}
	return (_t0_-_t1_)*dsigma_dt(x,y,mass_target,s_hat);
}

void Drell_Yan_Gen::Test(){
	cout << "Beginning Test" << endl;

	double tmp1=Monte_Carlo_test(0);//*convGeV2cm2/microbarn;
	double tmp2=Monte_Carlo_test(1);//*convGeV2cm2/microbarn;

	cout << "Monte_Carlo_test 1\n"; 

	cout << tmp1 << " " << tmp2 << endl;

	tmp1=Monte_Carlo_test_2(0);//*convGeV2cm2/microbarn;
	tmp2=Monte_Carlo_test_2(1);//*convGeV2cm2/microbarn;

	cout << "Monte_Carlo_test 2\n";

	cout << tmp1 << " " << tmp2 << endl;
/*
	cout << "Theta test" << endl;

	double E1,m1,m2,m3;
	E1 =8; m1=1; m2=0.9; m3=0.2;
	double t_low = t1(E1,m1,m2,m3);
	double t_high = t0(E1,m1,m2,m3);
	double E_low = E3Min(E1,m1,m2,m3);
	double E_high = E3Max(E1,m1,m2,m3);
	cout << "Elow=" << E_low << " theta= " << Theta_from_E3(E1,m1,m2,E_low*1.01,m3) << endl;
	cout << "E_highgh=" << E_high << " theta= " << Theta_from_E3(E1,m1,m2,E_high,m3) << endl;*/
}

double Drell_Yan_Gen::Q2_eval(double E1, double m_target, double x, double y){
	return pow(x*MASS_PROTON,2)+pow(y*m_target,2)+2*E1*m_target*x*y;
}

double Drell_Yan_Gen::Y_Min_from_QMIN(double E1, double m1, double m2, double x){
	return (sqrt(m2*(x*x*(E1-m1)*(E1+m1)+pow(QMIN,2)))-E1*x)/m2;
}

double Drell_Yan_Gen::Interaction_Cross_Section(){
	return (sig_p*target_n+sig_n*target_n);
}

//const string pdf_file_path, const string pdf_neutron_file_path, const string channel_name, double TOL, double energy_bins)
Drell_Yan_Gen::Drell_Yan_Gen(double mchi, double Beam_Energy, std::function<double(double, double, double, double, double)> dsig, double tar_p, double tar_n, production_channel& prodchan){
	if(2*BeamEnergy*MASS_PROTON+2*pow(MASS_PROTON,2)<4*mchi*mchi){
		cerr << "Dark matter mass too large to be produced by this beam through Drell_Yan\n" << endl;
		throw std::domain_error("Energy conservation violation.");
	}
	QMIN = prodchan.QMin();
	OFF_SHELL=true;
	target_n = tar_n;
	target_p = tar_p;
	mx = mchi;
	chan_name = prodchan.Production_Channel();

	
	BeamEnergy=Beam_Energy;
	TOLERANCE = prodchan.TOL();
	string proton_pdf_file_path = prodchan.Proton_PDF_File();
	string neutron_pdf_file_path = prodchan.Neutron_PDF_File();//This is not used at the moment.

	dsig_hat = dsig;

	if(proton_pdf_file_path!=""){
		//proton_pdf = cteqpdf();
		std::ifstream infile;
     	infile.open(proton_pdf_file_path);
		if (!infile) {
	     	cerr << "error: unable to open input file: " << proton_pdf_file_path << endl;
	     	throw -1;
     	}
		proton_pdf.setct11(proton_pdf_file_path);
		neutron_pdf.setct11(proton_pdf_file_path);
	}
	else{
		std::ifstream infile;
		infile.open(default_file_path);
		if(!infile){
			cerr << "error: unable to open input file: " << default_file_path << endl;
	     	throw -1;
		}
		cout << "Loading CT10NLO PDF.\nSee http://hep.pa.msu.edu/cteq/public/ct10.html for more information." << endl;
		proton_pdf.setct11(default_file_path);
		neutron_pdf.setct11(default_file_path);
	}

	//cout << proton_pdf.parton(2,0.3,2) << " " << proton_pdf.parton(-2,0.15,2) << endl;

	//This is not exactly right. Neglecting difference in proton and neutron mass.
	Efmax = E3Max(BeamEnergy, MASS_PROTON, MASS_PROTON, mx);
	//x,y->0. QMin actually places a limit on this.
	Efmin = mx;
	Efres=(double)(Efmax-Efmin)/prodchan.Energy_Bins();
	//This is going to evaluate the production cross section.
	Evaluate_Cross_Sections();
	//Test();
}

/*
void Drell_Yan_Gen::Evaluate_Cross_Sections(){
	pMax_n = 0;
	pMax_p = 0;
	sig_p=0;
	sig_n=0;

	sig_p = Monte_Carlo(0);
	sig_n = Monte_Carlo(1);
}
*/
double Drell_Yan_Gen::Monte_Carlo(int mode){
	double tot1=0;
	double tot2=0;
	int n=10000;
	double x,y;
	do{
		tot1=(tot1+2*tot2)/3.0;
		//tot1=tot2;
		tot2=0;
		for(int i=0; i<n; i++){
			x=Random::Flat(_XMIN_,_XMAX_);
			y=Random::Flat();
			tot2+=sig(x,y,mode);
		}
		tot2/=n;
		n*=2;		
	} while(2*std::abs(tot1-tot2)/(tot1+tot2)>TOLERANCE);

	return (tot1+2*tot2)/3.0;
}


void Drell_Yan_Gen::Evaluate_Cross_Sections(){
	vector<double> dsig_p;
	vector<double> dsig_n;
	vector<double> dsig_p_max;
	vector<double> dsig_n_max;
	sig_p=0;
	sig_n=0;
	double max;
	vector<double> Ekvals;
	double Ek;
	pMax_p = 0;
	pMax_n = 0;
	//Test();
	cout << "Generating Drell Yan cross sections with QMIN = " << QMIN << " GeV." << endl;
	for(Ek=mx; Ek<Efmax; Ek+=Efres){
		//cout << Ek << "\n";
		//cout << "Proton\n";
		max=0;
		dsig_p.push_back(Monte_Carlo(Ek,0,max));
		dsig_p_max.push_back(max);
		sig_p+=dsig_p.back();
		//cout << Ek << " " << dsig_p.back() << endl;
		if(max>pMax_p){
			pMax_p=max;
		}
		//cout << "Neutron\n";
		max=0;
		dsig_n.push_back(Monte_Carlo(Ek,1,max));
		dsig_n_max.push_back(max);
		sig_n+=dsig_n.back();
		if(max>pMax_n){
			pMax_n=max;
		}
		Ekvals.push_back(Ek);

	}
	sig_p*=Efres;
	sig_n*=Efres;

	//cout << sig_p << endl;
	//throw -1;
	cout << "sigma_proton = " << sig_p*convGeV2cm2/microbarn << " microbarn sigma_neutron = " << sig_n*convGeV2cm2/microbarn << " microbarn" << endl;
	dsigma_dE_proton = Linear_Interpolation2(Ekvals, dsig_p);
	dsigma_dE_neutron = Linear_Interpolation2(Ekvals, dsig_n);
	dsigma_dE_proton_max = Linear_Interpolation(dsig_p_max, Efmin, Ek);
	dsigma_dE_neutron_max = Linear_Interpolation(dsig_n_max, Efmin, Ek);

	//dsigma_dE_proton.Report();
	dsigma_dE_proton.Convert_to_CDF();
	//cout << "CDF" << endl;
	//dsigma_dE_proton.Report();
	dsigma_dE_proton.Invert();
	//cout << "Inverted" << endl;
	//dsigma_dE_proton.Report();
	//cout << "Neutron" << endl;
	//dsigma_dE_neutron.Report();
	dsigma_dE_neutron.Convert_to_CDF();
	//cout << "CDF" << endl;
	//dsigma_dE_neutron.Report();
	dsigma_dE_neutron.Invert();
	//cout << "Inverted" << endl;
	//dsigma_dE_neutron.Report();
}

double Drell_Yan_Gen::dsig(double Ek, double x, double y, double Q, int mode){
	double tot=0;	
	if(mode==0){
		for(int i = 0; i<4; i++){
			tot+=proton_pdf.parton(indices[i],x,Q)*proton_pdf.parton(qbar_indices[i],y,Q)*dsig_hat(Ek,x,y,p_qf[i],MASS_PROTON);
		}
	}
	else{
		for(int i = 0; i<4; i++){
			tot+=neutron_pdf.parton(indices[i],x,Q)*neutron_pdf.parton(qbar_indices[i],y,Q)*dsig_hat(Ek,x,y,n_qf[i],MASS_NEUTRON);
		}
	}
	return tot;
}

double Drell_Yan_Gen::sig(double x, double y, int mode){
	double Q;
	double tot=0;
	if(mode==0){
		
		Q=sqrt(Q2_eval(BeamEnergy,MASS_PROTON,x,y));
		if(Q<QMIN){
			return 0;
		}

		for(int i = 0; i<4; i++){
			tot+=proton_pdf.parton(indices[i],x,Q)*proton_pdf.parton(qbar_indices[i],y,Q)*sig_hat(x,y,p_qf[i],MASS_PROTON);
		}
	}
	else{
		
		Q=sqrt(Q2_eval(BeamEnergy,MASS_NEUTRON,x,y));
		if(Q<QMIN){
			return 0;
		}

		for(int i = 0; i<4; i++){
			tot+=neutron_pdf.parton(indices[i],x,Q)*neutron_pdf.parton(qbar_indices[i],y,Q)*sig_hat(x,y,n_qf[i],MASS_NEUTRON);
		}
	}
	return tot;
}

//This needs testing
bool Drell_Yan_Gen::Verify(double Ek, double mass_target, double x, double y, double &Q){
	if(E3Min(x*BeamEnergy,x*MASS_PROTON,y*mass_target,mx)>Ek||E3Max(x*BeamEnergy,x*MASS_PROTON,y*mass_target,mx)<Ek){
		return false;
	}
	Q=sqrt(Q2_eval(BeamEnergy, mass_target, x, y));
	if(Q<QMIN){
		return false;
	}
	return true;
}

double Drell_Yan_Gen::Monte_Carlo(double Ek, int mode, double& max){
	double tot1=0;
	double tot2=0;
	int n=10000;
	double x,y,hold;
	double mass_target,Q;
	if(mode==0)
		mass_target = MASS_PROTON;
	else
		mass_target = MASS_NEUTRON;
	//cout << "Ek=" << Ek << endl;
	do{
		tot1=(tot1+2*tot2)/3.0;
		//tot1=tot2;
		tot2=0;
		for(int i=0; i<n; i++){
			x=Random::Flat(_XMIN_,_XMAX_);
			y=Random::Flat();
			hold=0;
			if(!Verify(Ek,mass_target,x,y,Q)){
				continue;
			}
			if((hold=dsig(Ek,x,y,Q,mode))>max){
				max=hold;
			}
			//DEBUG
			
			if(hold<0){
				cout << "Drell_Yan producing negative cross section\n" << endl;
				cout << hold << " " << Ek << " " << x << " " << y << " " << Q << " " << mode << endl;
				for(int index = -2; index<3; index++){
					cout << index << " " << proton_pdf.parton(index,x,Q) << " " << proton_pdf.parton(index,y,Q) << endl;
				}
				continue;
			}
			tot2+=hold;
		}
		tot2/=n;
		n*=2;
		//cout << "Monte Carlo iteration " << tot1 << " " << tot2 << " " << n << " " << std::abs(tot1-tot2) << " " << 2*std::abs(tot1-tot2)/(tot1+tot2) << endl;
	} while(2*std::abs(tot1-tot2)/(tot1+tot2)>TOLERANCE);
	//cout << Ek << " " << (tot1+2*tot2)/3.0 << endl;
	return (tot1+2*tot2)/3.0;
}
/*
bool Drell_Yan_Gen::GenDM(std::list<Particle>& vec, std::function<double(Particle&)> det_int, Particle& part){
	Particle chi1(mx); 
	Particle chi2(mx);

	chi1.name = product_name;
	chi2.name = product_name;

	double p_p = sig_p*target_p;
	double p_n = sig_n*target_n;

	double mass_target,mode,x,y,hold,Q,t;
	double *pMax;

	if(Random::Flat()*(p_p+p_n)<p_p){
		mass_target = MASS_PROTON;
		mode = 0;
		pMax = &pMax_p;
	}
	else{
		mass_target = MASS_NEUTRON;
		mode = 1;
		pMax = &pMax_n;
	}

	while(true){
		x=Random::Flat();
		y=Random::Flat();

		if((Q=sqrt(Q2_eval(BeamEnergy,mass_target,x,y)))<QMIN){
			continue;
		}
		t=Random::Flat(t1(BeamEnergy*x,MASS_PROTON*x,mass_target*y,mx),t0(BeamEnergy*x,MASS_PROTON*x,mass_target*y,mx));
		if((hold=dsig(t,x,y,Q,mode))>Random::Flat()*(*pMax)){
			if(*pMax<hold){
				*pMax = hold;
			}
			break;
		}
	}

	double Ek = E3_from_t(x*BeamEnergy, MASS_PROTON*x, mass_target*y, mx, t);

	double theta = Theta_from_E3(x*BeamEnergy, x*MASS_PROTON, y*mass_target, Ek, mx);

	chi1.ThreeMomentumPolar(sqrt(Ek*Ek-mx*mx),theta,Random::Flat(0,2*pi));
	if(part.px!=0||part.py!=0){
		chi1.Rotate_y(part.Theta());
		chi1.Rotate_z(part.Phi());
	}
	chi2.ThreeMomentum(part.px*x-chi1.px,part.py*x-chi1.py,part.pz*x-chi1.pz);
	//cout << "theta " << theta << " x " << x << " y " << y << " Ek " << Ek << " mass_target " << mass_target <<  endl;
	//chi1.report(cout);
	//chi2.report(cout);

	vec.push_back(part);
 	double intersect1=det_int(chi1);
    double intersect2=det_int(chi2);

    //chi1.report(cout);
    //chi2.report(cout);

    if((intersect1)>0 || (intersect2)>0){;
        if(intersect1>0){
            vec.push_back(chi1);
        }
        if(intersect2>0){
            vec.push_back(chi2);
        }
        return true;
    }  
    else
        return false;

}*/

bool Drell_Yan_Gen::GenDM(std::list<Particle>& vec, std::function<double(Particle&)> det_int, Particle& part){
	Particle chi1(mx); 
	Particle chi2(mx);

	part.name = "Proton";

	//cout << "GENDM" << endl;

	chi1.name = product_name;
	chi2.name = product_name;

	double p_p = sig_p*target_p;
	double p_n = sig_n*target_n;
	//This is sloppy. Better way to do this?
	Linear_Interpolation2& dsig_i = dsigma_dE_proton;
	Linear_Interpolation& dsig_max_i = dsigma_dE_proton_max;
	double mass_target,mode;
	//cout << "p_p " << p_p << endl;
	//cout << "p_n " << p_n << endl;
	if(Random::Flat()*(p_p+p_n)<p_p){
		mass_target = MASS_PROTON;
		mode = 0;
	}	
	else{
		dsig_i = dsigma_dE_neutron;
		dsig_max_i=dsigma_dE_proton_max;
		mass_target = MASS_NEUTRON;
		mode = 1;
	}


	//Select a value of E_k from dsig_i.
	//cout << dsig_i.Min() << " " << dsig_i.Max() << endl;

	//Placeholder values XMIN XMAX, YMIN YMAX until I've figured this out.
	double x,y,Q,Ek;
	double XMIN,XMAX,YMAX,YMIN;
	YMAX=XMAX=1;
	while(true){
		//These lower limits could be tightened to increase efficiency, not sure how to do it yet.
		Ek=dsig_i(Random::Flat(dsig_i.Min(),dsig_i.Max()));
		//Ek=Random::Flat(Efmin,Efmax);
		
		XMIN = std::max((Ek+mx-mass_target)/BeamEnergy,0.0);
		if(XMIN>XMAX){
			continue;
		}

		x=Random::Flat(XMIN,XMAX);
		double YMIN_Q = Y_Min_from_QMIN(BeamEnergy,MASS_PROTON,mass_target,x);
		//cout << "YMIN ESTIMATE=" << YMIN_Q << endl;

		YMIN = std::max(YMIN_Q,std::max((Ek+mx-BeamEnergy*x)/mass_target,0.0));
		//cout << "Ek=" << Ek;
		//cout << " XMIN=" << XMIN << " XMAX=" << XMAX << " X=" << x <<  " YMIN=" << YMIN << " YMAX=" << YMAX << endl;
		if(YMIN>YMAX){
			//cout << "YMIN>YMAX" << endl; 
			continue;
		}
		y=Random::Flat(YMIN,YMAX);
		//cout << "Q=" << sqrt(Q2_eval( BeamEnergy,  mass_target,  x, y)) << " Qmin=" << sqrt(Q2_eval( BeamEnergy,  mass_target,  x, YMIN)) << " Qmax=" << sqrt(Q2_eval( BeamEnergy,  mass_target,  x, YMAX)) << endl;

		if(Verify(Ek,mass_target,x,y,Q)&&dsig(Ek, x, y, Q, mode)>Random::Flat(0,dsig_max_i(Ek)))
		//if(Verify(Ek,mass_target,x,y,Q)&&dsig(Ek, x, y, Q, mode)>Random::Flat(0,pMax))
		{
			break;
		}
	}

	double theta = Theta_from_E3(x*BeamEnergy, x*MASS_PROTON, y*mass_target, Ek, mx);

	chi1.ThreeMomentumPolar(sqrt(Ek*Ek-mx*mx),theta,Random::Flat(0,2*pi));
	if(part.px!=0||part.py!=0){
		chi1.Rotate_y(part.Theta());
		chi1.Rotate_z(part.Phi());
	}
	chi2.ThreeMomentum(part.px*x-chi1.px,part.py*x-chi1.py,part.pz*x-chi1.pz);
	//cout << "theta " << theta << " x " << x << " y " << y << " Ek " << Ek << " mass_target " << mass_target <<  endl;
	//chi1.report(cout);
	//chi2.report(cout);

	vec.push_back(part);
 	double intersect1=det_int(chi1);
    double intersect2=det_int(chi2);

    //chi1.report(cout);
    //chi2.report(cout);

    if((intersect1)>0 || (intersect2)>0){;
        if(intersect1>0){
            vec.push_back(chi1);
        }
        if(intersect2>0){
            vec.push_back(chi2);
        }
        return true;
    }  
    else
        return false;
}
