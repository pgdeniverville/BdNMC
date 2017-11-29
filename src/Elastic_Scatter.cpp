#include "Scatter.h"
#include "Kinematics.h"

#include <math.h>

using std::vector; using std::function;

using namespace elastic_scattering;

Elastic_Scatter::Elastic_Scatter(){
    chan_number=0;
}

void Elastic_Scatter::add_channel(Particle& end_state, function<double(double)>& cross_total, function<double(double)>& cross_max, function<double(double, double)>& dsig, double& num_density){
   cross_tot.push_back(cross_total);
   cross_maxima.push_back(cross_max);
   dsigma.push_back(dsig);
   number_density.push_back(num_density);
   recoil_particle.push_back(end_state);
   chan_number++;
}

//This creates Particle objects for the recoiling particles.
bool Elastic_Scatter::probscatter(std::shared_ptr<detector>& det, std::list<Particle>& partlist, std::list<Particle>::iterator& partit){
    Particle recoil(0);
	if(probscatter(det, *partit, recoil)&&(min_angle<=0||recoil.Theta()>min_angle)&&(max_angle>2*pi||recoil.Theta()<max_angle)){
        partit->Generate_Position();
        Link_Particles(*partit,recoil);
		partlist.insert(std::next(partit),nucleon);
        Particle partout(partit->m);
        partout.name = "Recoil_"+partit->name();
        partout.ThreeMomentum(partit->px-recoil.px,partit->py-recoil.py,partit->pz-recoil.pz);
        Link_Particles(*partit, partout);
        partlist.insert(std::next(partit),partout);
		return true;
    }
    return false
}

//part is the incoming particle, recoil is the target particle after scattering.
bool Elastic_Scatter::probscatter(std::shared_ptr<detector>& det, Particle &part, Particle &recoil){
    using namespace std::placeholders;
    double LXDet = convmcm*(det->Ldet(DM));
    double total=0; vector<double> prob;
    for(int i = 0; i < chan_number; i++){
       prob.push_back(LXDet*convGeV2cm2*cross_tot[i]*number_density[i]);
       total+=prob.back();
    }
    if(total > pMax*Random::Flat(0,1)){
        if(prob > pMax){
            pMax = prob
        }
        
        double u=Random::Flat(0,1);
        double p=0; 
        for(int i = 0; i < chan_number; i++){
            p+=prob[i]
            if(p/total > u){
                std::function<double(double)> Xsec = std::bind(disgma,part.E,_1);
                recoil = recoil_particle[i];
                if(scatmax(part.E) < scatmin(part.E, part.m, recoil.m)){
                    return false;
                }
                scatterevent(part,recoil,Xsec,cross_maxima[i]);
            }
        }

        return true;
    }

    return false
}

void Elastic_Scatter::scatterevent (Particle &part, Particle &recoil, std::function<double(double)> Xsec,Xstd::function<double(double)> Xmax){
    double Efmax = scatmax(part.E, part.m, recoil.m);
	double Efmin = scatmin(part.E, part.m, recoil.m); 
	double dsigmax = std::max(Xsec(Efmax),Xmax(part.E));
	double xe,thetaN,phiN,pN;
    while(true){
        xe = Random::Flat(0,1)*(Efmax-Efmin)+Efmin;//Recoil energy
        if(Xsec(xe)/dsigmax > Random::Flat(0,1)){
            thetaN = Ef_to_Theta(part.E,xe,part.m,recoil.m);
            phiN = Random::Flat(0,1)*2*pi;
            pN = sqrt(pow(xe,2)-pow(recoil.m,2));
            recoil.ThreeMomentum(pN*sin(thetaN)*cos(phiN),pN*sin(thetaN)*sin(phiN),cos(thetaN)*pN);
            recoil.Rotate_y(part.Theta());
            recoil.Rotate_z(part.Phi());
            break;
        }
    }
}

double Elastic_Scatter::scatmax(double Ei, double m1, double m2){
    return(std::min(Escatmax,E2fMax(Ei, m1,m2)));
}

double Elastic_Scatter::scatmax(double Ei, double m1, double m2){
    return(std::max(Escatmin,E2fMin(Ei,m1,m2)));
}
