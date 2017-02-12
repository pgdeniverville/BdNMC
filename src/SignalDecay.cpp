#include "SignalDecay.h"
#include "constants.h"
#include <cmath>
#include "Kinematics.h"
#include "decay.h"

using std::list;
using std::vector;
using std::cout; using std::endl;
using std::string;


double decay_probability(double t1, double t2, double lifetime){
    cout << "dec_prob_report\n" << "exp(-t1/tau)=" << exp(-t1/lifetime) << " exp(-t2/tau)=" << exp(-t2/lifetime) << " prob=" << exp(-t1/lifetime) - exp(-t2/lifetime) <<  endl;
    return abs(exp(-t2/lifetime) - exp(-t1/lifetime));
}

double generate_decay_time(double t1, double t2, double lifetime){
    return -log(exp(-t1/lifetime)-Random::Flat(0,1)/(exp(-t1/lifetime)-exp(-t2/lifetime)))*lifetime;
}

//Need to create local copies of the input vectors
SignalDecay::SignalDecay(double lifetime, std::vector<double> branching_ratios, std::vector<std::vector<Particle> > final_states){
    Branching_Ratios = branching_ratios;
    Final_States = final_states;
    Lifetime = lifetime;
    pMax=0;
    for(unsigned i = 0; i < Final_States.size();i++){
        if(Final_States[i].size()==0){
            Channel_Name.push_back("Undetermined Products");
            continue;
        }
        
        string tmp = Final_States[i][0].name;
        tmp += " ";
        for(unsigned j = 1; j < Final_States[i].size(); j++){
            cout << Final_States[i][j].name << endl;
            tmp += Final_States[i][j].name;
        }
        Channel_Name.push_back(tmp);
    }
}

bool SignalDecay::probscatter(std::shared_ptr<detector>& det, list<Particle>& partlist, list<Particle>::iterator& Parentit){
    double time1 = (*Parentit).Momentum()*det->cross_point[0]/(*Parentit).Speed()/speed_of_light;
    double time2 = (*Parentit).Momentum()*det->cross_point[1]/(*Parentit).Speed()/speed_of_light;
    if(time1>time2){
        double ttmp = time1;
        time1 = time2;
        time2 = ttmp;
    }
    cout << "Lifetime=" << Lifetime*1.0/sqrt(1-pow(Parentit->Speed(),2)) << " t1=" << time1 << " t2=" << time2 << endl;

    double prob = decay_probability(time1, time2, Lifetime*1.0/sqrt(1-pow(Parentit->Speed(),2)));
    double u = Random::Flat(0,1);
    if(prob>u*pMax){
        cout << "u=" << u << " pMax*u=" << u*pMax << endl; 
        if(prob>pMax){
            pMax=prob;
        }
    
        double timegen = generate_decay_time(time1, time2, Lifetime*1.0/sqrt(1-pow(Parentit->Speed(),2))); 
        Parentit->Set_Time(timegen);
        
        double br_chan = Random::Flat(0,1);
        unsigned i;
        
        for(i = 0; i<Branching_Ratios.size();i++){
            if(br_chan<Branching_Ratios[i])
                break;
            br_chan -= Branching_Ratios[i];
        }

        if(Final_States[i].size() == 2){
            Particle daughter1(Final_States[i][0]);
            Particle daughter2(Final_States[i][1]);
            

            TwoBodyDecay(*Parentit, daughter1, daughter2);

            partlist.insert(std::next(Parentit),daughter1);
            partlist.insert(std::next(Parentit),daughter2);
        }
        else{
            Particle product(*Parentit);

            product.name = Channel_Name[i];
            
            partlist.insert(std::next(Parentit),product);
        }

        return true;
    }

    return false;
}

bool SignalDecay::probscatter(std::shared_ptr<detector>& det, Particle &Parent){
    double time1 = Parent.Momentum()*det->cross_point[0]/Parent.Speed()/speed_of_light;
    double time2 = Parent.Momentum()*det->cross_point[1]/Parent.Speed()/speed_of_light;

    double prob = decay_probability(time1, time2, Lifetime);

    if(prob>Random::Flat(0,1)*pMax){
        if(prob>pMax)
            pMax=prob;
        return true;
    } 

    return false;
}
