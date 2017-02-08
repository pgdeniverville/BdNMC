#include "SignalDecay.h"
#include "constants.h"
#include <cmath>
#include "Kinematics.h"
#include "decay.h"

using std::list;
using std::vector;
using std::cout; using std::endl;

double decay_probability(double t1, double t2, double lifetime){
    return exp(-t1/lifetime) - exp(-t2/lifetime);
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
}

bool SignalDecay::probscatter(std::shared_ptr<detector>& det, list<Particle>& partlist, list<Particle>::iterator& Parentit){
    double time1 = (*Parentit).Momentum()*det->cross_point[0]/(*Parentit).Speed()/speed_of_light;
    double time2 = (*Parentit).Momentum()*det->cross_point[1]/(*Parentit).Speed()/speed_of_light;

    double prob = decay_probability(time1, time2, Lifetime);

    if(prob>Random::Flat(0,1)*pMax){
        if(prob>pMax)
            pMax=prob;
    
        unsigned i;

        double br_chan = Random::Flat(0,1);

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
            std::cerr << "This code does not currently support >2 body final states for decays." << endl;
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
