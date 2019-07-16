#include "SignalDecay.h"
#include "constants.h"
#include <cmath>
#include "Kinematics.h"
#include "decay.h"

using std::list;
using std::vector;
using std::cout; using std::endl;
using std::string;

using namespace std::placeholders;

double decay_probability(double t1, double t2, double lifetime){
    //cout << "dec_prob_report\n" << "exp(-t1/tau)=" << exp(-t1/lifetime) << " exp(-t2/tau)=" << exp(-t2/lifetime) << " prob=" << exp(-t1/lifetime) - exp(-t2/lifetime) <<  endl;
    return abs(exp(-t2/lifetime) - exp(-t1/lifetime));
}

double generate_decay_time(double t1, double t2, double lifetime){
    return -log(exp(-t1/lifetime)-Random::Flat(0,1)*(exp(-t1/lifetime)-exp(-t2/lifetime)))*lifetime;
}

//Need to create local copies of the input vectors
SignalDecay::SignalDecay(double lifetime, std::vector<double> branching_ratios, std::vector<std::vector<Particle> > final_states){
    Branching_Ratios = branching_ratios;
    Final_States = final_states;
    Lifetime = lifetime;
    cout << "V_lifetime: " << Lifetime << endl; 
    pMax=0;
    for(unsigned i = 0; i < Final_States.size();i++){
        if(Final_States[i].size()==0){
            Channel_Name.push_back("Undetermined Products");
            continue;
        }
        
        string tmp = Final_States[i][0].name;
        tmp += " ";
        for(unsigned j = 1; j < Final_States[i].size(); j++){
            //cout << Final_States[i][j].name << endl;
            tmp += Final_States[i][j].name;
        }
        Channel_Name.push_back(tmp);
    }
}


bool SignalDecay::probscatter(std::shared_ptr<detector>& det, list<Particle>& partlist, list<Particle>::iterator& Parentit){

    double time1 = (*Parentit).Momentum()*(Parentit->crossing[0])/(*Parentit).Speed()/speed_of_light;
    double time2 = (*Parentit).Momentum()*(Parentit->crossing[1])/(*Parentit).Speed()/speed_of_light;
    if(time1>time2){
        double ttmp = time1;
        time1 = time2;
        time2 = ttmp;
    }
    //cout << "Lifetime=" << Lifetime*1.0/sqrt(1-pow(Parentit->Speed(),2)) << " t1=" << time1 << " t2=" << time2 << endl;

    double prob = decay_probability(time1, time2, Lifetime*boost_calc(Parentit->m,Parentit->E));
    double u = Random::Flat(0,1);
    //cout << prob << " " << pMax << " " << u*pMax << endl;
    if(prob>u*pMax){
        Parentit->EVENT_SET=true;
        //cout << "u=" << u << " pMax*u=" << u*pMax << endl; 
        //cout << "Momentum " << Parentit->Momentum() << " speed " << Parentit->Speed() << " crossing1 " << Parentit->crossing[0] << " crossing2 " << Parentit->crossing[1] << " pos1 " << Parentit->crossing[0]*(Parentit->Momentum()) << endl;
        //cout << "Lifetime=" << Lifetime*1.0/sqrt(1-pow(Parentit->Speed(),2)) << " t1=" << time1 << " t2=" << time2 << " prob " << prob << endl;
        if(prob>pMax){    
            pMax=prob;
        }
    
        double timegen = generate_decay_time(time1, time2, Lifetime*boost_calc(Parentit->m,Parentit->E)); 
        Parentit->Increment_Time(timegen);
        
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
            daughter1.EVENT_SET=true; 
            daughter2.EVENT_SET=true;
            //Parentit->Generate_Position();
            TwoBodyDecay(*Parentit, daughter1, daughter2);
            double angular_spread = Angle_Spread(daughter1.px,daughter1.py,daughter1.pz,daughter2.px,daughter2.py,daughter2.pz);  
            if(angular_spread<min_angle||angular_spread>max_angle||daughter1.E<Escatmin||daughter2.E<Escatmin){
                return false;
            }
            partlist.insert(std::next(Parentit),daughter1);
            partlist.insert(std::next(Parentit),daughter2);
        }
        else{
            Particle product(*Parentit);
            product.EVENT_SET=true;
            product.name = Channel_Name[i];
            
            partlist.insert(std::next(Parentit),product);
        }

        return true;
    }

    return false;
}

bool SignalDecay::probscatter(std::shared_ptr<detector>& det, Particle &Parent){
    double time1 = Parent.Momentum()*(Parent.crossing[0])/Parent.Speed()/speed_of_light;
    double time2 = Parent.Momentum()*(Parent.crossing[1])/Parent.Speed()/speed_of_light;
    if(time1>time2){
        double ttmp = time1;
        time1 = time2;
        time2 = ttmp;
    }
    double prob = decay_probability(time1, time2, Lifetime*boost_calc(Parent.m,Parent.E));

    //Parent.report(cout);
    //cout << "MOM: " << Parent.Momentum() << " crossing1: " << Parent.crossing[0] << " crossing2: " << Parent.crossing[1] << " Speed: " << Parent.Speed() << endl; 
    //cout << "Lifetime: " << Lifetime << " Lifetime_Actual: " << Lifetime*Parent.E/Parent.m << " time1: " << time1 << " time2: " << time2 << " pos1: " << speed_of_light*time1 << " pos2: "<< speed_of_light*time2 << " prob: " << prob << endl;
    if(prob>Random::Flat(0,1)*pMax){
        //cout << "Momentum " << Parent.Momentum() << " speed " << Parent.Speed() << " crossing1 " << Parent.crossing[0] << " crossing2 " << Parent.crossing[1] << " pos1 " << Parent.crossing[0]*Parent.Momentum() << " pos2 " <<  Parent.crossing[1]*Parent.Momentum() << endl;
        //cout << "Lifetime=" << Lifetime*1.0/sqrt(1-pow(Parent.Speed(),2)) << " t1=" << time1 << " t2=" << time2 << " prob " << prob << endl; 
        if(prob>pMax){
            pMax=prob;
        }
        return true;
    } 

    return false;
}

/*************************
    SignalDecay_2
*************************/


SignalDecay_2::SignalDecay_2(double lifetime, std::vector<std::shared_ptr<DMGenerator> > dmgens){
    Lifetime = lifetime;
    Channels = dmgens;
    cout << dmgens[0]->Channel_Name() << endl;
    pMax=0;
    for(unsigned i = 0; i< Channels.size();i++){
        br_total+=(Channels[i])->BranchingRatio();
        Channel_Name.push_back((Channels[i])->Channel_Name());
    }
}

bool SignalDecay_2::probscatter(std::shared_ptr<detector>& det, list<Particle>& partlist, list<Particle>::iterator& Parentit){

    if(Parentit->E<Escatmin)
        return false;

//    cout << "Begin probscatter for " << Parentit->name << endl;

    double time1 = (*Parentit).Momentum()*(Parentit->crossing[0])/(*Parentit).Speed()/speed_of_light;
    double time2 = (*Parentit).Momentum()*(Parentit->crossing[1])/(*Parentit).Speed()/speed_of_light;
    if(time1>time2){
        double ttmp = time1;
        time1 = time2;
        time2 = ttmp;
    }

//    Parentit->report();
//    cout << "Energy=" << Parentit->E << " Lifetime=" << Lifetime*boost_calc(Parentit->m,Parentit->E) << " Speed=" << Parentit->Speed() << " Boost=" << boost_calc(Parentit->m,Parentit->E) << " t1=" << time1 << " t2=" << time2 << endl;
    double prob = decay_probability(time1, time2, Lifetime*boost_calc(Parentit->m,Parentit->E));
//    cout << "prob=" << prob << endl;
//    cout << "get_prob()=" << get_pMax() << endl;
    double u = Random::Flat(0,1);
//    cout << prob << " " << pMax << " " << u*pMax << endl;


    if(prob>u*pMax){
        Parentit->EVENT_SET=true;
 
        if(prob>pMax){
            pMax=prob;
        }

        double timegen = generate_decay_time(time1, time2, Lifetime*boost_calc(Parentit->m,Parentit->E)); 

        Parentit->Increment_Time(timegen);
        Parentit->EVENT_SET=true;

/*
        cout << "Decay!\n";
        cout << "decay timegen=" << timegen << endl;
        Parentit->report();
*/
        double br_chan = Random::Flat(0,br_total);
        unsigned i;

        for(i = 0; i<Channels.size();i++){
            if(br_chan<(Channels[i])->BranchingRatio()){
                break;
            }
            br_chan -= (Channels[i])->BranchingRatio();
        }

        std::function<bool(Particle&)> det_int = bind(&detector::Ldet,det,_1);

        Channels[i]->GenDM(partlist, det_int, *Parentit);
/*
        cout << "Reporting partlist" << endl;
        for(list<Particle>::iterator it = partlist.begin(); it != partlist.end(); it++){
            it->report();
        }
        cout << "End partlist\n\n";
*/
        //No cuts yet. Not sure how to implement them. I'll do post-processing for now.
/*
        if(Final_States[i].size() == 2){
            Particle daughter1(Final_States[i][0]);
            Particle daughter2(Final_States[i][1]);
            daughter1.EVENT_SET=true; 
            daughter2.EVENT_SET=true;
            //Parentit->Generate_Position();
            TwoBodyDecay(*Parentit, daughter1, daughter2);
            double angular_spread = Angle_Spread(daughter1.px,daughter1.py,daughter1.pz,daughter2.px,daughter2.py,daughter2.pz);  
            if(angular_spread<min_angle||angular_spread>max_angle||daughter1.E<Escatmin||daughter2.E<Escatmin){
                return false;
            }
            partlist.insert(std::next(Parentit),daughter1);
            partlist.insert(std::next(Parentit),daughter2);
        }
        else{
            Particle product(*Parentit);
            product.EVENT_SET=true;
            product.name = Channel_Name[i];
            
            partlist.insert(std::next(Parentit),product);
        }
*/
        return true;
    }

    return false;
}

bool SignalDecay_2::probscatter(std::shared_ptr<detector>& det, Particle &Parent){
    if(Parent.E<Escatmin)
        return false;

    double time1 = Parent.Momentum()*(Parent.crossing[0])/Parent.Speed()/speed_of_light;
    double time2 = Parent.Momentum()*(Parent.crossing[1])/Parent.Speed()/speed_of_light;
    if(time1>time2){
        double ttmp = time1;
        time1 = time2;
        time2 = ttmp;
    }
    double prob = decay_probability(time1, time2, Lifetime*boost_calc(Parent.m,Parent.E));

    //Parent.report(cout);
    //cout << "MOM: " << Parent.Momentum() << " crossing1: " << Parent.crossing[0] << " crossing2: " << Parent.crossing[1] << " Speed: " << Parent.Speed() << endl; 
    //cout << "Lifetime: " << Lifetime << " Lifetime_Actual: " << Lifetime*Parent.E/Parent.m << " time1: " << time1 << " time2: " << time2 << " pos1: " << speed_of_light*time1 << " pos2: "<< speed_of_light*time2 << " prob: " << prob << endl;
    if(prob>Random::Flat(0,1)*pMax){
        //cout << "Momentum " << Parent.Momentum() << " speed " << Parent.Speed() << " crossing1 " << Parent.crossing[0] << " crossing2 " << Parent.crossing[1] << " pos1 " << Parent.crossing[0]*Parent.Momentum() << " pos2 " <<  Parent.crossing[1]*Parent.Momentum() << endl;
        //cout << "Lifetime=" << Lifetime*1.0/sqrt(1-pow(Parent.Speed(),2)) << " t1=" << time1 << " t2=" << time2 << " prob " << prob << endl; 
        if(prob>pMax){
            //Parent.report();
            //cout << "Energy=" << Parent.E << " Lifetime=" << Lifetime*boost_calc(Parent.m,Parent.E) << " Speed=" << Parent.Speed() << " Boost=" << boost_calc(Parent.m,Parent.E) << " t1=" << time1 << " t2=" << time2 << endl;
            //cout << prob << " " << pMax << endl;
            pMax=prob;
        }
        return true;
    } 

    return false;
}