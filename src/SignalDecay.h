#ifndef GUARD_SignalDecay_h
#define GUARD_SignalDecay_h

#include "detector.h"
#include "Particle.h"
#include "Random.h"
#include "Scatter.h"

#include <iostream>
#include <vector>
#include <list>
#include <string>

/*
 * SignalDecay calculates the signal from decays inside of a defined decay volume.
 * Currently does not support chains of decays. Is there any way to add it?
 * 
 * Simplest hack is to create multiple SignalDecay objects, with different decay
 * end states based on the possible parent particles.
 *
 * Currently only supports isotropic two-body final states. Three body final
 * states may require a derived class.
 *
 * Could this instead be handled by Particle itself? It's a lot of overhead to
 * shove into the Particle object.
 *
 * While this inherits the Scatter class, it is not a perfect fit for it.
 * It does not use many of its functions. The Scatter class may need to be
 * revised.
 *
 * This does not make use of Escatmax, Escatmin, min_angle or max_angle.
 */


/*
 * An idea for class to handle each decay channel. Maybe not worth it?
class DecayChannel{
    public:
        DecayChannel(std::vector<Particle> &daughters, double branching);
        void SampleDecay(Particle &daughter);
        double Branching_Ratio(){return _Branching_Ratio;}

    private:
        std::vector<Particle> Daughters;
        double _Branching_Ratio;
};
*/

class SignalDecay: public Scatter{
    public:
        //Creating a local copy
        SignalDecay(double lifetime, std::vector<double> branching_ratios, std::vector<std::vector<Particle> > final_states);
        ~SignalDecay(){};

        bool probscatter(std::shared_ptr<detector>& det, std::list<Particle>& partlist, std::list<Particle>::iterator&);
        bool probscatter(std::shared_ptr<detector>& det, Particle &Parent);

        void set_Model_Parameters(double lifetime, std::vector<double> branching_ratios, std::vector<std::vector<Particle> >);

       //void Generate_Position(std::shared_ptr<detector>& det, Particle &DM, std::vector<Particle> &scat); 
        double get_pMax(){return pMax;}
    private:
        std::vector<std::vector<Particle> > Final_States; 
        std::vector<double> Branching_Ratios;
        double Lifetime;
        std::vector<std::string> Channel_Name;
        //std::vector<int> decay_counts;
};


//New and improved signal decay channel
class SignalDecay_2: public Scatter{
    public:
        //The DMGenerators should be 2 or 3 body decay. Abusing it to generate those end states and add them to the partlist!
        //At some point I could implement cuts on a per channel basis, would need to write some kind of cut class that took particle lists? 
        //REMEMBER: Set the lifetime of DMGenerator objects to 0, so they decay immediately. The decay time was already calculated by SignalDecay_2!
        SignalDecay_2(double lifetime, std::vector<std::shared_ptr<DMGenerator> >);
        ~SignalDecay_2(){};
        bool probscatter(std::shared_ptr<detector>& det, std::list<Particle>& partlist, std::list<Particle>::iterator&);
        bool probscatter(std::shared_ptr<detector>& det, Particle &Parent);

        //The branching ratio of the final states we're looking for may be smaller than one!
        double get_pMax(){return br_total*pMax;}

        void set_Model_Parameters(double lifetime, std::vector<std::string> chan_names, std::vector<std::vector<Particle> >);

       //void Generate_Position(std::shared_ptr<detector>& det, Particle &DM, std::vector<Particle> &scat); 

    private: 
        double Lifetime=0;
        double br_total=0;
        std::vector<std::string> Channel_Name;
        std::vector<std::shared_ptr<DMGenerator> > Channels;
        double min_energy;
        //std::vector<int> decay_counts;
};

#endif
