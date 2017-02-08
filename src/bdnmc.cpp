#include <iostream>
#include <fstream>
#include <vector>
#include <list>
#include <string>
#include <functional>
#include <cstdlib>
#include <exception>
#include <numeric>
#include <memory>

#include "constants.h"

#include "Distribution.h"
#include "sanfordwang.h"
#include "BurmanSmith.h"
#include "BMPT_dist.h"

#include "Integrator.h"
#include "detector.h"
#include "record.h"
#include "Particle.h"
#include "Random.h"
#include "decay.h"
#include "branchingratios.h"
#include "Parameter.h"
#include "DMgenerator.h"
#include "Scatter.h"
#include "Particle_List.h"
#include "partonsample.h"
#include "Proton_Brem_Distribution.h"
#include "Position_Distributions.h"


//Plotting stuff
//#include "DMNscattering.h"
//#include "DMscattering.h"

using std::cout;    using std::endl;
using std::vector;  using std::string;
using std::bind;    using std::function;
using std::list;    using std::vector;
using std::exception;
using std::cerr;

//const double mp = MASS_PROTON;
//const double mn = MASS_NEUTRON;
//const double me = MASS_ELECTRON;
//const double EDMres = 0.1;
//cm per meter
//const double cmpm = 100.0;

double t_delay_fraction(double tcut, double pos, double speed){
	double tdelay=pos/speed/speed_of_light-pos/speed_of_light;
	if(tdelay>tcut)
		return 1;
	else
		return 1 - (tcut-tdelay)/tcut;
}

int main(int argc, char* argv[]){

	using namespace std::placeholders;
    string parameter_file;

    if(argc==2)
        parameter_file = argv[1];
    else
        parameter_file = "parameter.dat";

    //-------------------
    // Read in parameters
    //-------------------
    std::ifstream parstream(parameter_file);
    Parameter * par;//should make this a unique_ptr
	try{
    par = new Parameter(parstream);
    } catch(int e_int){
        if(e_int==-1){
            cerr << "Could not open " << parameter_file << ". Terminating run." << endl; 
            parstream.close();      
            return -1;
        }
		if(e_int==-2){
			cerr << "No partic/le_list_file supplied for particle_list distribution. Terminating run." << endl; 
            parstream.close();      
            return -1;
		}
        else{
            cerr << "Something unexpected went wrong. Terminating run." << endl;
            parstream.close();      
            return -1;
        }
    }
    catch(exception &e){
       cerr << "Exception in " << e.what() << " thrown. Terminating run." << endl;
        parstream.close();      
       return -1; 
    }
    if(par->Integrity()!=0){
        cerr << "Missing parameters, terminating run. Please check previous error messages." << endl;
        parstream.close();
        return -1;
    }
	cout << "Parameter read successfully\n";
    //Initializing Random Number Generator
    if(par->Seed()>=0)
		Random(par->Seed());
	else
		Random();


	//Timing cut in seconds	
	double timing_cut = par->Timing_Cut();

    string output_filename = par->Output_File();
    if(output_filename == ""){
        output_filename = "../Events/events.dat";
    }
	string summary_filename = par->Summary_File();
    if(summary_filename == ""){
        summary_filename = "../Events/summary.dat";
    }


	string outmode; 
	std::unique_ptr<std::ofstream> comprehensive_out; 

	//if((outmode=par->Output_Mode())=="summary")
	if((outmode=par->Output_Mode())=="comprehensive"||outmode=="dm_detector_distribution"){
		comprehensive_out = std::unique_ptr<std::ofstream>(new std::ofstream(output_filename));
    	if(!comprehensive_out->is_open()){
        	cout << "Unable to open output file: " << output_filename << endl;
        	parstream.close();
        	return 1;
    	}
	}
    
	parstream.close();
    
	//Run Parameters
    int samplesize = par->Sample_Size();
    int repeat = 1;//The number of times to attempt a scattering. Currently does work with comprehensive.
    double POT = par->Protons_on_Target();
    double num_pi0 = par->Pi0_per_POT()*POT;//ratio should be set in parameter

    //Detector setup
    std::shared_ptr<detector> det = std::shared_ptr<detector>(par->Get_Detector());
    function<double(Particle)> det_int = bind(&detector::Ldet,det,_1);//should get rid of this eventually, just pass detector object references.
	//function<double(Particle)> det_int = [](Particle x){return 1.0;};
	
	//Production Mode
	vector<std::shared_ptr<DMGenerator> > Prod_list;//Should be a unique_ptr
	vector<std::shared_ptr<Distribution> > PartDist_list;//This is an abstract class from which all distributions inherit.
	vector<double> Vnum_list;
	vector<string> proddist_vec;

	double Vnumtot=0;//This is a scaling, indicates the total number of candidate particles produced.
                     //Many of these particles will not reach the detector. 
	
    double beam_energy = par->Beam_Energy();
	double target_n = par->Target_N_Num();
	double target_p = par->Target_P_Num();
    double target_p_cross = par->P_Cross();

    for(list<production_channel>::iterator proditer = prodlist->begin(); proditer!=prodlist->end(); proditer++){
		string prodchoice = proditer->Production_Channel(); 
		string proddist = proditer->Prod_Dist();
		std::shared_ptr<DMGenerator> DMGen;
		std::shared_ptr<Distribution> PartDist;
		double Vnum;
		cout << "Setting up distribution " << proddist << " for channel " << prodchoice << endl;
		//Don't delete these without good reason! They're used immediately!
		double p1,p2,p3;
        
        //Production Distribution setup. This can be handled in the main code for now.
        //Production is all Standard Model.
     	if(proddist=="particle_list"){
			bool set_pos = proditer->par_list_pos();
			if(set_pos){
				cout << "position off-sets were read from particle_list file\n";
			}
			std::shared_ptr<Particle_List> pl(new Particle_List(proditer->particle_list_file,set_pos));
			PartDist = pl;
        }
        //Meson Production Distributions.
        else if(proddist=="pi0_sanfordwang"||proddist=="k0_sanfordwang"){
			std::shared_ptr<sanfordwang> sw(new sanfordwang(proddist));
			sw->set_fit_parameters(*proditer);//This does nothing if no fit parameters have been set.
			PartDist = sw;
		}
		else if(proddist=="burmansmith"){
			std::shared_ptr<BurmanSmith> bs(new BurmanSmith(beam_energy,target_p));
			//beam_energy should be kinetic energy for this case
			PartDist = bs;
		}
		else if(proddist=="bmpt"){
			cout << "Energy = " << beam_energy << endl;
			cout << "Mass Number = " << target_n+target_p << endl;
			std::shared_ptr<BMPT> bmpt(new BMPT(beam_energy,target_n+target_p));
			cout << "BMPT burn-in complete\n";
			PartDist = bmpt;
		}
        //This should actually be a DMGen (Which needs to be renamed), with the distribution of Beam Protons as the distribution!
        else if(proddist=="proton_brem"||proddist=="proton_brem_baryonic"){
			if(proditer->ptmax()<0 || proditer->zmax() < 0 || proditer->zmax()<proditer->zmin() || proditer->zmin() < 0){
				cerr << "Invalid properties for production_distribution proton_brem." << endl;
				return -1;
			}
			std::shared_ptr<Proton_Brem_Distribution> pbd(new Proton_Brem_Distribution(beam_energy, kappa,mv,proditer->ptmax(),proditer->zmax(),proditer->zmin(),alD,proddist));
			//cout << "kappa = " << kappa << " mv = " << mv << " " << proditer->ptmax() << " " << proditer->zmax() << " " << proditer->zmin() << endl;
			Vnum = pbd->V_prod_rate()*POT;
			PartDist = pbd;
		}
		else if(outmode=="particle_list"){
			cerr << "Invalid or missing distribution " << proddist << " declared for particle_list outmode\n Terminating run.\n";
			return -1;
		}
        
           
    }
