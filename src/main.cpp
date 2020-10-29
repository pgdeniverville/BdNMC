#include <iostream>
#include <iomanip>
#include <fstream>
#include <vector>
#include <list>
#include <string>
#include <functional>
#include <cstdlib>
#include <exception>
#include <numeric>
#include <memory>
#include <climits>

#include "record.h"

#include "constants.h"

#include "Distribution.h"
#include "sanfordwang.h"
#include "BurmanSmith.h"
#include "BMPT_dist.h"
#include "Model.h"

#include "Integrator.h"
#include "detector.h"

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
#include "SignalDecay.h"


//Plotting stuff
//#include "DMNscattering.h"
//#include "DMscattering.h"

using std::cout;    using std::endl;
using std::vector;  using std::string;
using std::bind;    using std::function;
using std::list;    using std::vector;
using std::exception; using std::next;
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
    
/**************************************
	Making plots!
**************************************/
	/*
	double beam_energy2=400; double epsilon=1; double alphaD=1; double ptmax=1;
	std::string model = "proton_brem";
	for(double mv = 0.001; mv<=2; mv+=0.001){
		double zmin = std::max(3*mv/beam_energy2,0.1);
		double zmax = 1-zmin;
		Proton_Brem_Distribution pbd(beam_energy2,epsilon,mv,ptmax,zmax,zmin,alphaD,model);
		cout << mv << " " <<  pbd.V_prod_rate() << " " << zmin << " " << zmax << endl;
	}
	return 0;
*/

    //std::ofstream logging("log2.dat", std::ofstream::out);
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

    //Signal Channel
	string sigchoice = par->Signal_Channel();
    //Particles with this name will be checked for intersection with the
    //detector and passed to Scatter.probscatter.
    
    string sig_part_name = "DM";
    vector<string> sig_part_vec; 
	
	string outmode; 
	std::unique_ptr<std::ofstream> comprehensive_out; 

	//if((outmode=par->Output_Mode())=="summary")
	if((outmode=par->Output_Mode())=="comprehensive"||outmode=="dm_detector_distribution"){
		comprehensive_out = std::unique_ptr<std::ofstream>(new std::ofstream(output_filename));
    	if(!comprehensive_out->is_open()){
        	cout << "Unable to open output file: " << output_filename << endl;
        	parstream.close();
        	return -1;
    	}
	}
    
	parstream.close();
    
	//Run Parameters
    int samplesize = par->Sample_Size();
    int repeat = 1;//The number of times to attempt a scattering. Currently does work with comprehensive.
    double POT = par->Protons_on_Target();

	//Detector setup
    std::shared_ptr<detector> det = std::shared_ptr<detector>(par->Get_Detector());
    function<double(Particle&)> det_int = bind(&detector::Ldet,det,_1);//should get rid of this eventually, just pass detector object references.

    // return 0;
	//Production Mode
	vector<std::shared_ptr<DMGenerator> > DMGen_list;//Should be a unique_ptr
	vector<std::shared_ptr<Distribution> > PartDist_list;//This is an abstract class from which all distributions inherit.
	vector<double> Vnum_list;
	vector<string> proddist_vec;
	
	double Vnumtot=0;//This is a scaling that takes into account number of originator mesons produced.

	double beam_energy = par->Beam_Energy();
	
	std::shared_ptr<Model> model;
	std::shared_ptr<list<production_channel> > prodlist = par->Get_Production_List();	
	int chan_count = prodlist->size(); 
	
	std::shared_ptr<Scatter> SigGen;
	//Should add scatter angle to this as well at some point.	
	double max_scatter_energy = par->Max_Scatter_Energy();
	double min_scatter_energy = par->Min_Scatter_Energy();
	double min_angle = par->Min_Angle();
	double max_angle = par->Max_Angle();
    //	if((par->Model_Name())!="Dark_Photon"&&par->Model_Name()!="Dark_Photon_DM"){
	cout << "Setting up model " << par->Model_Name() << endl;
	if(par->Model_Name()=="Pseudoscalar_Mediator"){
		model=std::shared_ptr<Pseudoscalar>(new Pseudoscalar(*par));
	}
	else if(par->Model_Name()=="Axion_Dark_Photon"){
		model = std::shared_ptr<Axion_Dark_Photon>(new Axion_Dark_Photon(*par));
	}
	else if(par->Model_Name()=="Inelastic_Dark_Matter"){
		model = std::shared_ptr<Inelastic_Dark_Matter>(new Inelastic_Dark_Matter(*par));
	}
	//This will eventually replace Dark Photon
	else if(par->Model_Name()=="Dark_Photon" or par->Model_Name()=="Dark_Photon_DM"){
		model = std::shared_ptr<Kinetic_Mixing>(new Kinetic_Mixing(*par));
	}
    else if(par->Model_Name()=="Dark_Scalar"){
        model = std::shared_ptr<Scalar_Mediator>(new Scalar_Mediator(*par));
    }	
    model->Prepare_Model(*par);
	
	if(par->Output_Mode()=="particle_list"){
		return 0;
	}

	model->get_DMGen(DMGen_list);
	PartDist_list=model->get_Distribution();
	//This will be a list eventually, but SigGen is still a pointer to a single Signal Gen.
	model->get_first_SigGen(SigGen);
	model->get_Vnum(Vnum_list);
	Vnumtot=model->get_Vnumtot();
	sig_part_vec = model->get_sig_part_vec();
	
	cout << "Model assigned" << endl;
	//}
	//else{}//End of setup.

	if(Vnumtot<=0){
		cout << "No DM production expected. Terminating run.\n";
		return 0;
	}


	SigGen->set_angle_limits(min_angle, max_angle);
	SigGen->set_energy_limits(min_scatter_energy, max_scatter_energy);
	
    //Begin Run
	cout << "--------------------" << endl;	
    cout << "Run parameters:" << endl;	
    cout << "--------------------" << endl;
    cout << "Number of events to be generated = " << samplesize  << endl;
    //This will eventually be default.
    cout << "Model Name = " << par->Model_Name() << endl;
    model->Report_Model();


	for(int i = 0; i<chan_count; i++){
        cout << "Production-Channel " << i+1 << " = " << DMGen_list[i]->Channel_Name();
		if(DMGen_list[i]->query_off_shell()){
			cout << " in Off-Shell mode.\n";
		} 
		else
			cout << " in On-Shell mode.\n";
		cout << "Production Distribution " << i+1 << " = " << PartDist_list[i]->get_name() << endl;
		cout << "Branching Ratio "  << i+1 << " = " << DMGen_list[i]->BranchingRatio() << endl;
	   	cout <<	"Production events in channel " << i+1 << " = " << Vnum_list[i] << endl;
	}
    cout << "Signal Channel = " << sigchoice << endl;
 	cout << "Beam Energy = " << beam_energy << " GeV" << endl;

	cout << "Maximum Scattering Energy = " << max_scatter_energy << " GeV" << endl;
	cout << "Minimum Scattering Energy = " << min_scatter_energy << " GeV" <<  endl;
	cout << "Timing Cut = " << timing_cut << " s\n";

	std::unique_ptr<std::ofstream> summary_out;	
	summary_out = std::unique_ptr<std::ofstream>(new std::ofstream(summary_filename,std::ios::app));
		
	if(!summary_out->is_open()){
        cerr << "Unable to open output file: " << summary_filename << endl;
        parstream.close();
        return -1;
    }

	////////////////////
	//START OF BURN-IN//
	////////////////////
        
    double BURN_MAX = par->Burn_In();
	if(BURN_MAX<0){
        cout << "burn_max < 0 specified. Assuming default value burn_max=1000\n";
        BURN_MAX=1000;
    }
    if(outmode=="dm_detector_distribution"){
        BURN_MAX = 0;
        cout << "Detector_Mode selected.\nSkipping Burn-In.\n";
    }
    if(par->Weighted_Events()==true){	
    	BURN_MAX = 0;
    	cout << "Weighted Events selected.\nSkipping Burn-In.\n";
    }

    //This bit should be deprecated once everything moves over to model classes.
    if(sig_part_vec.empty()){
	    sig_part_vec.push_back(sig_part_name);
	}

    double BURN_OVERRIDE = par->Burn_Timeout();
    if(Vnumtot< par->Min_Event()){
    	cout << "Fewer than " << par->Min_Event() << " events possible, skipping BURN_IN and setting pMax=0 to skip run.\n";
    	BURN_MAX=0;
    }
    
	for(int i=0;BURN_MAX!=0&&i<chan_count; i++){
        int nburn = 0;
		if(Vnum_list[i]==0){
			cout << "Skipping Channel " << i+1 << ", no events expected.\n";
			continue;
		}
		else{
			cout << "Begin Channel " << i+1 << " Burn-In" << endl;
		}
		list<Particle>::iterator burniter;
		list<Particle>::iterator nextit;
        for(int burnattempt=0; (nburn < BURN_MAX)&&(burnattempt<BURN_MAX*BURN_OVERRIDE); burnattempt++){
            list<Particle> vecburn;
			Particle dist_part (0);
			PartDist_list[i]->Sample_Particle(dist_part);
			//cout << "det_int " << i << " = " << det_int(dist_part) << endl;
			if(DMGen_list[i]->GenDM(vecburn, det_int, dist_part)){
				//cout << "Got a hit!\n";
                burniter = vecburn.begin();
				for(nextit=next(burniter); burniter != vecburn.end(); burniter=nextit++){
//                    burniter->report();
                    if(std::find(sig_part_vec.begin(),sig_part_vec.end(), burniter->name)!=sig_part_vec.end()){
                        SigGen->probscatter(det, *burniter);
						nburn++;
						//burniter=next;
					}
				}
			}
		}


		if(nburn < BURN_MAX){
			cout << "Burn-In timed out, scattering number=" << nburn << " scattering goal=" << BURN_MAX << endl;
		}
		cout << "pMax at end of Burn-In = " << SigGen->get_pMax() << endl;
	}
    
	///////////////////
    //SIMULATION LOOP//
 	///////////////////
    cout << "Run " << par->Run_Name()  << " Start" << endl;


    if(outmode=="comprehensive"){
		*comprehensive_out << "Run " << par->Run_Name() << endl;
    }	

    long long int trials = 0;
    vector<long> trials_list(chan_count,0);
	vector<long double> scat_list(chan_count,0.0);
    int nevent=0;
    vector<int> NDM_list(chan_count,0);
//    vector<double> timing_efficiency(chan_count,0.0);
	
	//int escat=0;lso thought it could have used Leatherhead, or some of the other mutanimals. He probably would have fit right in with the Scale tail clan.
    bool scatter_switch;
    int trials_max = par->Max_Trials();

    if(trials_max>0){
	    cout << "Maximum Trials set to " << trials_max << endl;
    }

    if((SigGen->get_pMax()<=0 || SigGen->get_pMax()*Vnumtot<=par->Min_Event()) && outmode!="dm_detector_distribution" && !par->Weighted_Events()){
        cout << "pMax less than tolerance limit, skipping remainder of run\n";
    }
    else{
    	//cout << "Beginning loop\n" << endl;
        for(; (nevent < samplesize) && ((trials < trials_max)||(trials_max<=0)); trials++){
            int i;
            scatter_switch = false;
            double vrnd = Random::Flat(0.0,1)*Vnumtot;
            for(i=0; i<chan_count; i++){
                //cout << i << " vrnd=" <<  vrnd << " vs Vnum_list " << Vnum_list[i] <<  endl;
                if(vrnd<=Vnum_list[i]){
                    break;
                }
                else{
                    vrnd-=Vnum_list[i];
                }
            }
            trials_list[i]++; 
            list<Particle> vec;
            Particle dist_part (0);
            PartDist_list[i]->Sample_Particle(dist_part);
            
            list<Particle>::iterator iter;
            list<Particle>::iterator nextit;

            //cout << "nevent = " << nevent << endl;

            if(DMGen_list[i]->GenDM(vec, det_int, dist_part)){
                //Yes, this list is named vec.
                iter = vec.begin();
                for(nextit=next(iter) ; iter != vec.end();iter=nextit++){
                //The way this is structured means I can't do my usual repeat thing to boost stats.
                    if(std::find(sig_part_vec.begin(),sig_part_vec.end(), iter->name)!=sig_part_vec.end()){
                        NDM_list[i]++;
                        if(outmode=="dm_detector_distribution"){
                            *comprehensive_out << DMGen_list[i]->Channel_Name() << " " << det->Ldet(*iter) << " ";
                            iter->report(*comprehensive_out);
                            scatter_switch=true;
                            continue;;
                        }
                        //may need to replace this with a list<Particle> later
                        //cout << SigGen->get_pMax() << endl;
                        //If pMax=0, all events are accepted.
                        if(par->Weighted_Events()){
                        	SigGen->set_pMax(0);
                        }
                        if(SigGen->probscatter(det, vec, iter)){
                            //cout << "Scatter?\n"; 
                            //cout << "Made it inside the probscatter if\n";
                            //cout << "prob = " << SigGen->get_pMax() << endl;
                            double timing_prob_factor;
							if(timing_cut>0){
                                timing_prob_factor=t_delay_fraction(timing_cut,sqrt(pow(iter->end_coords[0],2)+pow(iter->end_coords[1],2)+pow(iter->end_coords[2],2)),iter->Speed());
                            }
                            else{
                            	timing_prob_factor=1;
                            }

                            if(par->Weighted_Events()){
                            	scat_list[i]+=SigGen->get_pMax()*timing_prob_factor;
                            }
                            else{
                            	scat_list[i]+=timing_prob_factor;
                            }
                            scatter_switch = true;	
                        }
                        else{
                            iter = vec.erase(iter);
                        }
                    }    
                }
            }
            if(scatter_switch&&outmode=="comprehensive"){
                *comprehensive_out << "event " << ++nevent;
                if(par->Weighted_Events()){
                	*comprehensive_out << " " << SigGen->get_pMax();
                }
                *comprehensive_out << endl;
                Record_Particles(*comprehensive_out, vec);
                *comprehensive_out << "endevent " << nevent << endl << endl;    
            }
            else if(scatter_switch){
                ++nevent;
	        }
        } 
    }
	cout << "Run complete\n";

    //This will be moved later.

    if(outmode=="summary"||outmode=="dm_detector_distribution"||
            outmode=="comprehensive"){
	    *summary_out << "Run " << par->Run_Name() << endl;
    }

	vector<double> signal_list(chan_count,0.0);
 	double signal = 0.0;
	int scattot = 0;
	int NDM = 0;    

    cout << chan_count << endl;

	for(int i=0; i<chan_count; i++){
		if(scat_list[i]==0)
			signal_list[i]=0;
		else{
			//cout << "Quick Diagnostic ";
			//cout << scat_list[i] << " " << trials << " " << Vnumtot << " " << SigGen->get_pMax() << " " << repeat << " " << par->Efficiency() << endl;
            if(par->Weighted_Events()){
                signal_list[i] = (double)scat_list[i]/(double)trials*Vnumtot/repeat*par->Efficiency();
            }
            else{
                signal_list[i] = (double)scat_list[i]/(double)trials*Vnumtot*SigGen->get_pMax()/repeat*par->Efficiency();
            }
		}
		scattot+=scat_list[i];
  		cout << DMGen_list[i]->Channel_Name() << ": " << signal_list[i];
		if(par->Weighted_Events()){
			cout << " Weight: " << scat_list[i] << " V_num: " << Vnum_list[i] << endl;
		}
		else{
			cout << " Events: " << scat_list[i] << " V_num: " << Vnum_list[i] << endl;			
		}
	
		if(outmode=="summary"||outmode=="dm_detector_distribution"||
                outmode=="comprehensive"){
			*summary_out << DMGen_list[i]->Channel_Name() << " ";
			model->Report(*summary_out);
			*summary_out << signal_list[i] << " " << sigchoice << " " << POT << " " << par->Efficiency() << " " << samplesize << " " << Vnum_list[i] << " " << Vnumtot << endl;
        }
        NDM+=NDM_list[i]; 
		signal+=signal_list[i];
 	}


    /*if(outmode=="dm_detector_distribution"){
        *comprehensive_out << "Total " << trials << " " << POT << " " << Vnumtot << " " << samplesize << " " << (double)NDM/(2*trials) << " " << endl;
    }*/

	if(outmode=="summary"||outmode=="comprehensive"){
		*summary_out << "Total ";
		model->Report(*summary_out);
		*summary_out << signal << " " << sigchoice << " " << POT << " " << par->Efficiency() << " " << samplesize << " " << endl;
    }
//    if(outmode=="summary"||outmode=="comprehensive"){
//        *summary_out << "Total " << mv  <<  " "  << mdm << " " << signal << " " << kappa << " " << alD << " " << sigchoice << " " << POT << " " << par->Efficiency() << " " << samplesize << " " << endl;
//    }
    else if(outmode=="dm_detector_distribution"){
    	cout << "Outputting summary" << endl;
    	*summary_out << "Total ";
		model->Report(*summary_out);
		*summary_out << signal << " " << sigchoice << " " << POT << " " << par->Efficiency() << " " << Vnumtot << " " << samplesize << " " << (double)NDM/(2*trials) << endl;
   
    	*comprehensive_out << "Total ";
		model->Report(*comprehensive_out);
		*comprehensive_out << signal << " " << sigchoice << " " << POT << " " << par->Efficiency() << " " << Vnumtot << " " << samplesize << " " << (double)NDM/(2*trials) << endl;
    }

    if(outmode =="comprehensive"){
		*comprehensive_out << endl << " Summary " << signal << " ";
		model->Report(*comprehensive_out);
		*comprehensive_out << sigchoice << " " << POT << " " <<    par->Efficiency() << " " << samplesize << " " << endl;
    }

    summary_out->close();

  	if(outmode=="dm_detector_distribution"||outmode=="comprehensive")
		comprehensive_out->close();

    cout << "Number of trials = " << trials << endl;	
	cout << "Number of candidates intersecting detector = " << NDM << endl;	
	cout << "Number of " << sigchoice <<  " = " << scattot << endl;
	for(int i=0; i<chan_count;i++)
		cout << "Number of events from channel " << i+1 << " " << DMGen_list[i]->Channel_Name() << " = " << scat_list[i] << " in " << trials_list[i] << " trials." << endl;  
	cout << "Number of sample events generated = " << nevent << endl;
	cout << "Acceptance = " << (double)NDM/(2*trials) << endl;	
 	if(!par->Weighted_Events()){
		cout << "Maximum event probability = "  << SigGen->get_pMax() << endl;
		cout << "Average event probability = "  << (double)scattot/NDM*SigGen->get_pMax() << endl;
	}
	else{
		cout << "Average event probability = " << scattot/(double)NDM << endl;
	}
	cout << "Predicted number of signal events = "  << signal << endl;
	for(int i=0; i<chan_count;i++)
		cout << "Predicted number of signal events from channel " << i+1 << " " << DMGen_list[i]->Channel_Name()  <<  " = "  << signal_list[i] << endl;
	
	if(outmode=="comprehensive"||outmode=="dm_detector_distribution"){
		cout << "--------------------" << endl;	
		cout << "Events stored in file " << output_filename << endl;	
		cout << "--------------------" << endl;
	}
	cout << "--------------------" << endl;	
	cout << "Summary appended to file " << summary_filename << endl;	
	cout << "--------------------" << endl; 
    delete(par);
    return 0;
 }
