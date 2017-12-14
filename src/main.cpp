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
#include "Model.h"

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
#include "Axion_Dark_Photon.h"
#include "SignalDecay.h"


//Plotting stuff
//#include "DMNscattering.h"
//#include "DMscattering.h"

using std::cout;    using std::endl;
using std::vector;  using std::string;
using std::bind;    using std::function;
using std::list;    using std::vector;
using std::exception;
using std::cerr;

//Temporary! I will clean this up in the big revision
const string dark_axion_signal_string = "Dark_Photon";

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

	
	//Model parameters
    double kappa = par->Epsilon();
    double alD = par->alD();
    double mv = par->MassDP();
    double mdm = par->MassDM();	

    //This will turn into a Model class later on.
    Axion_Dark_Photon adp(par);

	//Detector setup
    std::shared_ptr<detector> det = std::shared_ptr<detector>(par->Get_Detector());
    function<double(Particle&)> det_int = bind(&detector::Ldet,det,_1);//should get rid of this eventually, just pass detector object references.
	//function<double(Particle)> det_int = [](Particle x){return 1.0;};

	//Production Mode
	vector<std::shared_ptr<DMGenerator> > DMGen_list;//Should be a unique_ptr
	vector<std::shared_ptr<Distribution> > PartDist_list;//This is an abstract class from which all distributions inherit.
	vector<double> Vnum_list;
	vector<string> proddist_vec;
	
	double Vnumtot=0;//This is a scaling that takes into account number of originator mesons produced.

	double beam_energy = par->Beam_Energy();
	double target_n = par->Target_N_Num();
	double target_p = par->Target_P_Num();
    double target_p_cross = par->P_Cross();
	
	std::shared_ptr<Model> model;
	std::shared_ptr<list<production_channel> > prodlist = par->Get_Production_List();	
	int chan_count = prodlist->size(); 
	
	double EDMRES = par->EDM_RES();

	std::shared_ptr<Scatter> SigGen;
//Update this so it scales energy to Beam energy.
	double max_dm_energy = par->Max_DM_Energy();
	//Should add scatter angle to this as well at some point.	
	double max_scatter_energy = par->Max_Scatter_Energy();
	double min_scatter_energy = par->Min_Scatter_Energy();
	double min_angle = par->Min_Angle();
	double max_angle = par->Max_Angle();
	
	if(par->Model_Name()=="Pseudoscalar_Mediator"){
		std::shared_ptr<Pseudoscalar> mod(new Pseudoscalar(*par));
		//This will eventually get moved to the bottom, once the big else statements is moved into models.
		mod->get_DMGen(DMGen_list);
		mod->get_Distribution(PartDist_list);
		mod->get_first_SigGen(SigGen);
		mod->get_Vnum(Vnum_list);
		Vnumtot = mod->get_Vnumtot();
		model = mod;
	}
	else{
		for(list<production_channel>::iterator proditer = prodlist->begin(); proditer!=prodlist->end(); proditer++){
			string prodchoice = proditer->Production_Channel(); 
			string proddist = proditer->Prod_Dist();
			std::shared_ptr<DMGenerator> DMGen;
			std::shared_ptr<Distribution> PartDist;
			double Vnum;
			cout << "Setting up distribution " << proddist << " for channel " << prodchoice << endl;
			//Don't delete these without good reason! They're used immediately!
			double p1,p2,p3;
			if(proddist=="pi0_sanfordwang"||proddist=="k0_sanfordwang"){
				std::shared_ptr<sanfordwang> sw(new sanfordwang(proddist));
				sw->set_fit_parameters(*proditer);//This does nothing if no fit parameters have been set.
				PartDist = sw;
			}
			else if(proddist=="particle_list"){
				bool set_pos = proditer->par_list_pos();
				if(set_pos){
					cout << "position off-sets were read from particle_list file\n";
				}
				std::shared_ptr<Particle_List> pl(new Particle_List(proditer->particle_list_file,set_pos));
				PartDist = pl;
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
			else if(proddist=="parton_V"||proddist=="parton_V_baryonic"){
				std::string proton_file = proditer->Parton_V_Proton_File();
				std::string neutron_file = proditer->Parton_V_Neutron_File();
				std::shared_ptr<parton_sample> parsam(new parton_sample(proton_file,neutron_file,\
							target_p,target_n));
				Vnum = POT*parsam->production_proton_cross_section()*microbarn/target_p_cross*(target_p+\
						target_n);
				if(proddist=="parton_V"){
					Vnum*=pow(kappa,2);
				}
				else
					Vnum*=alD;
				PartDist = parsam;
			}
			else if(proddist=="proton_brem"||proddist=="proton_brem_baryonic"){
				if(proditer->ptmax()<proditer->ptmin() || proditer->zmax() < 0 || proditer->zmax()<proditer->zmin() || proditer->zmin() < 0){
					cerr << "Invalid properties for production_distribution proton_brem." << endl;
					return -1;
			}	
				std::shared_ptr<Proton_Brem_Distribution> pbd(new Proton_Brem_Distribution(beam_energy, kappa,mv,proditer->ptmax(),proditer->zmax(),proditer->zmin(),alD,proddist,proditer->ptmin()));
				//cout << "kappa = " << kappa << " mv = " << mv << " " << proditer->ptmax() << " " << proditer->zmax() << " " << proditer->zmin() << endl;
				Vnum = pbd->V_prod_rate()*POT;
				PartDist = pbd;
			}
			else if(outmode=="particle_list"){
				cerr << "Invalid or missing distribution " << proddist << " declared for particle_list outmode\n Terminating run.\n";
				return -1;
			}
	        cout << proddist << " " << prodchoice << endl; 
	        
	        //Signal_Decay HANDLING
	        if(par->Model_Name()=="Axion_Dark_Photon"||(par->Model_Name()=="Dark_Photon"&&sigchoice=="Signal_Decay")){
	            if(proddist=="proton_brem"){
	                DMGen = std::shared_ptr<DMGenerator>(new Do_Nothing_Gen("Dark_Photon_Bremsstrahlung", dark_axion_signal_string));
	                cout << DMGen->Channel_Name() << endl;
	            }
	            else if(prodchoice=="pi0_decay"||prodchoice=="eta_decay"){
	                Particle gamma(0);
	                gamma.name = "Photon";
	                Particle dp(mv);
	                dp.name = "Dark_Photon";
	                //Code repeat here. Need to reorganize to eliminate this.
	                if(proddist=="default"){
					    std::shared_ptr<sanfordwang> sw(new sanfordwang("pi0_sanfordwang"));
					    sw->set_fit_parameters(*proditer);
					    PartDist = sw; //std::shared_ptr<Distribution>(&sw);
	                }
	                if(prodchoice=="pi0_decay"){
	                    if(mv>mpi0){
	                        cerr << "Attempting to produce on-shell mediator with mass " << mv << " larger than pion mass " << mpi0 << ".";
	                        return -1;
	                    }
	                    DMGen = std::shared_ptr<DMGenerator>(new 
	                            Two_Body_Decay_Gen(brpi0toVgamma(mv,mdm,kappa,alD),
	                            MASS_PION_NEUTRAL,"Pion",dp,gamma));
	                    PartDist->set_mass(MASS_PION_NEUTRAL);
	                }
	                else if (prodchoice=="eta_decay"){
	                    if(mv>meta){
	                        cerr << "Attempting to produce on-shell mediator with mass " << mv << " larger than eta " << meta << ".";
	                        return -1;
	                    }
	                    DMGen = std::shared_ptr<DMGenerator>(new 
	                            Two_Body_Decay_Gen(bretatoVgamma(mv,mdm,kappa,alD),
	                            MASS_ETA,"Eta",dp,gamma));
	                    PartDist->set_mass(MASS_ETA);
	                }
	                if(proditer->Meson_Per_Pi0()<=0){
					    //default for eta
	                    Vnum = DMGen->BranchingRatio()*num_pi0/30.0;
	                }
				    else{
					    Vnum = DMGen->BranchingRatio()*num_pi0*
	                        (proditer->Meson_Per_Pi0());
	                }
	                DMGen->Set_Channel_Name(prodchoice);
	            }
	            else{
	                cerr << "Invalid production or distribution selected\n Axion_Dark_Photon is still in dev, so may react oddly.\n";
	                return -1;
	            }
	        }
	        else if(prodchoice=="pi0_decay"||prodchoice=="pi0_decay_baryonic"){	
				if(proddist=="default"){
					std::shared_ptr<sanfordwang> sw(new sanfordwang("pi0_sanfordwang"));
					sw->set_fit_parameters(*proditer);
					//sw->sample_momentum(p1,p2,p3);//Can I remove this line now?
					PartDist = sw; //std::shared_ptr<Distribution>(&sw);
				}
				PartDist->set_mass(mpi0);
				if(prodchoice=="pi0_decay"){
					DMGen = std::shared_ptr<DMGenerator>(new pion_decay_gen(mv, mdm, kappa, alD));
				}
				else{
					DMGen = std::shared_ptr<DMGenerator>(new pion_decay_gen_baryonic(mv, mdm, kappa, alD));
				}
				if(proditer->Meson_Per_Pi0()<=0)
					Vnum = DMGen->BranchingRatio()*num_pi0;
				else
					Vnum = DMGen->BranchingRatio()*num_pi0*(proditer->Meson_Per_Pi0());
			}
			else if(prodchoice=="eta_decay"||prodchoice=="eta_decay_baryonic"){
	            if(proddist=="default"){
					std::shared_ptr<sanfordwang> sw(new sanfordwang("k0_sanfordwang"));
					sw->set_fit_parameters(*proditer);
					PartDist = sw;
				}

				PartDist->set_mass(meta);
				
				if(prodchoice=="eta_decay"){ 
					DMGen = std::shared_ptr<DMGenerator>(new eta_decay_gen(mv, mdm, kappa, alD));
				}
				else{
					DMGen = std::shared_ptr<DMGenerator>(new eta_decay_gen_baryonic(mv, mdm, kappa, alD));
				}
				if(proditer->Meson_Per_Pi0()<=0)
					Vnum = DMGen->BranchingRatio()*num_pi0/30.0;
				else
					Vnum = DMGen->BranchingRatio()*num_pi0*(proditer->Meson_Per_Pi0());
			}
			else if(prodchoice=="omega_decay"||prodchoice=="rho_decay"||
	                prodchoice=="phi_decay"||prodchoice=="omega_decay_baryonic"||prodchoice=="phi_decay_baryonic"){
				if(proddist=="default"){
					std::shared_ptr<sanfordwang> sw(new sanfordwang("pi0_sanfordwang"));
					sw->set_fit_parameters(*proditer);
					sw->sample_momentum(p1,p2,p3);
					PartDist = sw; //std::shared_ptr<Distribution>(&sw);
				}
				if(prodchoice=="omega_decay"){
					PartDist->set_mass(momega);
					DMGen = std::shared_ptr<DMGenerator>(new omega_decay_gen(mv, mdm, kappa, alD));
				}
				else if(prodchoice=="rho_decay"){
					PartDist->set_mass(mrho);
					DMGen = std::shared_ptr<DMGenerator>(new rho_decay_gen(mv, mdm, kappa, alD));
				}
				else if(prodchoice=="phi_decay"){
					PartDist->set_mass(mphi);
					DMGen = std::shared_ptr<DMGenerator>(new phi_decay_gen(mv, mdm, kappa, alD));
				}
				else if(prodchoice=="omega_decay_baryonic"){
					PartDist->set_mass(momega);
					DMGen = std::shared_ptr<DMGenerator>(new omega_decay_gen_baryonic(mv, mdm, kappa, alD));
				}
				else if(prodchoice=="phi_decay_baryonic"){
					PartDist->set_mass(mphi);
					DMGen = std::shared_ptr<DMGenerator>(new phi_decay_gen_baryonic(mv, mdm, kappa, alD));
				}
				Vnum = DMGen->BranchingRatio()*num_pi0*(proditer->Meson_Per_Pi0());
			}
			else if(prodchoice=="parton_production_baryonic"||prodchoice=="parton_production"){
				PartDist->set_mass(mv);
				DMGen = std::shared_ptr<DMGenerator>(new parton_V_gen(mv, mdm, kappa, alD, prodchoice));
			}
			else if(prodchoice=="V_decay"||prodchoice=="Brem_V"||
	                prodchoice=="V_decay_baryonic"){
				PartDist->set_mass(mv);
				DMGen = std::shared_ptr<DMGenerator>(new V_decay_gen(mv,mdm,kappa,alD,proddist));
				Vnum *= DMGen->BranchingRatio();
			}
	        else if(prodchoice=="piminus_capture"){
	            PartDist = std::shared_ptr<Distribution>(new DoNothingDist());
	            if(proddist!="default"){
	                cerr << "No production distribution supported for piminus_capture";
	            }
				DMGen = std::shared_ptr<DMGenerator>(new piminus_capture_gen(mv,mdm,kappa,alD));
	            Vnum = DMGen->BranchingRatio()*num_pi0*(proditer->Meson_Per_Pi0());
	        }
			else{
				cerr << "Invalid Production Channel Selection: " << prodchoice  << "\n";
				return -1;
			}

			std::shared_ptr<list<production_distribution> > distmodlist = proditer-> Get_Dist_Mods_List();
			for(list<production_distribution>::iterator distiter = distmodlist->begin(); distiter!=distmodlist->end();distiter++){
				list<std::shared_ptr<Distribution> > distlist;
				if(distiter->name()=="position_offset"){
					std::shared_ptr<Distribution> tmpdist (new Position_Offset(distiter->get_offset(0),distiter->get_offset(1),distiter->get_offset(2),distiter->get_offset(3)));
					PartDist->Add_Dist(tmpdist);
				}
			}
			if(outmode=="particle_list"){
				Particle part(0);
				std::ofstream parstream(proditer->Part_List_File(),std::ios::out);
				cout << "--------------------" << endl;
				cout << "Run parameters:" << endl;
				cout << "--------------------" << endl;	
				cout << "Number of events to be generated = " << samplesize  << endl;
				cout << "Production Distribution = " << proddist << endl;
				cout << "Writing to " << proditer->Part_List_File() << endl;
				for(int num = 0; num < samplesize; num++){
					PartDist->Sample_Particle(part);
					parstream << part.px << " " << part.py << " " << part.pz << " " << part.E << " " << part.end_coords[0] << " " << part.end_coords[1] << " " << part.end_coords[2] << " " << part.end_coords[3] << endl;
				}
				parstream.close();
				return 0;
			}
			proddist_vec.push_back(proddist);
			DMGen_list.push_back(DMGen);
			PartDist_list.push_back(PartDist);
			Vnum_list.push_back(Vnum);
			Vnumtot+=Vnum;
		}//End of Production distribution loop. 
		if(Vnumtot==0){
			cout << "No DM production expected. Terminating run.\n";
			return 0;
	}
	//return 0;

	cout << "Preparing signal channel " << sigchoice << endl;
	if(sigchoice=="NCE_electron"){
		SigGen = std::unique_ptr<Scatter>(new Electron_Scatter(mdm, mv, alD, kappa,max_scatter_energy,min_scatter_energy));
	}
	else if(sigchoice=="NCE_nucleon"){
        SigGen = std::unique_ptr<Scatter>(new Nucleon_Scatter(mdm+EDMRES/100.0,max_dm_energy,EDMRES,mdm,mv,alD,kappa,max_scatter_energy,min_scatter_energy,"Kinetic_V",par->Coherent(),det));	
	}
    //Eventually these won't have to be separate channels, they'll just check the model!
	else if(sigchoice=="NCE_nucleon_baryonic"){
		SigGen = std::unique_ptr<Scatter>(new Nucleon_Scatter(mdm+EDMRES/100.0,max_dm_energy,EDMRES,mdm,mv,alD,kappa,max_scatter_energy,min_scatter_energy,"Baryonic_V",par->Coherent(),det));
	}
	else if(sigchoice=="Pion_Inelastic"||sigchoice=="Pion_Inelastic_Charged"||sigchoice=="Inelastic_Delta_to_Gamma"){
		//I might need some checking for allowed energies.
        if(sigchoice=="Pion_Inelastic"){
		    SigGen = std::unique_ptr<Scatter>(new Pion_Inelastic(mdm+EDMRES/100.0,max_dm_energy,EDMRES,mdm,mv,alD,kappa,max_scatter_energy,min_scatter_energy));
        }
        else if(sigchoice=="Inelastic_Delta_to_Gamma"){
		    SigGen = std::unique_ptr<Scatter>(new Pion_Inelastic(mdm+EDMRES/100.0,max_dm_energy,EDMRES,mdm,mv,alD,kappa,max_scatter_energy,min_scatter_energy,1));
        }
        //This covers both charged and neutral pion production
        else if(sigchoice=="Pion_Inelastic_Charged"){
		    SigGen = std::unique_ptr<Scatter>(new Pion_Inelastic(mdm+EDMRES/100.0,max_dm_energy,EDMRES,mdm,mv,alD,kappa,max_scatter_energy,min_scatter_energy,2));
        }
	}
	else if(sigchoice=="Inelastic_Nucleon_Scattering_Baryonic" || sigchoice=="Inelastic_Nucleon_Scattering"){
		if(par->Scatter_Dist_Filename()==""){
			cerr << "No scatter_dist_filename provided for " << sigchoice << endl;
			return -1;
		}
		if(outmode=="comprehensive"){
			cerr << sigchoice << " does not support comprehensive output. Use summary mode instead.\n";
			return -1;
		}
		SigGen = std::unique_ptr<Scatter>(new Inelastic_Nucleon_Scatter(mdm,mv,alD,kappa,sigchoice,par->Scatter_Dist_Filename()));
		min_scatter_energy=0;
		max_scatter_energy=1e9;
	}
    //I'm going to have to rewrite this at some point.
    //Going to shift more of the setup into the model class.
    //Probably also going to make a model class superclass. Pass names of particles/channels to it,
    //and it supplies details about those particles/channels.
    else if(sigchoice=="Signal_Decay"){
        double lifetime;
        vector<double> Branching_Ratios;
        vector<vector<Particle> > Final_States;
        if(par->Model_Name()=="Axion_Dark_Photon"){
            lifetime=adp.Lifetime();
            adp.Branching_Ratios(Branching_Ratios);
            adp.Final_States(Final_States);
            sig_part_name = dark_axion_signal_string;
            SigGen = std::unique_ptr<Scatter>(new SignalDecay(lifetime, Branching_Ratios, Final_States));
        }
        //More temporary stuff ugh
        else if(par->Model_Name()=="Dark_Photon"){
            cout << "Setting up signal decay!"  << endl;
            
            Particle electron(MASS_ELECTRON);
            electron.name = "Electron";
            
            Particle DM(mdm);
            DM.name="DM";

            Particle muon(MASS_MUON);
            muon.name = "Muon";

            Particle hadronic(0);
            hadronic.name = "Hadronic Stuff";
            
            //PLaceholder
            sig_part_name = dark_axion_signal_string;
            
            double GV = Gamma_V(mv,mdm,kappa,alD);
            lifetime=hbar/GV;
            cout << "Width: " << GV << " Lifetime " << lifetime << endl;
            double br = 0;
            if((br=Gamma_V_to_leptons(mv,kappa,MASS_ELECTRON)/GV)>0){
                Branching_Ratios.push_back(br);
                cout << "BR(V->e e) = " << br << endl;
                vector<Particle> vec;
                vec.push_back(electron);
                vec.push_back(electron);
                Final_States.push_back(vec);
            }
            
            if((br=Gamma_V_to_leptons(mv,kappa,MASS_MUON)/GV)>0){
                Branching_Ratios.push_back(br);
                cout << "BR(V->mu mu) = " << br << endl;
                vector<Particle> vec;
                vec.push_back(muon);
                vec.push_back(muon);
                Final_States.push_back(vec);
            }

            if((br=Gamma_V_to_hadrons(mv,kappa)/GV)>0){
                Branching_Ratios.push_back(br);
                cout << "BR(V->hadronic) = " << br << endl;
                vector<Particle> vec;
                vec.push_back(hadronic);
                Final_States.push_back(vec);
            }
            
            if((br=GammaV_to_dm_dm(mv,mdm,kappa,alD)/GV)>0){
                Branching_Ratios.push_back(br);
                cout << "BR(V->DM DM) = " << br << endl;
                vector<Particle> vec;
                vec.push_back(DM);
                vec.push_back(DM);
                Final_States.push_back(vec);
            }

            SigGen = std::unique_ptr<Scatter>(new SignalDecay(lifetime, Branching_Ratios, Final_States));
    	}
        else{
            cerr << "No model declared for Signal_Decay.\n";
            return -1;
        }

    }
    else{
		cerr << "Invalid Channel Selection: " << sigchoice << endl;
		return -1;
	}
	    
	}
	

	SigGen->set_angle_limits(min_angle, max_angle);
	SigGen->set_energy_limits(min_scatter_energy, max_scatter_energy);
	
    //Begin Run
	cout << "--------------------" << endl;	
    cout << "Run parameters:" << endl;	
    cout << "--------------------" << endl;	
    cout << "Number of events to be generated = " << samplesize  << endl;	
    cout << "Dark photon mass = " << mv << " GeV" << endl;	
    cout << "Dark matter mass = " << mdm << " GeV" << endl;	
    cout << "alphaD = " << alD << endl;	
    cout << "kappa = " << kappa << endl;
	for(int i = 0; i<chan_count; i++){
        cout << "Production-Channel " << i+1 << " = " << DMGen_list[i]->Channel_Name();
		if(DMGen_list[i]->query_off_shell()){
			cout << " in Off-Shell mode.\n";
		} 
		else
			cout << " in On-Shell mode.\n";
		cout << "Production Distribution " << i+1 << " = " << proddist_vec[i] << endl;
		cout << "Branching Ratio "  << i+1 << " = " << DMGen_list[i]->BranchingRatio() << endl;
	   	cout <<	"V's produced in channel " << i+1 << " = " << Vnum_list[i] << endl;
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
        return 1;
    }

	////////////////////
	//START OF BURN-IN//
	////////////////////
        
/*
    //TESTING DETECTOR
    cout << "TESTING DETECTOR" << endl;
    Particle test_0(0.01);
    test_0.ThreeMomentum(0.0,0.0,1);
    Particle test_1(0.01);
    test_1.ThreeMomentum(0.04,0.04,1);
    Particle test_2(0.01);
    test_2.ThreeMomentum(0.07,0.07,1);
    Particle test_3(0.01);
    test_3.ThreeMomentum(0.1,0.04,1);

    det_int(test_0);
    cout << "TEST_1\n";
    det_int(test_1);
    cout << "TEST 2\n";
    det_int(test_2);
    cout << "TEST_3\n";
    det_int(test_3);
    //END DETECTOR TEST
*/
    double BURN_MAX = par->Burn_In();
	if(BURN_MAX<0){
        cout << "burn_max < 0 specified. Assuming default value burn_max=1000\n";
        BURN_MAX=1000;
    }
   
    if(outmode=="dm_detector_distribution"){
        BURN_MAX = 0;
        cout << "Detector_Mode selected.\nSkipping Burn-In.\n";
    }

    

    double BURN_OVERRIDE = par->Burn_Timeout();
	for(int i=0;BURN_MAX!=0&&i<chan_count; i++){
        int nburn = 0;
		if(Vnum_list[i]==0){
			cout << "Skipping Channel " << i+1 << ", no events expected.\n";
			continue;
		}
		else{
			cout << "Begin Channel " << i+1 << " Burn-In" << endl;
		}
		for(int burnattempt=0; (nburn < BURN_MAX)&&(burnattempt<BURN_MAX*BURN_OVERRIDE); burnattempt++){
            list<Particle> vecburn;
			Particle dist_part (0);
			PartDist_list[i]->Sample_Particle(dist_part);
			//cout << "det_int " << i << " = " << det_int(dist_part) << endl;
			if(DMGen_list[i]->GenDM(vecburn, det_int, dist_part)){
				for(list<Particle>::iterator burniter = vecburn.begin(); burniter != vecburn.end(); burniter++){
                    if(burniter->name==sig_part_name){
                        SigGen->probscatter(det, *burniter);
						nburn++;
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

    int trials = 0;
    vector<long> trials_list(chan_count,0);
	vector<int> scat_list(chan_count,0);
    int nevent=0;
    vector<int> NDM_list(chan_count,0);
    vector<double> timing_efficiency(chan_count,0.0);
	
	//int escat=0;lso thought it could have used Leatherhead, or some of the other mutanimals. He probably would have fit right in with the Scale tail clan.
    bool scatter_switch;
    int trials_max = par->Max_Trials();
	//if(SigGen->get_pMax()*Vnumtot<=1){
    if(SigGen->get_pMax()<=0){
        cout << "pMax less than tolerance limit, skipping remainder of run\n";
    }
    else{
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
            if(DMGen_list[i]->GenDM(vec, det_int, dist_part)){
                //Yes, this list is named vec.  
                for(list<Particle>::iterator iter = vec.begin(); iter != vec.end();iter++){
                //The way this is structured means I can't do my usual repeat thing to boost stats. 
                    if(iter->name.compare(sig_part_name)==0){
                        //iter->report(logging);
                        NDM_list[i]++;
                        if(outmode=="dm_detector_distribution"){
                            *comprehensive_out << DMGen_list[i]->Channel_Name() << " " << det->Ldet(*iter) << " ";
                            iter->report(*comprehensive_out);
                            scatter_switch=true;
                            continue;;
                        }
                        //may need to replace this with a list<Particle> later
                        //cout << SigGen->get_pMax() << endl;
                        if(SigGen->probscatter(det, vec, iter)){
                            //cout << "Scatter?\n"; 
                            scat_list[i]++;
                            if(timing_cut>0){
                                timing_efficiency[i]+=t_delay_fraction(timing_cut,sqrt(pow(iter->end_coords[0],2)+pow(iter->end_coords[1],2)+pow(iter->end_coords[2],2)),iter->Speed());
                            }
                            else{
                                timing_efficiency[i]+=1;
                            }
                            scatter_switch = true;	
                        }
                        else{
                            iter = vec.erase(iter);
                            iter--;
                        }
                    }    
                }
                
            }
            if(scatter_switch&&outmode=="comprehensive"){
                *comprehensive_out << "event " << ++nevent << endl;
                Record_Particles(*comprehensive_out, vec);
                *comprehensive_out << "endevent " << nevent << endl << endl;    
            }
            else if(scatter_switch)
                ++nevent;
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

	for(int i=0; i<chan_count; i++){
		if(scat_list[i]==0)
			signal_list[i]=0;
		else
			signal_list[i] = (double)scat_list[i]/(double)trials*Vnumtot*SigGen->get_pMax()/repeat*par->Efficiency()*timing_efficiency[i]/scat_list[i];
		scattot+=scat_list[i];
  		cout << DMGen_list[i]->Channel_Name() << ": " << (double)scat_list[i]/(double)trials_list[i]*Vnum_list[i]*SigGen->get_pMax()/repeat*par->Efficiency()*timing_efficiency[i]/scat_list[i];
		cout << " Timing_Efficiency: " << timing_efficiency[i]/scat_list[i] << " ";
		cout << "Events: " <<scat_list[i] << " Trials: " << trials_list[i] << " V_num: " << Vnum_list[i] << " pMax: " << SigGen->get_pMax() << " repeat: " << repeat << " efficiency: "  << par->Efficiency() << endl;;
		if(outmode=="summary"||outmode=="dm_detector_distribution"||
                outmode=="comprehensive"){
            *summary_out << DMGen_list[i]->Channel_Name() << " " << mv  <<  " "  << mdm << " " << signal_list[i] << " " << kappa << " " << alD << " " << sigchoice << " " << POT << " " << par->Efficiency() << " " << samplesize << " " << Vnum_list[i] << " " << Vnumtot << endl;
        }
        NDM+=NDM_list[i]; 
		signal+=signal_list[i];
 	}
	if(outmode=="summary"||outmode=="dm_detector_distribution"||
            outmode=="comprehensive"){
		*summary_out << "Total " << mv  <<  " "  << mdm << " " << signal << " " << kappa << " " << alD << " " << sigchoice << " " << POT << " " << par->Efficiency() << " " << samplesize << " " << endl;
        if(par->Model_Name()=="Axion_Dark_Photon"){
            adp.Report(*summary_out, signal);
        }
    }

//    if(outmode=="summary"||outmode=="comprehensive"){
//        *summary_out << "Total " << mv  <<  " "  << mdm << " " << signal << " " << kappa << " " << alD << " " << sigchoice << " " << POT << " " << par->Efficiency() << " " << samplesize << " " << endl;
//    }
    else if(outmode=="dm_detector_distribution"){
        *summary_out << "Total " << mv  <<  " "  << mdm << " " << trials << " " << kappa << " " << alD << " " << sigchoice << " " << POT << " " << Vnumtot << " " << samplesize << " " << (double)NDM/(2*trials) << " " << endl;
    }

    summary_out->close();
	
  	if(outmode=="dm_detector_distribution"||outmode=="comprehensive")
		comprehensive_out->close();

    cout << "Number of trials = " << trials << endl;	
	cout << "Number of candidates intersecting detector = " << NDM << endl;	
	cout << "Number of " << sigchoice <<  " = " << scattot << endl;
	for(int i=0; i<chan_count;i++)
		cout << "Number of scatterings from channel " << i+1 << " " << DMGen_list[i]->Channel_Name() << " = " << scat_list[i] << " in " << trials_list[i] << " trials." << endl;  
	cout << "Number of sample events generated = " << nevent << endl;
	cout << "Acceptance = " << (double)NDM/(2*trials) << endl;	
	cout << "Maximum event probability = "  << SigGen->get_pMax() << endl;
	cout << "Average event probability = "  << (double)scattot/NDM*SigGen->get_pMax() << endl;
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
