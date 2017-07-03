#include "Parameter.h"
#include <vector>
#include <map>
#include <stdlib.h>
#include <exception>
#include <time.h>
#include <locale>
#include "constants.h"

#include "parameter_keys.h"

using std::exception;
using std::string;
using std::map;
using std::cout; using std::cerr;
using std::endl;
using std::vector;
using std::list;


//All defaults should be stored here at some point.

int split_string(string &hold, string &key, string &val){
    if(hold.length()==0||hold[0]=='#')
        return 0;
    
    std::size_t space_pos;
    if((space_pos = hold.find(" "))==std::string::npos)
        return -1;
    
    key = hold.substr(0,space_pos);
    val = hold.substr(space_pos+1, hold.find(" ", space_pos+1));
    return 1;
}

void parse_parameter_file(const string &filename, map<string, string> &parammap){
	std::ifstream instream(filename);
	if(instream.is_open()){
		string hold, key, val;
        int line_num=0;
        int error_state = 0;
		//Must be used if you are using different parsing modes (like materials!)
        //int last_pos = instream.tellg();
        while(getline(instream, hold)){
            line_num++;
			if((error_state=split_string(hold, key, val))==0);
            else if(error_state==-1){
                cerr << line_num << ": Improperly formatted input: ";
                cerr << hold << "\n";
            }
            else if(error_state!=1){
                cerr << line_num << ": Unknown error code: " << error_state << "\n";
            }
			else{
                parammap[key] = val;
            }
            //last_pos=instream.tellg(); 
        }

		instream.close();
	}
	else
		std::cerr << "Could not open " << filename << endl;	

}//end of parse_parameter_file



bool production_channel::query_dist_param(const string &key, double &var){
    try{
        if(dist_param_map.count(key)==1){
            var = stod(dist_param_map[key]);
			return true;
		}
    }
    catch(exception& e){
		cerr << "Exception: " << e.what() << endl;
    	cerr << "Invalid model parameter: " << key << endl;
    }
	return false;
}

production_distribution::production_distribution(){
	dist_mod = "";
}
double production_distribution::get_offset(int i){
	if(i>3||i<0){
		std::cerr << "Out of domain for get_offset\n";
	}	
	return offset_coords[i];
}

production_distribution parse_distribution_mod(std::ifstream &instream, string &hold, string &key, string &val, int &line_num, int &error_state, int &last_pos){
	production_distribution tmpdist;
	tmpdist.dist_mod=val;
	while(getline(instream, hold)){
		line_num++;
		//cout << "Position: " << instream.tellg() << endl;
		//cout << "Distribution parser " << line_num << ": " << hold << endl;
        if((error_state=split_string(hold, key, val))==0){
            last_pos = instream.tellg();
			continue;
		}
        else if(error_state==-1){
            cerr << line_num << ": Improperly formatted input: ";
            cerr << hold << "\n";
            cerr << "Terminating distribution parsing.\n";
            break;
        }
        else if(error_state!=1){
            cerr << "Unknown error code: " << error_state << "\n";
            cerr << "Terminating distribution parsing.\n";
            break;
        }
        try{
			if(key == "x-offset"){
				tmpdist.offset_coords[0]=stoi(val);
			}
			else if(key == "y-offset"){
				tmpdist.offset_coords[1]=stoi(val);
			}
			else if(key == "z-offset"){
				tmpdist.offset_coords[2]=stoi(val);
			}
			else if(key == "t-offset"){
				tmpdist.offset_coords[3]=stoi(val);
			}
			else{
				line_num--;
				break;
			}
        }
        catch(exception& e){
            cerr << "Exception: " << e.what() << " encountered while parsing line " << line_num << ":" << endl;
            cerr << hold << endl;
        }
		last_pos = instream.tellg();
		//cout << "last_pos = " <<  last_pos << endl;
    }
	//cout << "end of dist parser position = " << last_pos << endl;
    return tmpdist;
}

production_channel::production_channel(){
	prod_chan="";
   	prod_dist="default"; 
	particle_list_file="";
   	sanfordwang_file="";
   	parton_V_n_file=""; 
	parton_V_p_file=""; 
	meson_per_pi0=-1; 
	PTMAX=-1; 
	PTMIN=0;
    ZMIN=-1;
   	ZMAX=-1;
	particle_list_position=false;
	dist_mods = std::unique_ptr<list<production_distribution> > (new list<production_distribution>);
}

string lowercase(const string &str){
	std::locale loc;
	string str2(str);
	for(std::string::size_type i=0; i<str.length(); ++i)
		str2[i] = tolower(str[i], loc);
	return str2; 
}

production_channel parse_production_channel(std::ifstream &instream, string &hold, string &key, string &val, int &line_num, int &error_state, int &last_pos){
	production_channel tmpprod;
	tmpprod.prod_chan=val;
	while(getline(instream, hold)){
		line_num++;
		//cout << "Position: " << instream.tellg() << endl;
		//cout << line_num << ": " << hold << endl; 
        if((error_state=split_string(hold, key, val))==0){
            last_pos = instream.tellg();
			continue;
		}
        else if(error_state==-1){
            cerr << line_num << ": Improperly formatted input: ";
            cerr << hold << "\n";
            cerr << "Terminating production channel parsing.\n";
            break;
        }
        else if(error_state!=1){
            cerr << "Unknown error code: " << error_state << "\n";
            cerr << "Terminating production channel parsing.\n";
            break;
        }
        try{
			if(key==production_distribution_key)
				tmpprod.prod_dist=val;
			else if(key==particle_list_file_key)
				tmpprod.particle_list_file=val;
			else if(key==meson_per_pi0_key)
				tmpprod.meson_per_pi0=stod(val);
			else if(key==parton_V_neutron_file_key)
				tmpprod.parton_V_n_file=val;
			else if(key==parton_V_proton_file_key)
				tmpprod.parton_V_p_file=val;
			else if(key==zmax_key)
				tmpprod.ZMAX=stod(val);
			else if(key==zmin_key)
				tmpprod.ZMIN=stod(val);
			else if(key==ptmax_key)
				tmpprod.PTMAX=stod(val);
            else if(key==ptmin_key)
                tmpprod.PTMIN=stod(val);
			else if(key==part_list_pos_key){
				if(lowercase(val)=="true")
					tmpprod.particle_list_position=true;
			}
			else if(key == sanford_wang_key || key == distribution_parameter_key){
				parse_parameter_file(val, tmpprod.dist_param_map);
			}
			else if(key == dist_mod_key){
				//cout << hold << endl;
				(tmpprod.dist_mods)->push_back(parse_distribution_mod(instream, hold, key, val, line_num, error_state, last_pos));
				instream.seekg(last_pos);
                continue;

			//	cout << hold << endl;
			}
			else{
				line_num--;
				break;
			}
        }
        catch(exception& e){
            cerr << "Exception: " << e.what() << " encountered while parsing line " << line_num << ":" << endl;
            cerr << hold << endl;
        }
		last_pos = instream.tellg();
		//cout << "last_pos = " <<  last_pos << endl;
    }
	//cout << "end of prod parser position = " << last_pos << endl;
    return tmpprod;
}

struct _material{
    double n_num, p_num, e_num, n_density, mass;
    string name;
};

/*
struct _material{
    double n_num, p_num, e_num, n_density, mass;
    string name;
};
*/
_material parse_material(std::ifstream &instream, string &hold, string &key, string &val, int &line_num, int &error_state, int &last_pos){
    _material tmpmat;
    tmpmat.n_num=0; tmpmat.p_num=0; tmpmat.e_num=0; tmpmat.n_density=0;
    tmpmat.name = val;
    while(getline(instream, hold)){
        line_num++;
        if((error_state=split_string(hold, key, val))==0){
        	last_pos=instream.tellg();
			continue;
		}
        else if(error_state==-1){
            cerr << line_num << ": Improperly formatted input: ";
            cerr << hold << "\n";
            cerr << "Terminating material parsing.\n";
            break;
        }
        else if(error_state!=1){
            cerr << "Unknown error code: " << error_state << "\n";
            cerr << "Terminating material parsing.\n";
            break;
        }
        try{
			if(key=="number_density")
				tmpmat.n_density=stod(val);
			else if(key =="proton_number")
				tmpmat.p_num = stod(val);
			else if(key =="neutron_number")
				tmpmat.n_num=stod(val);
			else if(key =="electron_number")
				tmpmat.e_num=stod(val);
			else if(key =="mass")
				tmpmat.mass=stod(val);
			else{
				break;
			}
        }
        catch(exception& e){
            cerr << "Exception: " << e.what() << " encountered while parsing line " << line_num << ":" << endl;
            cerr << hold << endl;
        }
        last_pos=instream.tellg();
    }
    return tmpmat;
}
/*
 *  Builds a parameter object from the parameter stream supplied.
 *
 */
Parameter::Parameter(std::ifstream &instream){
	prodlist = std::shared_ptr<list<production_channel> > (new list<production_channel>); 
	if(instream.is_open()){
        string hold, key, val;
        vector<_material> mat_vec;
        int line_num=0;
        int error_state = 0;
        int last_pos = instream.tellg();
        integrity=0;
        while(getline(instream, hold)){
			line_num++;
			//cout << "Position: " << instream.tellg() << endl;
			//cout << "parameter parser line " << line_num << ": " << hold << endl;
            if((error_state=split_string(hold, key, val))==0);
            else if(error_state==-1){
                cerr << line_num << ": Improperly formatted input: ";
                cerr << hold << "\n";
            }
            else if(error_state!=1){
                cerr << "Unknown error code: " << error_state << "\n";
            }
            else if(key == "material"){
                mat_vec.push_back(parse_material(instream, hold, key, val, line_num, error_state, last_pos));
                instream.seekg(last_pos);
                continue;
            }
			else if(key == "production_channel"){
				prodlist->push_back(parse_production_channel(instream, hold, key, val, line_num, error_state, last_pos));
				instream.seekg(last_pos);
				continue;
			}	
            else{
                keymap[key] = val;
            }
            last_pos=instream.tellg(); 
        }
/*        for(std::map<string, string>::iterator it=keymap.begin(); it!=keymap.end(); ++it)
            cout << it->first << " => " << it->second << endl;
        for(std::vector<_material>::iterator it=mat_vec.begin(); it!=mat_vec.end(); ++it)
            cout << it-> name << " " << it->n_density << " " << it->p_num << " " << it->e_num << " " << it->n_num << " " << it->mass << endl;*/

        
        Build_Detector(keymap);
        for(vector<_material>::iterator it = mat_vec.begin(); it != mat_vec.end(); ++it){
            det->add_material(it->n_density, it->p_num, it->n_num, it->e_num,it->mass,it->name);
        }
        Set_String(modelkey,model_name,keymap,"Hidden_Sector_DM");
        Set_Model_Parameters(keymap);    
        Set_Double(POTkey, POT, keymap);
		Set_Double(pi0_ratio_key,pi0ratio,keymap,0.9);
		Set_Integer("samplesize", samplesize, keymap);
        //Production Channel
		if(prodlist->size()==0){
            cerr << "No production channel selected\n";
            integrity=-1;
        }
		//Signal Channel
        if(keymap.count(signalkey)==1){
            sig_channel = keymap[signalkey];
        }
        else{
            cerr << "No signal channel selected\n";
            integrity=-1;
        }
        Set_Bool(coherent_key, coherent, keymap, false); 
		Set_String(output_mode_key, output_mode, keymap, out_mode_def);	
		Set_String(run_key, run_name, keymap, std::to_string(time(NULL)));	
		for(list<production_channel>::iterator iter = prodlist->begin(); iter!=prodlist->end(); iter++)
		{
			if(iter->prod_dist=="particle_list"||output_mode=="particle_list"){
				if(iter->particle_list_file==""){
					integrity=-1;
					cerr << "particle_list_file undefined for production_distribution particle_list.\n"; 
				}
			}
			else if(iter->prod_dist=="parton_V"||iter->prod_dist=="parton_V_baryonic"){
				if(iter->parton_V_n_file==""&&iter->parton_V_p_file==""){
					cerr << "parton_V production distribution files not supplied.\n"; 
					integrity = -1;		
				}
			}
		}
		
		if(prodlist->front().prod_dist==""&&output_mode=="particle_list"){
			integrity = -1;
			cerr << "Must select a production distribution for this run mode\n";
		}
		Set_Double(beam_energy_key, beam_energy, keymap, 8.0);
		Set_Double(max_scatter_energy_key, max_scatter_energy, keymap, 1e9);
		Set_Double(min_scatter_energy_key, min_scatter_energy, keymap, 0.0);
		Set_Double(max_dm_energy_key, max_dm_energy, keymap, beam_energy);
		Set_Double(dm_energy_res_key, edmres, keymap, 0.1);
		Set_Double(angle_lower_limit_key, angle_lower_limit, keymap, 0);
		Set_Double(angle_upper_limit_key, angle_upper_limit, keymap, 2.1*pi);
		Set_Double(timing_key, timing_cut, keymap, 0);
        Set_String(scatter_dist_filename_key, scatter_dist_filename, keymap, "");

		Set_Target_Parameters(keymap);

		Set_Double(p_cross_key, p_cross, keymap, 0);

		//Burn-in details
		Set_Integer(burn_in_key, burn_in, keymap, 1000);
		Set_Integer(burn_timeout_key, burn_timeout, keymap, 20000);
		Set_Integer(repeat_key, repeat, keymap, 1);
		Set_Integer(seed_key, seed, keymap, -1);
		
		Set_Integer(max_trials_key,max_trials,keymap,-1);

		Set_Double(efficiency_key, efficiency, keymap, 1.0);
		Set_Double(meson_per_pi0_key, meson_per_pi0, keymap, -1.0);
		//output file
        if(keymap.count(outkey)==1){
            output_file = keymap[outkey];
        }
        else{
            output_file = "";
        }
		if(keymap.count(sumkey)==1){
            summary_file = keymap[sumkey];
        }
        else{
            summary_file = "";
        }
    }
    else{
        std::cerr << "Parameter file not open.\n";
        throw -1; 
    }
}//end of Parameter::Parameter
//This parses a distribution parameter file


void Parameter::Build_Detector(map<string, string> &keymap){
    if(keymap.count(detkey)==1){
            if(keymap[detkey]=="sphere"){
                for(vector<string>::iterator it=spherearr.begin(); it!=spherearr.end(); ++it){
                    if(keymap.count(*it)!=1){
                        cerr << "Missing detector_sphere parameter: " << *it << endl;
                        integrity = -1;
                    } 
                }
                if(integrity==0){
                    try{
                    det= std::make_shared<detector_sphere>(stod(keymap["x-position"]),stod(keymap["y-position"]),
                            stod(keymap["z-position"]),stod(keymap["radius"]));
                    }
                    catch(exception& e){
                        cerr << "Exception: " << e.what() << " encountered constructing detector object\n";
                        integrity = -1;
                    }
                }
            }
            else if(keymap[detkey]=="cylinder"){
                for(vector<string>::iterator it = cylinderarr.begin(); it!=cylinderarr.end(); ++it){
                    if(keymap.count(*it)!=1){
                        cerr << "Missing detector_cylinder parameter: " << *it << endl;
                        integrity = -1;
                    } 
                }
                if(integrity==0){
                    try{
                    det = std::make_shared<detector_cylinder>(stod(keymap["x-position"]),stod(keymap["y-position"]),
                            stod(keymap["z-position"]),stod(keymap["length"]),stod(keymap["radius"]),
                            stod(keymap["det-theta"]),stod(keymap["det-phi"]));
                    }
                    catch(exception& e){
                        cerr << "Exception: " << e.what() << " encountered constructing detector object\n";
                    }
                }
            }
            else if(keymap[detkey]=="cuboid"){
                for(vector<string>::iterator it = cuboidarr.begin(); it!=cuboidarr.end(); ++it){
                    if(keymap.count(*it)!=1){
                        cerr << "Missing detector_cuboid parameter: " << *it << endl;
                        integrity = -1;
                    } 
                }
                if(integrity==0){
                    try{
                        det = std::make_shared<detector_cuboid>(stod(keymap["x-position"]),stod(keymap["y-position"]),
                            stod(keymap["z-position"]),stod(keymap["width"]),stod(keymap["height"]),stod(keymap["length"]),
                            stod(keymap["det-phi"]),stod(keymap["det-theta"]),stod(keymap["det-psi"]));
                    }
                    catch(exception& e){
                        cerr << "Exception: " << e.what() << " encountered constructing detector object\n";
                    }
                }
            }
            else{
                cerr << "Invalid detector shape declared: " << keymap[detkey] << endl;
                integrity = -1;
            }
             
        }
        else{
            cerr << "No detector shape declared." << endl;
            integrity = -1;
        }
    
}//end of Parameter::Build_Detector

void Parameter::Set_Bool(const string &key, bool &var, map<string, string> &keymap, const bool &def){
    try{
        if(keymap.count(key)==1){
            if(keymap[key]=="true"){
               var = true; 
            }
        }
        else{
            var = def;
        }
    }
    catch(exception& e){
            cerr << "Invalid model parameter: " << key << endl;
            integrity = -1;
    }
}

//These should be replaced with a single overloaded function.
void Parameter::Set_Double(const string &key, double &var, map<string, string> &keymap){
    try{
        if(keymap.count(key)==1){
            var = stod(keymap[key]);
        }
        else{
            cerr << key << " not set." << endl;
            integrity = -1;
        }
    }
    catch(exception& e){
            cerr << "Invalid model parameter: " << key << endl;
            integrity = -1;
    }
}

void Parameter::Set_Double(const string &key, double &var, map<string, string> &keymap, const double &def){
    try{
        if(keymap.count(key)==1){
            var = stod(keymap[key]);
        }
        else{
            var = def;
        }
    }
    catch(exception& e){
            cerr << "Invalid model parameter: " << key << endl;
            integrity = -1;
    }
}

void Parameter::Set_Integer(const string &key, int &var, map<string, string> &keymap){
    try{
        if(keymap.count(key)==1){
            var = stoi(keymap[key]);
        }
        else{
            cerr << key << " not set." << endl;
            integrity = -1;
        }
    }
    catch(exception& e){
            cerr << "Invalid model parameter: " << key << endl;
            integrity = -1;
    }
}

void Parameter::Set_Integer(const string &key, int &var, map<string, string> &keymap, const int &def){
    try{
        if(keymap.count(key)==1){
            var = stoi(keymap[key]);
        }
        else{
            var = def;
        }
    }
    catch(exception& e){
            cerr << "Invalid model parameter: " << key << endl;
            integrity = -1;
    }
}

void Parameter::Set_String(const string &key, string &var, map<string, string> &keymap, const string &def){
	if(keymap.count(key)==1)
		var = keymap[key];
	else
		var = def;
}

const string e_num_target_key = "e_num_target";
const string n_num_target_key = "n_num_target";
const string p_num_target_key = "p_num_target";
const string n_dens_target_key = "n_density_target";
const string target_mass_key = "atomic_mass_target";
const string target_length_key = "target_length";

void Parameter::Set_Target_Parameters(map<string, string> &keymap){
	Set_Double(e_num_target_key, target_e_num, keymap, 0);
	Set_Double(n_num_target_key, target_n_num, keymap, 0);
	Set_Double(p_num_target_key, target_p_num, keymap, 0);
	Set_Double(n_dens_target_key, target_n_density, keymap, 0);
	Set_Double(target_length_key, target_length, keymap, 0);
	Set_Double(target_mass_key, target_mass, keymap, 0);
}

void Parameter::Set_Model_Parameters(map<string, string> &keymap){
    Set_Double(dmkey, mass_dm, keymap);
    Set_Double(Vkey, mass_V, keymap);
    Set_Double(epskey, epsilon, keymap);
    Set_Double(alkey, alpha_D, keymap);
}

bool Parameter::Query_Map(const string key, double& val){
    if(keymap.count(key)!=1){
        val = 0.0;
        return false;
    }
    else{
        val = stod(keymap[key]);
        return true;
    }
}
