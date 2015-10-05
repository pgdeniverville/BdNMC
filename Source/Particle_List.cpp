
#include "Distribution.h"
#include "Particle_List.h"
#include "Kinematics.h"
#include <vector>
#include <exception>
#include <cmath>

using std::string;		using std::cerr;
using std::vector;
using std::cout;
using std::endl;
using std::stod;

Particle_List::Particle_List(const string& infile){
	//filename = infile;
	part_count = 0;
	instream.open(infile, std::ifstream::in);
	if(instream.is_open()){
		if(load_particle_batch()==0){
			std::cerr << "File empty: " << infile << endl;
			throw -1;	
		}
	}
	else{
		std::cerr << "Could not open " << infile << endl;
		throw -1;
	}
	iter=partlist.begin();
}
//Loads particles into the list until the entire file is read, or the list is of length batch_size.
//If batch_size is negative, the entire file is read.
int Particle_List::load_particle_batch(){
	int i=0;
	int error_state;
	/*if(instream.eof()){
		cout << "Reload\n"; 
		instream.seekg(0);
	}*/
	string hold;
	while(std::getline(instream, hold)){
		if((error_state=parse_line(hold))==0){
			i++;
			continue;
		}
		else if(error_state==1){
			continue;
		}
		else if(error_state==-1){
			cerr << "Format required is p_x p_y p_z" << endl;
			cerr << "Line is improperly formatter: " << hold << endl;
		}
		else{
			cerr << "Unknown error code: " << error_state << "\n";
			cerr << "Line is improperly formatted: " << hold << "\n";
		}
	}
	return i;
}

int Particle_List::parse_line(string& line){
	if(line.length()==0||line[0]=='#')
		return 1;//Line can be skipped.
	std::size_t space_pos=0;
	std::size_t next;
	vector<string> sublist;
	while((next=line.find(' ',space_pos))!=std::string::npos){
		sublist.push_back(line.substr(space_pos, next-space_pos));
		space_pos=line.find(' ',space_pos)+1;
	}
	sublist.push_back(line.substr(space_pos));
	double px, py, pz, E;
	try{
		px = stod(sublist[0]);
		py = stod(sublist[1]);
		pz = stod(sublist[2]);
		E = stod(sublist[3]);
		//partlist.push_back(std::shared_ptr<_3momentum>(new _3momentum(px,py,pz,E)));
		partlist.push_back(_3momentum(px,py,pz,E));
	}
	catch(std::exception e){
		return -1;
	}
	return 0;
}

void Particle_List::sample_momentum(double &pmom, double &theta, double &phi){
	Get_Particle(pmom, theta, phi);
}

void Particle_List::Get_Particle(double &pmom, double &theta, double &phi2){
	if(iter==partlist.end()){
		iter = partlist.begin();
	}
	//theta=calculatetheta((*iter)->px, (*iter)->py, (*iter)->pz, (*iter)->E);
	//phi2=phi((*iter)->px, (*iter)->py, (*iter)->pz, (*iter)->E);
	//pmom = sqrt(pow((*iter)->px,2)+pow((*iter)->py,2)+pow((*iter)->pz,2));	
	theta=calculatetheta(iter->px, iter->py, iter->pz, iter->E);
	phi2=phi(iter->px, iter->py, iter->pz, iter->E);
	pmom = sqrt(pow(iter->px,2)+pow(iter->py,2)+pow(iter->pz,2));	
	iter++;
}
