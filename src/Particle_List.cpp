
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

Particle_List::Particle_List(const string& infile, bool set_pos){
	POS = set_pos;
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
//Note that batch_size is not yet implemented for some reason.
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
			if(!POS)
				cerr << "Format required is p_x p_y p_z" << endl;
			else
				cerr << "Format required is p_x p_y p_z x y z t" << endl;
			cerr << "Line is improperly formatted: " << hold << endl;
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
		double x,y,z,t;
		x=0;y=0;z=0;t=0;
		if(POS){
			if(sublist.size()>=7){
				x = stod(sublist[4]);
				y = stod(sublist[5]);
				z = stod(sublist[6]);
				if(sublist.size()>=8)
					t = stod(sublist[7]);
			}
		}
		partlist.push_back(_part(px,py,pz,E,x,y,z,t));
	}
	catch(std::exception e){
		return -1;
	}
	return 0;
}

void Particle_List::sample_particle(Particle &part){
	if(iter==partlist.end()){
		iter = partlist.begin();
	}
	part.ThreeMomentum(iter->px,iter->py,iter->pz);
	if(POS){
		part.Set_Origin(iter->x,iter->y,iter->z);
		part.Set_Creation_Time(iter->t);
 		//part.Set_Position(iter->x,iter->y,iter->z);
		part.Set_Time(iter->t);
	}
	iter++;
}

