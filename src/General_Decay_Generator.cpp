#include "DMgenerator.h"
#include "Random.h"

bool General_Decay_Generator::GenDM(std::list<Particle>& vec, std::function<double(Particle&)> det_int, Particle& part){
	double u = Random::Flat(0,total_width);

	for(std::vector<std::shared_ptr<DMGenerator>>::size_type i = 0; i < decay_channels.size(); i++){
		if(u < partial_widths[i]){
			return decay_channels[i]->GenDM(vec, det_int, part);
		}
		else
			u-=partial_widths[i];
	}

	return false;
}