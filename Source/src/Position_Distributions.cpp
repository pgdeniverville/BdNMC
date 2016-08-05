#include "Position_Distributions.h"

void Position_Offset::sample_particle(Particle &part){
//	part.Set_Position(part.origin_coords[0]+X,part.origin_coords[1]+Y,part.origin_coords[2]+Z);
	part.Set_Origin(part.origin_coords[0]+X,part.origin_coords[1]+Y,part.origin_coords[2]+Z);
	part.Set_Creation_Time(part.origin_coords[3]+T);
	part.Set_Time(part.origin_coords[3]+T);
}
