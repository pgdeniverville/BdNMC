#include "record.h"
#include <iomanip>

using std::list;        using std::cout;
using std::ostream;     using std::endl;
using std::setw;

const int MAXWIDTH = 15;

ostream& Record_Particles(ostream& outstr, list<Particle> &partlist){
    for(list<Particle>::iterator it = partlist.begin(); it != partlist.end(); it++){
        outstr << " " << it->name << " " << it->px << " " << it->py << " " << it->pz << " " << it->E;
        //if(it->EVENT_SET){
            outstr << " " << it->origin_coords[0] << " " << it->origin_coords[1] << " " << it->origin_coords[2] << " " << it->origin_coords[3];
        //}
        outstr << endl;
    }
    outstr << endl;

    return outstr;
}
