#include "record.h"
#include <iomanip>

using std::list;        using std::cout;
using std::ostream;     using std::endl;
using std::setw;

const int MAXWIDTH = 15;

ostream& Record_Particles(ostream& outstr, list<Particle> &partlist){
    for(list<Particle>::iterator it = partlist.begin(); it != partlist.end(); it++){
        outstr << std::setfill(' ') << std::left;
        outstr << setw(MAXWIDTH) << it->name << setw(MAXWIDTH) << it->px << setw(MAXWIDTH) << it->py << setw(MAXWIDTH) << it->pz << setw(MAXWIDTH) << it->E;
        //if(it->EVENT_SET){
            outstr << setw(MAXWIDTH) << it->origin_coords[0] << setw(MAXWIDTH) << it->origin_coords[1] << setw(MAXWIDTH) << it->origin_coords[2] << setw(MAXWIDTH) << it->origin_coords[3];
        //}
        outstr << endl;
    }
    outstr << endl;

    return outstr;
}
