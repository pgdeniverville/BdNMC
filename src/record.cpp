#include "record.h"

#include <iomanip>

using std::list;        using std::cout;
using std::ostream;     using std::endl;
using std::setw;

const int MAXWIDTH = 15;

ostream& Record_Particles(ostream& outstr, list<Particle> &partlist){
    for(list<Particle>::iterator it = partlist.begin(); it != partlist.end(); it++){
        it->report(outstr);
    }

    return outstr;
}
