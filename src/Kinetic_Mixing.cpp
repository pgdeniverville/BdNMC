#include "Model.h"

bool Kinetic_Mixing::Set_Model_Parameters(Parameter& par){
    if(par.Query_Map(g_chi_key,gchi)&&(par.Query_Map(dark_matter_mass_key,mchi))&&par.Query_Map(dark_mediator_mass_key, ma)&&par.Query_Map(g_quark_key,gq)&&par.Query_Map(g_electron_key,gae)){
        return true;
    }
    return false;
}