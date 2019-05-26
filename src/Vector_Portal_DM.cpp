/* 
*This file will hold all the standard Dark Photon stuff.
*I will be migrating everything from main to here eventually.
*/

#include "Model.h"
#include "branchingratios.h"
#include "Scatter.h"
#include <cmath>
#include "decay.h"
#include "DMgenerator.h"
#include "Model.h"
#include "Kinematics.h"
#include "SignalDecay.h"

using std::string;
using std::vector;
using std::cout;
using std::endl;
using std::bind;
using std::function;

using namespace std::placeholders;