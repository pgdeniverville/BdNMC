#include "branchingratios.h"
#include <iostream>
#include <functional>

using std::cout;
using std::endl;
using namespace std::placeholders;

int main(){
	cout << d2N_proton_brem_to_V(70, 0.2, 1, 0.2, 0.5) << endl;
}
