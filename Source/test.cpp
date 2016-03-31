#include "constants.h"
#include "Kinematics.h"

#include <iostream>
#include <functional>


using std::cout;
using std::endl;
using namespace std::placeholders;

int main(){
	cout << Enmin(1,0.01,0.938) << endl;
	cout << Enmax(1,0.01,0.938) << endl;
}
