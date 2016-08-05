#include <time.h>
#include <iostream>
#include <ios>
#include <string>
using std::cout; using std::endl;


int main(){
	std::string now = std::to_string(time(NULL));
	cout << now << endl;
}
