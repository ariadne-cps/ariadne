#include <ariadne.h>

#include <iostream>
#include <string>
#include <sstream>

using namespace Ariadne;
using namespace std;

int main() {
    interval<double> ivl1(1.1,2.2);
    interval<double> ivl2;
    interval<double> ivl3(2.1,3.2);

    string input("[1.1,2.2] ");
    stringstream iss(input);
	
    iss >> ivl2;
    if(!(ivl1.lower()==ivl2.lower() && ivl1.upper()==ivl2.upper())) {
	cout << "test_interval: FAILED construction from stream\n";
	return 1;
    }

    cout << "test_interval: PASS\n";
    
    return 0;
}
