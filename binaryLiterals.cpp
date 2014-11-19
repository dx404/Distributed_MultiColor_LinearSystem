/**
 *  g++     -std=c++1y binaryLiterals.cpp 
 *  clang++ -std=c++1y binaryLiterals.cpp 
 *
 */
#include <iostream>
#include <vector>
using namespace std;

int main(int argc, char *argv[]) {
	vector<int> num = {
		0b11'111, 0b11'111, // binary, digit separators
		037,				// octary 
		31,					// decimal
		0x1f, 0x1f			// hexadecimal 	
	};
	for (auto x : num) {
		cout << x << endl;  // decimal output
	}
	return 0;
}