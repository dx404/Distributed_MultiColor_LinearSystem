// data generation
/**
 * 
 	g++ -std=c++0x dataGen.cpp -o dataGen.out
 	./dataGen.out 10 10 10 > dataSheet.txt
 */
#include <iostream>
#include <random>
using namespace std;

int main (int argc, char *argv[]) {
	if (argc < 4) {
		cerr << "Please Input the dimensions for i, j, k: " << endl;
		return 97; 
	}
	int i_bound = atoi(argv[1]);
	int j_bound = atoi(argv[2]);
	int k_bound = atoi(argv[3]);
	int totalSize = i_bound * j_bound * k_bound;
	
	cout << i_bound << " " << j_bound << " " << k_bound << " " << endl;
	cout << endl;

	uniform_real_distribution<double> unif(-2, 2);
	random_device rand_dev;
    mt19937 rand_engine(rand_dev());


    // x intialized to zero
    // the Vector b
    for (int i = 0; i < totalSize; i++) {
    	cout << unif(rand_engine) << endl;
    }

	return 0;
}