// normal 

#include <iostream>
#include <vector>
#include <cmath>

using namespace std;

const vector<double> coeff = {
	 1, 1, 1, 
	 1, 1, 1, 
	 1, 1, 1, 

	 1, 1, 1, 
	 1, -26, 1, 
	 1, 1, 1,

	 1, 1, 1, 
	 1, 1, 1, 
	 1, 1, 1,
};

double GS_one(vector<double> &x, const vector<double> &b, int i_bound, int j_bound, int k_bound) {
	double avgNormChange = 0;
	for (int i = 0; i < i_bound; i++) 
	for (int j = 0; j < j_bound; j++) 
	for (int k = 0; k < k_bound; k++) {
		int linearID = i*j_bound*k_bound + j*k_bound + k;
		double val = 0;
		for (int di = -1; di <= 1; di++) 
		for (int dj = -1; dj <= 1; dj++)
		for (int dk = -1; dk <= 1; dk++) 
			if ((di | dj | dk) && 
				i+di >= 0 && i+di < i_bound &&  
				j+dj >= 0 && j+dj < j_bound &&  
				k+dk >= 0 && k+dk < k_bound ){
				int linearAdjID = (i+di)*j_bound*k_bound + (j+dj)*k_bound + (k+dk);
				val += x[linearAdjID]; // add coefficent later
			}
		val = (b[linearID] - val) / (-26.0);
		avgNormChange += (val - x[linearID]) * (val - x[linearID]);
		x[linearID] = val;
	}
	avgNormChange = sqrt(avgNormChange) / (i_bound * j_bound * k_bound); 
	return avgNormChange;
}

int main(int argc, char *argv[]) {
	int i_bound, j_bound, k_bound;
	cin >> i_bound >> j_bound >> k_bound;
	int totalSize = i_bound * j_bound * k_bound;
	cout << endl;
	vector<double> x = {
2.04207	,
2.04307	,
2.04607	,
2.04707	,
2.05807	,
2.05907	,
2.06207	,
2.06307	};
	vector<double> b = {
1.4208	,
1.4308	,
1.4608	,
1.4708	,
1.5808	,
1.5908	,
1.6208	,
1.6308	};


	int n = 0;
	double err = 1;	
	err = GS_one (x, b, i_bound, j_bound, k_bound);

	cout << endl;
	for (auto val : x) {
		cout << val << endl;
	}
	cout << "err : " << err << endl;
	cout << "n : " << n << endl;
	return 0;
}
