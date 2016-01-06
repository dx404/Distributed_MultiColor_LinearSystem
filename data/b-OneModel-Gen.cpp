#include <iostream>
#include <cstdio>
#include <vector>
#include <cmath>
using namespace std;

inline double cal_b (int i, int j, int k, int N) {
	int countOnSurface =
			(i == 0 || i == N - 1)
		  + (j == 0 || j == N - 1)
		  + (k == 0 || k == N - 1);
	vector<int> coeff = {0.0, -9.0, -15.0, -19.0}; // inner, surface, edge, corner
	return coeff[countOnSurface];
}

inline int getStoreIndex (int i, int j, int k, int N) {
	return i*N*N + j*N + k;
}

//
int main(int argc, char *argv[]) {
	int N =  atoi(argv[1]);
	cout << N << " " << N << " " << N << " " << endl;
	cout << endl;
	for (int i = 0; i < N; i++)
	for (int j = 0; j < N; j++)
	for (int k = 0; k < N; k++) {
		double b_ijk = cal_b(i, j, k, N);
		cout << b_ijk << endl;
	}
	return 0;
}

