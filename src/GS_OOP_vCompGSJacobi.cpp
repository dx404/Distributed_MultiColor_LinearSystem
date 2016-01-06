/*
 * hello.cpp
 *
 *  Created on: Feb 10, 2015
 *      Author: duozhao
 */

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
	vector<int> coeff = {0, -9, -15, -19}; // inner, surface, edge, corner
	return coeff[countOnSurface];
}

inline int getStoreIndex (int i, int j, int k, int N) {
	return i*N*N + j*N + k;
}

int main(int argc, char *argv[]) {
	int d =  4; // 2 << N;

	int N = 1 << d;
	
	int T = 10;
	vector<double> data (N*N*N, 0); //--> 1
	vector<double> NextData (N*N*N, 0); //--> 1
	for (int t = 0; t < T; t++){
		for (int i = 0; i < N; i++)
		for (int j = 0; j < N; j++)
		for (int k = 0; k < N; k++) {
			int mix_ijk = getStoreIndex(i,j,k,N);
			double val = 0;
			for (int di=-1; di <=1; di++)
			for (int dj=-1; dj <=1; dj++)
			for (int dk=-1; dk <=1; dk++)
				if ( (di | dj | dk) &&
				  i+di>=0 && i+di<N &&
				  j+dj>=0 && j+dj<N &&
				  k+dk>=0 && k+dk<N ){
					int mixIndex = getStoreIndex(i+di, j+dj, k+dk, N);
					val += data[mixIndex];
			}
			val = (cal_b(i, j, k, N) - val) / (-26.0);
			data[mix_ijk] = val;
		}
	}
	double error = 0;
	for (int i = 0; i < N; i++)
	for (int j = 0; j < N; j++)
	for (int k = 0; k < N; k++) {
		int mix_ijk = getStoreIndex(i,j,k,N);
		double val = data[mix_ijk];
		error += (val - 1.0) * (val - 1.0);
//		printf("(%d, %d, %d) : %f \n", i, j, k, val);
	}
	error = sqrt(error/(N*N*N));
	cout << error << endl;
	return 0;
}



