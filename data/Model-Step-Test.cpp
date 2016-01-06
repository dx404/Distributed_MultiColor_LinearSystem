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

//
int main(int argc, char *argv[]) {
	int d =  atoi(argv[1]); // 2 << N;
	int ds = atoi(argv[2]); // 2 << ds;
	int T =  1 << atoi(argv[3]);
	int N = 1 << d;
	int sub_N = 1 << ds;
	int core_N = 1 << (d-ds);
	vector<double> data (N*N*N, 0); //--> 1
//	vector<double> NextData (N*N*N, 0); //--> 1
	for (int t = 0; t < T; t++){
		for (int i = 0; i < sub_N; i++)
		for (int j = 0; j < sub_N; j++)
		for (int k = 0; k < sub_N; k++) {
			for (int core_i = 0; core_i < core_N; core_i++)
			for (int core_j = 0; core_j < core_N; core_j++)
			for (int core_k = 0; core_k < core_N; core_k++) {
				int i_global = core_i*sub_N+i, j_global = core_j*sub_N+j, k_global = core_k*sub_N+k;
				int mix_ijk = getStoreIndex(i_global, j_global, k_global, N);
				double val = 0;
				for (int di=-1; di <=1; di++)
				for (int dj=-1; dj <=1; dj++)
				for (int dk=-1; dk <=1; dk++)
					if ( (di | dj | dk) &&
							i_global+di>=0 && i_global+di<N &&
							j_global+dj>=0 && j_global+dj<N &&
							k_global+dk>=0 && k_global+dk<N ){
						int mixIndex = getStoreIndex(i_global+di, j_global+dj, k_global+dk, N);
						val += data[mixIndex];
				}
				val = (cal_b(i_global, j_global, k_global, N) - val) / (-26.0);
				data[mix_ijk] = val;
			}
		}
//		copy(NextData.begin(), NextData.end(), data.begin());
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
	cout << N << ", " << core_N << ", " << sub_N << ", " << T << ", "<< error << endl;
	return 0;
}

