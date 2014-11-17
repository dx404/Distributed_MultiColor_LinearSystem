#include <iostream>
#include <sstream>
#include <vector>

using namespace std;

const int blockBitWidth = 4;
inline int zmix (int i, int j, int k) {
    int zVal = 0;
    for (int mask = (1 << (blockBitWidth-1)); mask; mask >>= 1) {
        zVal = (i & mask) ? (zVal<<1) + 1 : zVal<<1;
        zVal = (j & mask) ? (zVal<<1) + 1 : zVal<<1;
        zVal = (k & mask) ? (zVal<<1) + 1 : zVal<<1;
    }
	return zVal;
}




inline int zmix2 (int i, int j, int k) {
	// return seqMix(i, j, k);
    int zVal = 0;

    stringstream ss_i;
	ss_i << (bitset<32>) i;
    int i_comb;
    ss_i >> oct >> i_comb;

    stringstream ss_j;
	ss_j << (bitset<32>) j;
    int j_comb;
    ss_j >> oct >> j_comb;

    stringstream ss_k;
	ss_k << (bitset<32>) k;
    int k_comb;
    ss_k >> oct >> k_comb;


    
	return (i_comb << 2) + (j_comb << 1) + k_comb;
}

inline vector<int> zDeMix(int zVal) {

	int i = 0, j = 0, k = 0;
	for (int mask = (1 << (3*blockBitWidth-1)); mask; mask >>= 3) {
        i = (zVal &  mask)     ? (i<<1) + 1 : i<<1;
        j = (zVal & (mask>>1)) ? (j<<1) + 1 : j<<1;
        k = (zVal & (mask>>2)) ? (k<<1) + 1 : k<<1;
    }

    
	return {i, j, k};
}

inline vector<int> zDeMix2 (int zVal) {
	// return seqMix(i, j, k);
    unsigned long int mask = 011111111111;
    bitset<32> bs;

    stringstream ss_i;
    ss_i << oct << ((zVal>>2) & mask);
    ss_i >> oct >> bs;
    int i = bs.to_ulong();

    stringstream ss_j;
    ss_j << oct << ((zVal>>1) & mask);
    ss_j >> oct >> bs;
    int j = bs.to_ulong();

    stringstream ss_k;
    ss_k << oct << (zVal & mask);
    ss_k >> oct >> bs;
    int k = bs.to_ulong();

	return {i, j, k};
}




int main() {
	cout << zmix (4, 5, 6) << endl;
	cout << zmix2 (4, 5, 6) << endl;
	zDeMix2 (458);
	vector<int> dm = zDeMix(458);
	vector<int> dm2 = zDeMix2(458);

	cout << dm[0]  << dm[1]  << dm[2]  << endl;
	cout << dm2[0] << dm2[1] << dm2[2] << endl;
 	return 0;
}