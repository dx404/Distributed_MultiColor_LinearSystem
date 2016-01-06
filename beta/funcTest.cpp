#include <iostream>
#include <vector>
#include <string>
#include <cmath>
#include <fstream>
using namespace std;

int zmix(int i, int j, int k, int bitWidth) {
    int zVal = 0;
    for (int mask = (1 << (bitWidth-1)); mask; mask >>= 1) {
        zVal = (i & mask) ? (zVal<<1) + 1 : zVal<<1;
        zVal = (j & mask) ? (zVal<<1) + 1 : zVal<<1;
        zVal = (k & mask) ? (zVal<<1) + 1 : zVal<<1;
    }
	return zVal;
}

vector<int> zDeMix(int zVal, int bitWidth) {
	int i = 0, j = 0, k = 0;
	for (int mask = (1 << (3*bitWidth-1)); mask; mask >>= 3) {
        i = (zVal &  mask)     ? (i<<1) + 1 : i<<1;
        j = (zVal & (mask>>1)) ? (j<<1) + 1 : j<<1;
        k = (zVal & (mask>>2)) ? (k<<1) + 1 : k<<1;
    }
	return {i, j, k};
}

vector<vector<double > > haloInit (int iBlockSize, int jBlockSize, int kBlockSize, int halo_width) {
	vector<vector<double > > halo(27);
	for (int di = -1; di <= 1; di++) 
	for (int dj = -1; dj <= 1; dj++)
	for (int dk = -1; dk <= 1; dk++) 
		if (di | dj | dk) { 
			int haloLinearID = (1+di)*9 + (1+dj)*3 + (1+dk);
			int haloWidth_i = di ? halo_width : iBlockSize;
			int haloWidth_j = dj ? halo_width : jBlockSize;
			int haloWidth_k = dk ? halo_width : kBlockSize;
			int haloSize = haloWidth_i * haloWidth_j * haloWidth_k;
			halo[haloLinearID] = vector<double> (haloSize);
		}
	return halo;
}

void prepareSendHalo (
	vector<vector<double> > &sendHalo, const vector<double> &x, 
	int iBlockSize, int jBlockSize, int kBlockSize, int halo_width, int blockBitWidth) {
	for (int di = -1; di <= 1; di++) 
	for (int dj = -1; dj <= 1; dj++)
	for (int dk = -1; dk <= 1; dk++) 
		if (di | dj | dk) {
			int haloLinearID = (1+di) * 9 + (1+dj) * 3 + (1+dk);
			int iHaloStart = (di != 1) ? halo_width : iBlockSize;
			int jHaloStart = (dj != 1) ? halo_width : jBlockSize;
			int kHaloStart = (dk != 1) ? halo_width : kBlockSize;

			int iHaloWidth = di ? halo_width : iBlockSize;
			int jHaloWidth = dj ? halo_width : jBlockSize;
			int kHaloWidth = dk ? halo_width : kBlockSize;

			for (int i_sub = iHaloStart; i_sub < iHaloStart + iHaloWidth; i_sub++) 
			for (int j_sub = jHaloStart; j_sub < jHaloStart + jHaloWidth; j_sub++) 
			for (int k_sub = kHaloStart; k_sub < kHaloStart + kHaloWidth; k_sub++) { 
				int zID = zmix(i_sub, j_sub, k_sub, blockBitWidth);
				int subHaloID = (i_sub-iHaloStart) * jHaloWidth * kHaloWidth + (j_sub-jHaloStart) * kHaloWidth + (k_sub-kHaloStart);
				sendHalo[haloLinearID][subHaloID] = x[zID];
			}
	}
}

int main() {
	vector<double> data = 
{1.0008,
1.0108	,
1.0208	,
1.0308	,
1.0408	,
1.0508	,
1.0608	,
1.0708	,
1.0808	,
1.0908	,
1.1008	,
1.1108	,
1.1208	,
1.1308	,
1.1408	,
1.1508	,
1.1608	,
1.1708	,
1.1808	,
1.1908	,
1.2008	,
1.2108	,
1.2208	,
1.2308	,
1.2408	,
1.2508	,
1.2608};

	vector<vector<double>> halo = haloInit(1, 1, 1, 1);
	prepareSendHalo (halo, data, 1, 1, 1, 1, 2); 
	




















	return 0;
}