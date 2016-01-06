#include <iostream>
#include <vector>
#include <string>
#include <cmath>
#include <sstream>
#include <fstream>
#include <unordered_map>
using namespace std;
#define max3(x, y, z) max(max((x), (y)), (z))
#define ceil_div(x, y) ((x) % (y)) ? (((x)/(y))+1) : ((x)/(y))

class Capsule {
public:
	int halo_width;
	int iBlockSize, jBlockSize, kBlockSize, blockSize; // physical
	int iExtBlockSize, jExtBlockSize, kExtBlockSize; //
	int iExtBitWidth, jExtBitWidth, kExtBitWidth;
	int mixExtBitWidth; // include halo Extension
	int iPadBlockSize, jPadBlockSize, kPadBlockSize, padBlockSize; // include halo Extension
	unsigned int iMask, jMask, kMask;
	vector<unsigned int> i_interleave, j_interleave, k_interleave;

	int iMax_interleave, jMax_interleave, kMax_interleave;

	vector<double> x_block, b_block;
	vector<vector<double> > A_block;
	vector<vector<int> > Acol_block;

	double normChange; // for error computation

	Capsule(int iW, int jW, int kW) {
		iBlockSize = iW;
		jBlockSize = jW;
		kBlockSize = kW;
		initCommonStruct();
		prepareIndexMapping();
		x_block = vector<double> (padBlockSize, 0);
		b_block = vector<double> (padBlockSize, 0);
		A_block = vector<vector<double> > (padBlockSize, vector<double> (27, 3.14));
		Acol_block = vector<vector<int> > (padBlockSize, vector<int> (27, -1));
		halo_width = 0;
		normChange = 0;

		for (int i = 0; i < padBlockSize; i++) {
			for (int j = 0; j < 27; j++) {
				A_block[i][j] += 0.1*(i+2)*(j-7);
			}
		}

	}

	void initCommonStruct(){
		blockSize = iBlockSize * jBlockSize * kBlockSize;

		iExtBlockSize = iBlockSize + 2*halo_width;
		jExtBlockSize = jBlockSize + 2*halo_width;
		kExtBlockSize = kBlockSize + 2*halo_width;

		iExtBitWidth = ceil(log2(iExtBlockSize));
		jExtBitWidth = ceil(log2(jExtBlockSize));
		kExtBitWidth = ceil(log2(kExtBlockSize));
		mixExtBitWidth = iExtBitWidth + jExtBitWidth + kExtBitWidth;

		iPadBlockSize = 1 << iExtBitWidth;
		jPadBlockSize = 1 << jExtBitWidth;
		kPadBlockSize = 1 << kExtBitWidth;
		padBlockSize = iPadBlockSize * jPadBlockSize * kPadBlockSize;

	}

	void prepareIndexMapping(){
		iMask = jMask = kMask = 0; // to mask the bits in each dimension
		i_interleave.resize(iPadBlockSize);
		j_interleave.resize(jPadBlockSize);
		k_interleave.resize(kPadBlockSize);
		i_interleave[0] = j_interleave[0] = k_interleave[0] = 0;

		int s = iExtBitWidth - 1;
	    for (; s >= jExtBitWidth; s--) {
	        iMask = (iMask << 1) + 01;
	    }
	    for (; s >= kExtBitWidth; s--) {
	        iMask = (iMask << 2) + 02;
	        jMask = (jMask << 2) + 01;
	    }
	    for (; s >= 0; s--) {
	        iMask = (iMask << 3) + 04;
	        jMask = (jMask << 3) + 02;
	        kMask = (kMask << 3) + 01;
	    }
	    unsigned int y = ~iMask;
	    for (int x = 1; x < iPadBlockSize; x++) {
	        y = (y+1) | ~iMask;
	        i_interleave[x] = y & iMask;
	    }
	    y = ~jMask;
	    for (int x = 1; x < jPadBlockSize; x++) {
	        y = (y+1) | ~jMask;
	        j_interleave[x] = y & jMask;
	    }
	    y = ~kMask;
	    for (int x = 1; x < kPadBlockSize; x++) {
	        y = (y+1) | ~kMask;
	        k_interleave[x] = y & kMask;
	    }

	    iMax_interleave = i_interleave[iExtBlockSize-1];
	    jMax_interleave = j_interleave[jExtBlockSize-1];
	    kMax_interleave = k_interleave[kExtBlockSize-1];

	}

	inline int mortonIndex (int i, int j, int k) {
		return i_interleave[i] | j_interleave [j] | k_interleave[k];
	}

	/*
	 * di is either 1 or -1, the direction
	 * w is the non-update width
	 */
	void addr_GS_block_interior_update_009 (int di, int w) {
			int start, end;
			if (di == 1) {
				start = 0;
				end   = x_block.size();
			}
			else {
				start = x_block.size() - 1;
				end   = -1;
			}

			int adjIDArray[27], iPart[3], jPart[3], kPart[3];
			// update interior

			int iLowMorton = i_interleave[w];
			int jLowMorton = j_interleave[w];
			int kLowMorton = k_interleave[w];

			int iHighMorton = i_interleave[iExtBlockSize - w];
			int jHighMorton = j_interleave[jExtBlockSize - w];
			int kHighMorton = k_interleave[kExtBlockSize - w];

			for (int zID = start; zID != end; zID += di) { // reverse iterator
				int iStep = 0, jStep = 0, kStep = 0;
				// out of user-specified boundary
				if ((zID & iMask) < iLowMorton ||
					(zID & jMask) < jLowMorton ||
					(zID & kMask) < kLowMorton ||
					(iHighMorton && (zID & iMask) >= iHighMorton) ||
					(jHighMorton && (zID & jMask) >= jHighMorton) ||
					(kHighMorton && (zID & kMask) >= kHighMorton) ) {
					continue;
				}

				if (zID & iMask){
					iPart[iStep++] = ((zID & iMask)  - 1)  & iMask;
				}
				iPart[iStep++] = zID & iMask;
				if ((zID & iMask) != iMax_interleave){
					 iPart[iStep++] = ((zID | (~iMask)) + 1) & iMask;
				}

				if (zID & jMask){
					jPart[jStep++] = ((zID & jMask)  - 1) & jMask;
				}
				jPart[jStep++] = zID & jMask;
				if ((zID & jMask) != jMax_interleave){
					jPart[jStep++] = ((zID | (~jMask)) + 1) & jMask;
				}

				if (zID  & kMask){
					kPart[kStep++] = ((zID & kMask)  - 1) & kMask;
				}
				kPart[kStep++] = zID & kMask;
				if ((zID & kMask) != kMax_interleave){
				   kPart[kStep++] = ((zID | (~kMask)) + 1) & kMask;
				}

			   int t = 0;
			   for (int i = 0; i < iStep; i++)
			   for (int j = 0; j < jStep; j++)
			   for (int k = 0; k < kStep; k++) {
					int adjID = iPart[i] | jPart[j] | kPart[k];
					if (adjID != zID)
						adjIDArray[t++] = adjID;
				}

				double val = 0;
				for (int s = 0; s < t; s++) {
						val += x_block[adjIDArray[s]] * A_block[adjIDArray[s]][s];
						Acol_block[zID][s] = adjIDArray[s];
				}

			}
		}

	void GS_indirect_update_010() {
		for (int i = 0; i < padBlockSize; i++) {
			double val = 0;
			for (int j = 0; j < 27; j++) {
				int nID = Acol_block[i][j];
				if (nID == -1) break;
				if (i != nID) val += A_block[i][nID] * x_block[nID];
			}
			val = (b_block[i] - val) / (-26.0 * A_block[i][13]);
			x_block[i] = val;
		}
	}

	void GS_block_interior_update_009 (int di, int w) {
		int start, end;
		if (di == 1) {
			start = 0;
			end   = x_block.size();
		}
		else {
			start = x_block.size() - 1;
			end   = -1;
		}

		int adjIDArray[27], iPart[3], jPart[3], kPart[3];
		// update interior

		int iLowMorton = i_interleave[w];
		int jLowMorton = j_interleave[w];
		int kLowMorton = k_interleave[w];

		int iHighMorton = i_interleave[iExtBlockSize - w];
		int jHighMorton = j_interleave[jExtBlockSize - w];
		int kHighMorton = k_interleave[kExtBlockSize - w];

		for (int zID = start; zID != end; zID += di) { // reverse iterator
			int iStep = 0, jStep = 0, kStep = 0;
			// out of user-specified boundary
			if ((zID & iMask) < iLowMorton ||
				(zID & jMask) < jLowMorton ||
				(zID & kMask) < kLowMorton ||
				(iHighMorton && (zID & iMask) >= iHighMorton) ||
				(jHighMorton && (zID & jMask) >= jHighMorton) ||
				(kHighMorton && (zID & kMask) >= kHighMorton) ) {
				continue;
			}

			if (zID & iMask){
				iPart[iStep++] = ((zID & iMask)  - 1)  & iMask;
			}
			iPart[iStep++] = zID & iMask;
			if ((zID & iMask) != iMax_interleave){
				 iPart[iStep++] = ((zID | (~iMask)) + 1) & iMask;
			}

			if (zID & jMask){
				jPart[jStep++] = ((zID & jMask)  - 1) & jMask;
			}
			jPart[jStep++] = zID & jMask;
			if ((zID & jMask) != jMax_interleave){
				jPart[jStep++] = ((zID | (~jMask)) + 1) & jMask;
			}

			if (zID  & kMask){
				kPart[kStep++] = ((zID & kMask)  - 1) & kMask;
			}
			kPart[kStep++] = zID & kMask;
			if ((zID & kMask) != kMax_interleave){
			   kPart[kStep++] = ((zID | (~kMask)) + 1) & kMask;
			}

		   int t = 0;
		   for (int i = 0; i < iStep; i++)
		   for (int j = 0; j < jStep; j++)
		   for (int k = 0; k < kStep; k++) {
				int adjID = iPart[i] | jPart[j] | kPart[k];
				if (adjID != zID)
					adjIDArray[t++] = adjID;
			}

			double val = 0;
			for (int s = 0; s < t; s++) {
					val += x_block[adjIDArray[s]] * A_block[adjIDArray[s]][s];
			}

			val = (b_block[zID] - val) / (-26.0 * A_block[zID][13]);
			normChange += (val - x_block[zID]) * (val - x_block[zID]);
			x_block[zID] = val;
		}
	}

	// other block_wise

	double GS_Regular_block_update_02 (int innerCubeSize) {
		double normChange = 0;

		for (int i_block = 0; i_block < jPadBlockSize/innerCubeSize; i_block++)
		for (int j_block = 0; j_block < kPadBlockSize/innerCubeSize; j_block++)
		for (int k_block = 0; k_block < kPadBlockSize/innerCubeSize; k_block++)


		for (int i_sub = i_block*innerCubeSize; i_sub < (i_block+1)*innerCubeSize; i_sub++)
		for (int j_sub = j_block*innerCubeSize; j_sub < (j_block+1)*innerCubeSize; j_sub++)
		for (int k_sub = k_block*innerCubeSize; k_sub < (k_block+1)*innerCubeSize; k_sub++) {
			int zID =  (i_sub)*jPadBlockSize*kPadBlockSize
					+ (j_sub)*kPadBlockSize
					+ (k_sub);
			if (i_sub >= iExtBlockSize     || j_sub >= jExtBlockSize     || k_sub >= kExtBlockSize )
				continue;

			double val = 0;
			for (int di = -1; di <= 1; di++)
			for (int dj = -1; dj <= 1; dj++)
			for (int dk = -1; dk <= 1; dk++) if ((di | dj | dk) &&
					i_sub+di >= 0 && i_sub+di < iExtBlockSize &&
					j_sub+dj >= 0 && j_sub+dj < jExtBlockSize &&
					k_sub+dk >= 0 && k_sub+dk < kExtBlockSize ){
					int zAdjID = (i_sub+di)*jPadBlockSize*kPadBlockSize
							+ (j_sub+dj)*kPadBlockSize
							+ (k_sub+dk);
					int colID = (1+di)*9 + (1+dj)*3 + (1+dk);
//					val += x_block[zAdjID] * A_block[zAdjID][colID]; // add coefficent later
					val += x_block[zAdjID]; // add coefficent later
				}

//			val = (b_block[zID] - val) / (-26.0 * A_block[zID][13]);
			val = (b_block[zID] - val) / (-26.0 );
			normChange += (val - x_block[zID]) * (val - x_block[zID]);
			x_block[zID] = val;
		}
		return this->normChange = normChange;
	}

};


int main(int argc, char *argv[]) {
	string dataSheetFileName(argv[1]);
	fstream dataSheet (dataSheetFileName);
	int iBlockSize, jBlockSize, kBlockSize, blockSize; // physical
	dataSheet >> iBlockSize >> jBlockSize >> kBlockSize;
	Capsule *cube = new Capsule(iBlockSize, jBlockSize, kBlockSize);


	for (int i = 0; i < iBlockSize; i++) {
		for (int j = 0; j < jBlockSize; j++) {
			for (int k = 0; k < kBlockSize; k++) {
				int zID = cube->mortonIndex(i, j, k);
				dataSheet >> cube->b_block[zID];
			}
		}
	}
	/*
	for (int i = 0; i < iBlockSize; i++) {
		for (int j = 0; j < jBlockSize; j++) {
			for (int k = 0; k < kBlockSize; k++) {
				int zID =
						i * cube->iPadBlockSize * cube->jPadBlockSize
						+ j * cube->kPadBlockSize + k;
//				cout << i << "\t" << j << "\t" << k << "\t" << zID << endl;
				dataSheet >> cube->b_block[zID];
			}
		}
	}
	 */
	int subCubeSize = atoi(argv[2]);
	int t0 = clock();
	cube->addr_GS_block_interior_update_009(1, subCubeSize);
	cube->GS_indirect_update_010();
//	cube->GS_block_interior_update_009(1, subCubeSize);
//	cube->GS_Regular_block_update_02(subCubeSize);
//	cube->GS_Z_block_update_02(subCubeSize);
	int t1 = clock();
//	cout << endl;
//	cout << cube->padBlockSize << endl;
//	cout << oct << cube->iMask<< endl;
//	cout << cube->jMask<< endl;
//	cout << cube->kMask<< endl;
//
//	cout << cube->i_interleave.size() << endl;
//	cout << cube->j_interleave.size() << endl;
//	cout << cube->k_interleave.size() << endl;

	double t = 1.0 * (t1-t0) / CLOCKS_PER_SEC;
	cout << "Time Elapsed, " <<subCubeSize << "," << t << endl;

	delete cube;
	dataSheet.close();
	return 0;
}
