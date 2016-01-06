// read regular
// distribute
// Z-order load
// Z-order send
// tags
/*
mpicxx-mp -O3 -std=c++1y ./GS_OOP_v7.cpp
mpiexec-mp -np 10 ./a.out 2 2 2 1 dataSheet_128_128_128.txt 16
*/
#include <iostream>
#include <vector>
#include <string>
#include <cmath>
#include <sstream>
#include <fstream>
#include <mpi.h>
#include <map>
#include <unordered_map>

using namespace std;

#define send_data_tag 2001
#define return_data_tag 2002
#define exchange_data_tag 2003

#define max3(x, y, z) max(max((x), (y)), (z))
#define bitWidthDefault 32

int INIT_HALO_WIDTH = 1;

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

map<tuple<int, int, int>, int> memoZix;
map<int, tuple<int, int, int> > memoDeZix;

unordered_map<int, int> memoZixL;
unordered_map<int, int> memoDeZixL;

vector<int> memoInter (1024);
unordered_map<int, int> memoInterR;


void memoInterFill(int n) {
   for (int i = 0; i < n; i++) {
       int x = i, y = 0, k = 0;
       while (x) {
           y += (x & 1) << (3*k++);
           x >>= 1;
       }
       memoInter[i] = y;
       memoInterR[y] = i;
   }
}


class Capsule {
public:
	static int i_bound, j_bound, k_bound, totalSize;
	static int numCore_i, numCore_j, numCore_k, coresOnDuty, root_process;
	static int halo_width;

	int core_i, core_j, core_k, core_linear;
	int iStart, jStart, kStart; 
	int iBlockSize, jBlockSize, kBlockSize;
	int iExtBlockSize, jExtBlockSize, kExtBlockSize;
	int blockBitWidth, blockSize;

	vector<double> x_block, b_block; 
	vector<double> x_block2, b_block2; 

	int halo_noneUpdate_width;
	int halo_valid_width;

	double normChange;

private: 
	static MPI::Status status;

    int bitMask;
    int octMask_i, octMask_j, octMask_k;
    int i_bound_inter, j_bound_inter, k_bound_inter;

    int i_update_low_zmix, i_update_high_zmix;
    int j_update_low_zmix, j_update_high_zmix;
    int k_update_low_zmix, k_update_high_zmix;

	vector<vector<double> > sendHalo, recvHalo; 


public: 
	Capsule (int in_core_i, int in_core_j, int in_core_k) :
		core_i(in_core_i), core_j(in_core_j), core_k(in_core_k)
	{ 
		coresOnDuty  = numCore_i * numCore_j * numCore_k;
		root_process = numCore_i * numCore_j * numCore_k;
		core_linear = core_i * numCore_j * numCore_k + core_j * numCore_k + core_k;

		iBlockSize = (i_bound % numCore_i) ? i_bound / numCore_i + 1: i_bound / numCore_i;
		jBlockSize = (j_bound % numCore_j) ? j_bound / numCore_j + 1: j_bound / numCore_j;
		kBlockSize = (k_bound % numCore_k) ? k_bound / numCore_k + 1: k_bound / numCore_k;

		iStart = core_i * iBlockSize - halo_width;
		jStart = core_j * jBlockSize - halo_width;
		kStart = core_k * kBlockSize - halo_width;

		iExtBlockSize = iBlockSize + 2*halo_width;
		jExtBlockSize = jBlockSize + 2*halo_width; 
		kExtBlockSize = kBlockSize + 2*halo_width;
		blockBitWidth = ceil(max3(log2(iExtBlockSize), log2(jExtBlockSize), log2(kExtBlockSize)));
		blockSize = 1 << (3 * blockBitWidth); 

		bitMask = (1 << (3*blockBitWidth)) - 1;
	    octMask_i = 04444444444 & bitMask;
	    octMask_j = 02222222222 & bitMask;
	    octMask_k = 01111111111 & bitMask;

		i_bound_inter = memoInter[iExtBlockSize-1] << 2;
	    j_bound_inter = memoInter[jExtBlockSize-1] << 1;
	    k_bound_inter = memoInter[kExtBlockSize-1];
	    
		x_block = vector<double> (blockSize, 0);
		b_block = vector<double> (blockSize, 0);

		sendHalo = haloInit ();
		recvHalo = haloInit ();

		halo_noneUpdate_width = halo_width ? 1 : 0;
		halo_valid_width = halo_width;
	}


	vector<double> extractFrom (vector<double> &global_block) { // also applies to b
		vector<double> local_block (blockSize);
		for (int i_sub = 0; i_sub < iExtBlockSize; i_sub++) 
		for (int j_sub = 0; j_sub < jExtBlockSize; j_sub++) 
		for (int k_sub = 0; k_sub < kExtBlockSize; k_sub++) {
			if (iStart + i_sub < 0 || iStart + i_sub >= i_bound || 
				jStart + j_sub < 0 || jStart + j_sub >= j_bound || 
				kStart + k_sub < 0 || kStart + k_sub >= k_bound ) 
				continue;
			int srcPos = 
				  (iStart + i_sub) * j_bound * k_bound 
				+ (jStart + j_sub) * k_bound 
				+ (kStart + k_sub) ;
			int targetPos = zmix(i_sub, j_sub, k_sub);
			local_block[targetPos] = global_block[srcPos];
		}
		return local_block;
	}

	void prepareSendHalo () {
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
					int zID = zmix(i_sub, j_sub, k_sub);
					int subHaloID = (i_sub-iHaloStart) * jHaloWidth * kHaloWidth + (j_sub-jHaloStart) * kHaloWidth + (k_sub-kHaloStart);
					sendHalo[haloLinearID][subHaloID] = x_block[zID];
				}
		}
	}

	
	
	double GS_Z_block_update007_exterior () {
		double normChange = 0;

		int zAdjIDArray[27];
		int i_comp[3], j_comp[3], k_comp[3];

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
					int zID = zmix(i_sub, j_sub, k_sub);
					int iStep = 0, jStep = 0, kStep = 0;

					if ((zID & octMask_i) < i_update_low_zmix || 
						(zID & octMask_j) < j_update_low_zmix || 
						(zID & octMask_k) < k_update_low_zmix || 
						(zID & octMask_i) >= i_update_high_zmix || 
						(zID & octMask_j) >= j_update_high_zmix || 
						(zID & octMask_k) >= k_update_high_zmix ) {
						continue;
					}
					
				    if (zID & octMask_i){
				        i_comp[iStep++] = ((zID & octMask_i)  - 1)  & octMask_i;
				    }
				    i_comp[iStep++] = zID & octMask_i;
				    if ((zID & octMask_i) != i_bound_inter){
				         i_comp[iStep++] = ((zID | (~octMask_i)) + 1) & octMask_i;
				    }
				    
				    if (zID & octMask_j){
				        j_comp[jStep++] = ((zID & octMask_j)  - 1) & octMask_j;
				    }
				    j_comp[jStep++] = zID & octMask_j;
				    if ((zID & octMask_j) != j_bound_inter){
				        j_comp[jStep++] = ((zID | (~octMask_j)) + 1) & octMask_j;
				    }
				    
				    if (zID  & octMask_k){
				        k_comp[kStep++] = ((zID & octMask_k)  - 1) & octMask_k;
				    }
				    k_comp[kStep++] = zID & octMask_k;
				    if ((zID & octMask_k) != k_bound_inter){
				       k_comp[kStep++] = ((zID | (~octMask_k)) + 1) & octMask_k;
				    }

				   int t = 0;
				   for (int i = 0; i < iStep; i++)
				   for (int j = 0; j < jStep; j++)
				   for (int k = 0; k < kStep; k++) {
				   		int nzID = i_comp[i] | j_comp[j] | k_comp[k];
				    	if (nzID != zID)
				        	zAdjIDArray[t++] = nzID;
				    }

					double val = 0;
					for (int s = 0; s < t; s++) {
							val += x_block[zAdjIDArray[s]];
					}

					val = (b_block[zID] - val) / (-26.0);
					normChange += (val - x_block[zID]) * (val - x_block[zID]);
					x_block[zID] = val;
				}
		}
		normChange = sqrt(normChange) / (iExtBlockSize * jExtBlockSize * kExtBlockSize);
		return normChange;
	}

	double GS_Z_block_update007_interior (int di) { // di = 1 or - 1 end is exclusive
		double normChange = 0;
		int start = (di == 1) ?   0 : x_block.size();
		int end   = (di == -1)?  -1 : x_block.size() - 1;

		int zAdjIDArray[27];
		int i_comp[3], j_comp[3], k_comp[3];
		// update interior
		int boundary_noneUpdate_width = 2 * halo_width;
	    i_update_low_zmix  = memoInter[boundary_noneUpdate_width] << 2;
	   	j_update_low_zmix  = memoInter[boundary_noneUpdate_width] << 1;
	   	k_update_low_zmix  = memoInter[boundary_noneUpdate_width];

	    i_update_high_zmix = memoInter[iExtBlockSize - boundary_noneUpdate_width] << 2;
	    j_update_high_zmix = memoInter[jExtBlockSize - boundary_noneUpdate_width] << 1;
	    k_update_high_zmix = memoInter[kExtBlockSize - boundary_noneUpdate_width];


		for (int zID = start; zID != end; zID += di) { // reverse iterator
			int iStep = 0, jStep = 0, kStep = 0;

			if ((zID & octMask_i) < i_update_low_zmix || 
				(zID & octMask_j) < j_update_low_zmix || 
				(zID & octMask_k) < k_update_low_zmix || 
				(zID & octMask_i) >= i_update_high_zmix || 
				(zID & octMask_j) >= j_update_high_zmix || 
				(zID & octMask_k) >= k_update_high_zmix ) {
				continue;
			}
			
		    if (zID & octMask_i){
		        i_comp[iStep++] = ((zID & octMask_i)  - 1)  & octMask_i;
		    }
		    i_comp[iStep++] = zID & octMask_i;
		    if ((zID & octMask_i) != i_bound_inter){
		         i_comp[iStep++] = ((zID | (~octMask_i)) + 1) & octMask_i;
		    }
		    
		    if (zID & octMask_j){
		        j_comp[jStep++] = ((zID & octMask_j)  - 1) & octMask_j;
		    }
		    j_comp[jStep++] = zID & octMask_j;
		    if ((zID & octMask_j) != j_bound_inter){
		        j_comp[jStep++] = ((zID | (~octMask_j)) + 1) & octMask_j;
		    }
		    
		    if (zID  & octMask_k){
		        k_comp[kStep++] = ((zID & octMask_k)  - 1) & octMask_k;
		    }
		    k_comp[kStep++] = zID & octMask_k;
		    if ((zID & octMask_k) != k_bound_inter){
		       k_comp[kStep++] = ((zID | (~octMask_k)) + 1) & octMask_k;
		    }

		   int t = 0;
		   for (int i = 0; i < iStep; i++)
		   for (int j = 0; j < jStep; j++)
		   for (int k = 0; k < kStep; k++) {
		   		int nzID = i_comp[i] | j_comp[j] | k_comp[k];
		    	if (nzID != zID)
		        	zAdjIDArray[t++] = nzID;
		    }

			double val = 0;
			for (int s = 0; s < t; s++) {
					val += x_block[zAdjIDArray[s]];
			}

			val = (b_block[zID] - val) / (-26.0);
			normChange += (val - x_block[zID]) * (val - x_block[zID]);
			x_block[zID] = val;
			// cout << "::" << val << endl;
		}

		normChange = sqrt(normChange) / (iExtBlockSize * jExtBlockSize * kExtBlockSize);
		return normChange;

	}


	void haloSend () {
		if (halo_width == 0) 
			return;
		for (int di = -1; di <= 1; di++) 
		for (int dj = -1; dj <= 1; dj++)
		for (int dk = -1; dk <= 1; dk++) 
		if ((di | dj | dk) && 
			core_i + di >= 0 && core_i + di < numCore_i &&
			core_j + dj >= 0 && core_j + dj < numCore_j &&
			core_k + dk >= 0 && core_k + dk < numCore_k ) 
		{
			int destCoreLinear = (core_i + di) * numCore_j * numCore_k 
							   + (core_j + dj) * numCore_k 
							   + (core_k + dk);

			int srcHaloID  = (1+di)*9 + (1+dj)*3 + (1+dk);
			// int destHaloID = (1-di)*9 + (1-dj)*3 + (1-dk); // implicit
			int haloSize = sendHalo[srcHaloID].size();
			MPI::COMM_WORLD.Send (&sendHalo[srcHaloID][0], haloSize , MPI_DOUBLE, destCoreLinear, send_data_tag);
			
			
		}
	}

	void haloRecv() {
		if (halo_width == 0) 
			return;
		for (int di = -1; di <= 1; di++) 
		for (int dj = -1; dj <= 1; dj++)
		for (int dk = -1; dk <= 1; dk++) 
		if ((di | dj | dk) && 
			core_i + di >= 0 && core_i + di < numCore_i &&
			core_j + dj >= 0 && core_j + dj < numCore_j &&
			core_k + dk >= 0 && core_k + dk < numCore_k ) 
		{
			int srcCoreLinear = (core_i + di) * numCore_j * numCore_k 
								+ (core_j + dj) * numCore_k 
								+ (core_k + dk);
			//int srcHaloID  = (1-di)*9 + (1-dj)*3 + (1-dk); // implicit
			int destHaloID = (1+di)*9 + (1+dj)*3 + (1+dk); // Me
			
			int haloSize = recvHalo[destHaloID].size(); 
			if (haloSize) {
				MPI::COMM_WORLD.Recv(&recvHalo[destHaloID][0], haloSize, MPI_DOUBLE, srcCoreLinear, send_data_tag);
			}
		}
	}

	void haloExchange() {
		for (int di = -1; di <= 1; di++) 
		for (int dj = -1; dj <= 1; dj++)
		for (int dk = -1; dk <= 1; dk++) 
		if ((di | dj | dk) && 
			core_i + di >= 0 && core_i + di < numCore_i &&
			core_j + dj >= 0 && core_j + dj < numCore_j &&
			core_k + dk >= 0 && core_k + dk < numCore_k ) 
		{
			int adjCoreLinear = (core_i + di) * numCore_j * numCore_k 
							  + (core_j + dj) * numCore_k 
							  + (core_k + dk);

			int myHaloID  = (1+di)*9 + (1+dj)*3 + (1+dk);
			int haloSize = (int) sendHalo[myHaloID].size();
			if (haloSize) MPI::COMM_WORLD.Sendrecv (
				&sendHalo[myHaloID][0], haloSize, MPI_DOUBLE, adjCoreLinear, exchange_data_tag,
				&recvHalo[myHaloID][0], haloSize, MPI_DOUBLE, adjCoreLinear, exchange_data_tag
			);
			
		}
	}

	void updateFromRecvHalo () {
		for (int di = -1; di <= 1; di++) 
		for (int dj = -1; dj <= 1; dj++)
		for (int dk = -1; dk <= 1; dk++) 
			if (di | dj | dk) {
				int haloLinearID = (1+di) * 9 + (1+dj) * 3 + (1+dk);
				
				int iHaloStart = (di == -1) ? 0 : (di == 0) ? halo_width : iBlockSize + halo_width;
				int jHaloStart = (dj == -1) ? 0 : (dj == 0) ? halo_width : jBlockSize + halo_width;
				int kHaloStart = (dk == -1) ? 0 : (dk == 0) ? halo_width : kBlockSize + halo_width;

				int iHaloWidth = di ? halo_width : iBlockSize;
				int jHaloWidth = dj ? halo_width : jBlockSize;
				int kHaloWidth = dk ? halo_width : kBlockSize;

				for (int i_sub = iHaloStart; i_sub < iHaloStart + iHaloWidth; i_sub++) 
				for (int j_sub = jHaloStart; j_sub < jHaloStart + jHaloWidth; j_sub++) 
				for (int k_sub = kHaloStart; k_sub < kHaloStart + kHaloWidth; k_sub++) { 
					int zID = zmix(i_sub, j_sub, k_sub);
					int subHaloID = (i_sub - iHaloStart) * jHaloWidth * kHaloWidth + (j_sub - jHaloStart) * kHaloWidth + (k_sub - kHaloStart);
					x_block[zID] = recvHalo[haloLinearID][subHaloID];
				}
		}
	}

	void sendBackToRoot () {
		MPI::COMM_WORLD.Send (&x_block[0], blockSize , MPI_DOUBLE, root_process, send_data_tag);
		MPI::COMM_WORLD.Send (&this->normChange, 1 , MPI_DOUBLE, root_process, send_data_tag);
	}

	double syncBackAtRoot(vector<double> &x) {
		double normChange = 0.0;
		for (int i_sub = 0; i_sub < iExtBlockSize; i_sub++) 
		for (int j_sub = 0; j_sub < jExtBlockSize; j_sub++) 
		for (int k_sub = 0; k_sub < kExtBlockSize; k_sub++) {
			if (iStart + i_sub < 0 || iStart + i_sub >= i_bound || 
				jStart + j_sub < 0 || jStart + j_sub >= j_bound || 
				kStart + k_sub < 0 || kStart + k_sub >= k_bound ) 
			continue;
			int srcPos = 
				  (iStart + i_sub) * j_bound * k_bound 
				+ (jStart + j_sub) * k_bound 
				+ (kStart + k_sub) ;
			int targetPos = zmix(i_sub, j_sub, k_sub);
			normChange +=  (x_block[targetPos] - x[srcPos]) * (x_block[targetPos] - x[srcPos]);
			x[srcPos] = x_block[targetPos];
		}
		return normChange;
	}

	void betaSend() {
		int n = 1, my_id;
		int cores = numCore_i * numCore_j * numCore_k;
		// MPI_Send (&n, 1 , MPI_INT, (core_linear+1) % cores, send_data_tag, MPI_COMM_WORLD);
		MPI::COMM_WORLD.Send (&n, 1 , MPI_INT, (core_linear + 1) % cores, send_data_tag);

	}

private:
		// inner 
	vector<vector<double > > haloInit () {
		vector<vector<double > > vacantHalo(27);
		for (int di = -1; di <= 1; di++) 
		for (int dj = -1; dj <= 1; dj++)
		for (int dk = -1; dk <= 1; dk++) 
			if (di | dj | dk) { 
				int haloLinearID = (1+di)*9 + (1+dj)*3 + (1+dk);
				int haloWidth_i = di ? halo_width : iBlockSize;
				int haloWidth_j = dj ? halo_width : jBlockSize;
				int haloWidth_k = dk ? halo_width : kBlockSize;
				int haloSize = haloWidth_i * haloWidth_j * haloWidth_k;
				vacantHalo[haloLinearID] = vector<double> (haloSize, 0);
			}
		return vacantHalo;
	}

	inline void zNeighbors(vector<int> &zNeighborIndex, int centerIndex) {	    
	    vector<int> i_comp, j_comp, k_comp;
	    if (centerIndex & octMask_i)
	        i_comp.push_back(((centerIndex & octMask_i)  - 1) & octMask_i);
	    i_comp.push_back(centerIndex & octMask_i);
	    if ((centerIndex & octMask_i) != i_bound_inter)
	        i_comp.push_back(((centerIndex | (~octMask_i)) + 1) & octMask_i);
	        
	    if (centerIndex & octMask_j)
	        j_comp.push_back(((centerIndex & octMask_j)  - 1) & octMask_j);
	    j_comp.push_back(centerIndex & octMask_j);
	    if ((centerIndex & octMask_j) != j_bound_inter)
	        j_comp.push_back(((centerIndex | (~octMask_j)) + 1) & octMask_j);
	    
	    if (centerIndex  & octMask_k)
	        k_comp.push_back(((centerIndex & octMask_k)  - 1) & octMask_k);
	    k_comp.push_back(centerIndex & octMask_k);
	    if ((centerIndex & octMask_k) != k_bound_inter)
	        k_comp.push_back(((centerIndex | (~octMask_k)) + 1) & octMask_k);
	    
	    for (int i : i_comp)
	    for (int j : j_comp)
	    for (int k : k_comp) {
	    	if ((i | j | k) != centerIndex)
	        	zNeighborIndex.push_back(i | j | k);
	    }
	    	
	}


	inline int zmix (int i, int j, int k) {
		return seqMix(i, j, k);
		// int linearID = i * jExtBlockSize * kExtBlockSize + j * kExtBlockSize + k;
		// if (memoZixL.count(linearID)) return memoZixL[linearID];
	    
		return (memoInter[i]<<2) + (memoInter[j]<<1) + memoInter[k];


	    int zVal = 0;
	    for (unsigned long int mask = (1 << (blockBitWidth-1)); mask; mask >>= 1) {
	        zVal = (i & mask) ? (zVal<<1) + 1 : zVal<<1;
	        zVal = (j & mask) ? (zVal<<1) + 1 : zVal<<1;
	        zVal = (k & mask) ? (zVal<<1) + 1 : zVal<<1;
	    }
	    // int zVal2 = zmix2(i, j, k);
	    // memoDeZix[zVal] = {i, j, k};
	    // memoZix[{i, j, k}] = zVal;
	    // return seqMix(i, j, k);
		return zVal;
	}

	inline int zDeMix2 (int zVal) {
		return  memoInterR[(zVal>>2) & octMask_k] * jExtBlockSize * kExtBlockSize
				+ memoInterR[(zVal>>1) & octMask_k] * kExtBlockSize
				+ memoInterR[zVal & octMask_k]; 
	}

	inline vector<int> zDeMix(int zVal) {
		return seqDeMix(zVal); 
		// if (memoDeZixL.count(zVal)) { 
		// 	int seqVal = memoDeZixL[zVal];
		// 	int k =  seqVal % kExtBlockSize;
		// 	int j = (seqVal /= kExtBlockSize) % jExtBlockSize;
		// 	int i =  seqVal  / jExtBlockSize;
		// 	return {i, j, k};
		// }

		// return  {memoInterR[(zVal>>2) & octMask], memoInterR[(zVal>>1) & octMask], memoInterR[zVal & octMask]}; 



		int i = 0, j = 0, k = 0;
		for (unsigned long int mask = (1 << (3*blockBitWidth-1)); mask; mask >>= 3) {
	        i = (zVal &  mask)     ? (i<<1) + 1 : i<<1;
	        j = (zVal & (mask>>1)) ? (j<<1) + 1 : j<<1;
	        k = (zVal & (mask>>2)) ? (k<<1) + 1 : k<<1;
	    }
	    /*
	    vector<int> ijk2 = zDeMix2(zVal);
	    if (ijk2[0] != i || ijk2[1] != j || ijk2[2] != k) {
	    	printf("---> %i %i %i <--- \n", i, j, k);
	    }*/

	    // memoZix[{i, j, k}] = zVal; 
	    // memoDeZix[zVal] = {i, j, k};
	    // return seqDeMix(zVal); 
		return {i, j, k};
	}

	inline 
	int seqMix (int i, int j, int k) {
		return i * jExtBlockSize * kExtBlockSize + j * kExtBlockSize + k;
	}

	inline 
	vector<int> seqDeMix (int seqVal) {
		int k =  seqVal % kExtBlockSize;
		int j = (seqVal /= kExtBlockSize) % jExtBlockSize;
		int i =  seqVal  / jExtBlockSize;
		return {i, j, k};
	}


} ;

MPI::Status Capsule::status;
int Capsule::i_bound = 1;
int Capsule::j_bound = 1; 
int Capsule::k_bound = 1; 
int Capsule::totalSize = 1;
int Capsule::numCore_i = 1;
int Capsule::numCore_j = 1;
int Capsule::numCore_k = 1;
int Capsule::coresOnDuty = 1;
int Capsule::root_process = 1;
int Capsule::halo_width = 1;

int main (int argc, char *argv[]) {
	if (argc < 7) {
		cerr << "Please input : numCore_i, numCore_j, numCore_k, halo_width, inputFileName, iterationSteps" << endl;
		return 107;
	}

	MPI::Status status;
	int an_id, sender, ierr;

	MPI::Init(argc, argv);
	int my_id 	  = MPI::COMM_WORLD.Get_rank();
	int num_procs = MPI::COMM_WORLD.Get_size();

	int numCore_i  = Capsule::numCore_i  = atoi(argv[1]);
	int numCore_j  = Capsule::numCore_j  = atoi(argv[2]);
	int numCore_k  = Capsule::numCore_k  = atoi(argv[3]);
	int halo_width = Capsule::halo_width = atoi(argv[4]);
	INIT_HALO_WIDTH = halo_width;
	string dataSheetFileName = string(argv[5]);
	int iterationSteps = atoi(argv[6]);

	int num_procs_used = 1 + numCore_i * numCore_j * numCore_k;  // one for monitor and root process
	if (num_procs < num_procs_used) {
		cerr << "Not Enough Proccesses. " << endl;
		return 108;
	}
	int root_process = num_procs_used - 1;

	fstream dataSheet (dataSheetFileName);
	dataSheet >> Capsule::i_bound >> Capsule::j_bound >> Capsule::k_bound;
	int totalSize = Capsule::totalSize = Capsule::i_bound * Capsule::j_bound * Capsule::k_bound;

	memoInterFill(1024);

	if (my_id == root_process) {

		vector<double> x (totalSize, 0.0);
		vector<double> b (totalSize);
		
		for (int i = 0; i < x.size(); i++) {
			x[i] = 4 + 0.001 * i + 0.00007;
		}

		for (double &entry : b) {
			dataSheet >> entry;
		}
		double t0 = MPI::Wtime();
		
		for (int core_i = 0; core_i < numCore_i; core_i++)
		for (int core_j = 0; core_j < numCore_j; core_j++)
		for (int core_k = 0; core_k < numCore_k; core_k++) {
			vector<int> core_pack = {core_i, core_j, core_k};
			Capsule mainCull (core_i, core_j, core_k);
			vector<double> x_block = mainCull.extractFrom(x);
			vector<double> b_block = mainCull.extractFrom(b);

			MPI::COMM_WORLD.Send (&core_pack[0], 3 , MPI_INT, mainCull.core_linear, send_data_tag);
			MPI::COMM_WORLD.Send (&x_block[0], x_block.size(), MPI_DOUBLE, mainCull.core_linear, send_data_tag);
			MPI::COMM_WORLD.Send (&b_block[0], b_block.size(), MPI_DOUBLE, mainCull.core_linear, send_data_tag);

		}

		double normChange = 0.0, max_elapsed_secs = 0.0;
		for (int core_i = 0; core_i < numCore_i; core_i++)
		for (int core_j = 0; core_j < numCore_j; core_j++)
		for (int core_k = 0; core_k < numCore_k; core_k++) {
			double err = 0.0, elapsed_secs = 0.0;
			vector<int> core_pack = {core_i, core_j, core_k};
			Capsule mainCull (core_i, core_j, core_k);
			MPI::COMM_WORLD.Recv(&mainCull.x_block[0], mainCull.blockSize, MPI_DOUBLE, mainCull.core_linear, send_data_tag);
			MPI::COMM_WORLD.Recv(&err, 1, MPI_DOUBLE, mainCull.core_linear, send_data_tag);
			MPI::COMM_WORLD.Recv(&elapsed_secs, 1, MPI_DOUBLE, mainCull.core_linear, send_data_tag);
			normChange += err; 
			printf("core_linear: %d, err : %f\n", mainCull.core_linear, err);
			max_elapsed_secs = max (max_elapsed_secs, elapsed_secs);
			mainCull.syncBackAtRoot(x);
		}
		normChange /= (num_procs_used - 1.0);

		double t1 = MPI::Wtime();
		double elapsed_secs = t1 - t0;
		
		// printf("Core::%02i, seconds_elapsed : %f:: Error: %f, Steps : %i (%f, %f) \n", my_id, elapsed_secs, normChange, iterationSteps, t0, t1);
		printf("*#06#: %d, %d, %d: %d, %d, %d: %f,%f,%d\n",
				Capsule::i_bound, Capsule::j_bound, Capsule::k_bound,
				numCore_i, numCore_j, numCore_k,
				max_elapsed_secs, normChange, iterationSteps);
//		double err = 0;
//		for (double val : x) {
//			err += (val-1.0)*(val-1.0);
//		}
//		err = sqrt(err);
//		printf("Total Error: %f\n", err);
	}
	else if (my_id < root_process) { 
		vector<int> core_pack (3);
		MPI::COMM_WORLD.Recv(&core_pack[0], 3, MPI_INT, root_process, send_data_tag);
		Capsule capsule (core_pack[0], core_pack[1], core_pack[2]);
		
		MPI::COMM_WORLD.Recv (&capsule.x_block[0], capsule.blockSize, MPI_DOUBLE, root_process, send_data_tag);
		MPI::COMM_WORLD.Recv (&capsule.b_block[0], capsule.blockSize, MPI_DOUBLE, root_process, send_data_tag);
		
		double err = 1, max_err = 1;
		double t0 = MPI::Wtime();

		for (int i = 0; i < iterationSteps; i++) {
			err =  capsule.GS_Z_block_update007_exterior();
			if (capsule.halo_width == 1) {
				capsule.prepareSendHalo();
				capsule.haloExchange();
				capsule.halo_width = INIT_HALO_WIDTH;
			}
			else {
				capsule.halo_width--;
			}
			err += capsule.GS_Z_block_update007_interior(i & 1 ?  - 1 : 1); // fwd/bwd
//			err += capsule.GS_Z_block_update007_interior(1); // fwd
			
		}

		double t1 = MPI::Wtime();
    	double elapsed_secs = t1 - t0;
    	capsule.normChange = err;
		printf("Core: %02i, seconds_elapsed : %f:: Error: %f, Steps : %i (%f, %f) \n", my_id, elapsed_secs, err, iterationSteps, t0, t1);
		capsule.sendBackToRoot();
		MPI::COMM_WORLD.Send (&elapsed_secs, 1 , MPI_DOUBLE, root_process, send_data_tag);
		
	}
	MPI::Finalize();
	return 0;
}
