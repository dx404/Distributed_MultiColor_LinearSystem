
/*
 * For a more generic dimensional in each  direction
 * mpicxx-mp -O3 -std=c++1y ./GS_OOP_v9.cpp
 * mpiexec-mp -np 10 ./a.out 2 2 2 1 dataSheet_128_128_128.txt 16
*/
#include <iostream>
#include <vector>
#include <string>
#include <cmath>
#include <sstream>
#include <fstream>
#include <mpi.h>
#include <unordered_map>

using namespace std;

#define send_data_tag 2001
#define return_data_tag 2002
#define exchange_data_tag 2003

#define max3(x, y, z) max(max((x), (y)), (z))
#define ceil_div(x, y) ((x) % (y)) ? (((x)/(y))+1) : ((x)/(y))

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


class Capsule {
public:
	//Geometry
	static int iGlobalBound, jGlobalBound, kGlobalBound, globalSize;
	static int numCore_i, numCore_j, numCore_k;
	static int iBlockSize, jBlockSize, kBlockSize, blockSize; // physical
	static int halo_width;
	static int iExtBlockSize, jExtBlockSize, kExtBlockSize; //
	static int iExtBitWidth, jExtBitWidth, kExtBitWidth;
	static int mixExtBitWidth; // include halo Extension
	static int iPadBlockSize, jPadBlockSize, kPadBlockSize, padBlockSize; // include halo Extension
	static unsigned int iMask, jMask, kMask;
	static vector<unsigned int> i_interleave, j_interleave, k_interleave;

	int core_i, core_j, core_k, core_linear;
	int iStart, jStart, kStart; // for a global start
	int iExtStart, jExtStart, kExtStart; // for a global start
	int iMax_interleave, jMax_interleave, kMax_interleave;

	vector<double> x_block, b_block; 
	vector<vector<double> > sendHalo, recvHalo;

	int halo_noneUpdate_width;
	int halo_valid_width;

	double normChange; // for error computation

	static void initCommonStruct(){
		// [0, i_bound-1] by [0, j_bound-1] by [0, k_bound-1]
		iBlockSize = ceil_div(i_bound, numCore_i);
		jBlockSize = ceil_div(j_bound, numCore_j);
		kBlockSize = ceil_div(k_bound, numCore_k);
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

	static void prepareIndexMapping(){
		iMask = jMask = kMask = 0; // to mask the bits in each dimension
		i_interleave.resize(iExtBitWidth);
		j_interleave.resize(jExtBitWidth);
		k_interleave.resize(kExtBitWidth);
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
	    for (int x = 1; x < iExtBitWidth; x++) {
	        y = (y+1) | ~iMask;
	        i_interleave[x] = y & iMask;
	    }
	    y = ~jMask;
	    for (int x = 1; x < jExtBitWidth; x++) {
	        y = (y+1) | ~jMask;
	        j_interleave[x] = y & jMask;
	    }
	    y = ~kMask;
	    for (int x = 1; x < kExtBitWidth; x++) {
	        y = (y+1) | ~kMask;
	        k_interleave[x] = y & kMask;
	    }

	    iMax_interleave = i_interleave[iExtBlockSize-1];
	    jMax_interleave = j_interleave[jExtBlockSize-1];
	    kMax_interleave = k_interleave[kExtBlockSize-1];

	}

	inline static int mortonIndex (int i, int j, int k) {
		return i_interleave[i] | j_interleave [j] | k_interleave[k];
	}

private: 
	static MPI::Status status;

    int i_bound_inter, j_bound_inter, k_bound_inter;

    int i_update_low_zmix, i_update_high_zmix;
    int j_update_low_zmix, j_update_high_zmix;
    int k_update_low_zmix, k_update_high_zmix;



public: 
	Capsule (int in_core_i, int in_core_j, int in_core_k) :
		core_i(in_core_i), core_j(in_core_j), core_k(in_core_k)
	{ 
		core_linear = core_i * numCore_j * numCore_k + core_j * numCore_k + core_k;

		iStart = core_i * iBlockSize;
		jStart = core_j * jBlockSize;
		kStart = core_k * kBlockSize;

		iExtStart = iStart - halo_width;
		jExtStart = jStart - halo_width;
		kExtStart = kStart - halo_width;
	    
		x_block = vector<double> (padBlockSize, 0);
		b_block = vector<double> (padBlockSize, 0);

		sendHalo = haloInit ();
		recvHalo = haloInit ();

		halo_noneUpdate_width = halo_width ? 1 : 0;
		halo_valid_width = halo_width;
	}


	void loadFromGlobal (vector<double> &local_block, vector<double> &global_block) {
		for (int iLocal = 0; iLocal < iExtBlockSize; iLocal++)
		for (int jLocal = 0; jLocal < jExtBlockSize; jLocal++)
		for (int kLocal = 0; kLocal < kExtBlockSize; kLocal++) {
			int iGlobal = iExtStart + iLocal;
			int jGlobal = jExtStart + jLocal;
			int kGlobal = kExtStart + kLocal;
			if (iGlobal < 0 || iGlobal >= iGlobalBound ||
				jGlobal < 0 || jGlobal >= jGlobalBound ||
				kGlobal < 0 || kGlobal >= kGlobalBound )
				continue;
			int globalPos =
				  iGlobal * jGlobalBound * kGlobalBound
				+ jGlobal * kGlobalBound
				+ kGlobal ;
			int localPos = mortonIndex(iLocal, jLocal, kLocal);
			local_block[localPos] = global_block[globalPos];
		}
	}

	void syncToGlobal(vector<double> &local_block, vector<double> &global_block) {
		for (int iLocal = 0; iLocal < iExtBlockSize; iLocal++)
		for (int jLocal = 0; jLocal < jExtBlockSize; jLocal++)
		for (int kLocal = 0; kLocal < kExtBlockSize; kLocal++) {
			int iGlobal = iExtStart + iLocal;
			int jGlobal = jExtStart + jLocal;
			int kGlobal = kExtStart + kLocal;
			if (iGlobal < 0 || iGlobal >= iGlobalBound ||
				jGlobal < 0 || jGlobal >= jGlobalBound ||
				kGlobal < 0 || kGlobal >= kGlobalBound )
				continue;
			int globalPos =
				  iGlobal * jGlobalBound * kGlobalBound
				+ jGlobal * kGlobalBound
				+ kGlobal ;
			int localPos = mortonIndex(iLocal, jLocal, kLocal);
			global_block[globalPos] = local_block[localPos];
		}
	}


	void prepareSendHalo () {
		for (int di = -1; di <= 1; di++) 
		for (int dj = -1; dj <= 1; dj++)
		for (int dk = -1; dk <= 1; dk++) if (di | dj | dk) {
			int haloLinearID = (1+di) * 9 + (1+dj) * 3 + (1+dk);
			int iHaloStart = (di != 1) ? halo_width : iBlockSize;
			int jHaloStart = (dj != 1) ? halo_width : jBlockSize;
			int kHaloStart = (dk != 1) ? halo_width : kBlockSize;

			int iHaloWidth = di ? halo_width : iBlockSize;
			int jHaloWidth = dj ? halo_width : jBlockSize;
			int kHaloWidth = dk ? halo_width : kBlockSize;

			for (int iLocal = iHaloStart; iLocal < iHaloStart + iHaloWidth; iLocal++)
			for (int jLocal = jHaloStart; jLocal < jHaloStart + jHaloWidth; jLocal++)
			for (int kLocal = kHaloStart; kLocal < kHaloStart + kHaloWidth; kLocal++) {
				int localMassPos = mortonIndex(iLocal, jLocal, kLocal);
				int innerHaloPos =
							(iLocal-iHaloStart) * jHaloWidth * kHaloWidth
							  + (jLocal-jHaloStart) * kHaloWidth
							  + (kLocal-kHaloStart);
				sendHalo[haloLinearID][innerHaloPos] = x_block[localMassPos];
			}
		}
	}
	
	void GS_block_exterior_update_009 (int di, int w) {
		int adjIDArray[27], iPart[3], jPart[3], kPart[3];

		int iLowMorton = i_interleave[halo_width];
		int jLowMorton = j_interleave[halo_width];
		int kLowMorton = k_interleave[halo_width];

		int iHighMorton = i_interleave[iExtBlockSize - w];
		int jHighMorton = j_interleave[jExtBlockSize - w];
		int kHighMorton = k_interleave[kExtBlockSize - w];

		for (int di = -1; di <= 1; di++) 
		for (int dj = -1; dj <= 1; dj++)
		for (int dk = -1; dk <= 1; dk++)  if (di | dj | dk) {
			int haloLinearID = (1+di) * 9 + (1+dj) * 3 + (1+dk);
			int iHaloStart = (di != 1) ? halo_width : iBlockSize;
			int jHaloStart = (dj != 1) ? halo_width : jBlockSize;
			int kHaloStart = (dk != 1) ? halo_width : kBlockSize;

			int iHaloWidth = di ? halo_width : iBlockSize;
			int jHaloWidth = dj ? halo_width : jBlockSize;
			int kHaloWidth = dk ? halo_width : kBlockSize;

			for (int iLocal = iHaloStart; iLocal < iHaloStart + iHaloWidth; iLocal++)
			for (int jLocal = jHaloStart; jLocal < jHaloStart + jHaloWidth; jLocal++)
			for (int kLocal = kHaloStart; kLocal < kHaloStart + kHaloWidth; kLocal++) {
				int localMassPos = mortonIndex(iLocal, jLocal, kLocal);
				int iStep = 0, jStep = 0, kStep = 0;
				if ((localMassPos & iMask) < i_update_low_zmix ||
					(localMassPos & jMask) < j_update_low_zmix ||
					(localMassPos & kMask) < k_update_low_zmix ||
					(localMassPos & iMask) >= i_update_high_zmix ||
					(localMassPos & jMask) >= j_update_high_zmix ||
					(localMassPos & kMask) >= k_update_high_zmix ) {
					continue;
				}

				if (zID & iMask){
					iPart[iStep++] = ((zID & iMask)  - 1)  & iMask;
				}
				iPart[iStep++] = zID & iMask;
				if ((zID & iMask) != i_bound_inter){
					 iPart[iStep++] = ((zID | (~iMask)) + 1) & iMask;
				}

				if (zID & jMask){
					jPart[jStep++] = ((zID & jMask)  - 1) & jMask;
				}
				jPart[jStep++] = zID & jMask;
				if ((zID & jMask) != j_bound_inter){
					jPart[jStep++] = ((zID | (~jMask)) + 1) & jMask;
				}

				if (zID  & kMask){
					kPart[kStep++] = ((zID & kMask)  - 1) & kMask;
				}
				kPart[kStep++] = zID & kMask;
				if ((zID & kMask) != k_bound_inter){
				   kPart[kStep++] = ((zID | (~kMask)) + 1) & kMask;
				}

			   int t = 0;
			   for (int i = 0; i < iStep; i++)
			   for (int j = 0; j < jStep; j++)
			   for (int k = 0; k < kStep; k++) {
					int nzID = iPart[i] | jPart[j] | kPart[k];
					if (nzID != zID)
						adjIDArray[t++] = nzID;
				}

				double val = 0;
				for (int s = 0; s < t; s++) {
						val += x_block[adjIDArray[s]];
				}

				val = (b_block[zID] - val) / (-26.0);
				normChange += (val - x_block[zID]) * (val - x_block[zID]);
				x_block[zID] = val;
			}
		}
		normChange = sqrt(normChange) / (iExtBlockSize * jExtBlockSize * kExtBlockSize);
		return normChange;
	}

	/*
	 * di is either 1 or -1, the direction
	 * w is the non-update width
	 */
	void GS_block_interior_update_009 (int di, int w) {
		int start, end;
		if (di == 1) {
			start = 0;
			end   = x_block.size();
		}
		else {
			start =  x_block.size() - 1;
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
				(zID & iMask) >= iHighMorton ||
				(zID & jMask) >= jHighMorton ||
				(zID & kMask) >= kHighMorton ) {
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
			if ((zID & jMask) != j_bound_inter){
				jPart[jStep++] = ((zID | (~jMask)) + 1) & jMask;
			}

			if (zID  & kMask){
				kPart[kStep++] = ((zID & kMask)  - 1) & kMask;
			}
			kPart[kStep++] = zID & kMask;
			if ((zID & kMask) != k_bound_inter){
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
					val += x_block[adjIDArray[s]];
			}

			val = (b_block[zID] - val) / (-26.0);
			normChange += (val - x_block[zID]) * (val - x_block[zID]);
			x_block[zID] = val;
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
		for (int dk = -1; dk <= 1; dk++) if (di | dj | dk) {
			int haloLinearID = (1+di) * 9 + (1+dj) * 3 + (1+dk);

			int iHaloStart = (di == -1) ? 0 : (di == 0) ? halo_width : iBlockSize + halo_width;
			int jHaloStart = (dj == -1) ? 0 : (dj == 0) ? halo_width : jBlockSize + halo_width;
			int kHaloStart = (dk == -1) ? 0 : (dk == 0) ? halo_width : kBlockSize + halo_width;

			int iHaloWidth = di ? halo_width : iBlockSize;
			int jHaloWidth = dj ? halo_width : jBlockSize;
			int kHaloWidth = dk ? halo_width : kBlockSize;

			for (int iLocal = iHaloStart; iLocal < iHaloStart + iHaloWidth; iLocal++)
			for (int jLocal = jHaloStart; jLocal < jHaloStart + jHaloWidth; jLocal++)
			for (int kLocal = kHaloStart; kLocal < kHaloStart + kHaloWidth; kLocal++) {
				int localMassPos = mortonIndex(iLocal, jLocal, kLocal);
				int innerHaloPos =
						  (iLocal - iHaloStart) * jHaloWidth * kHaloWidth
						+ (jLocal - jHaloStart) * kHaloWidth
						+ (kLocal - kHaloStart);
				x_block[localMassPos] = recvHalo[haloLinearID][innerHaloPos];
			}
		}
	}

	void sendBackToRoot () {
		MPI::COMM_WORLD.Send (&x_block[0], blockSize , MPI_DOUBLE, root_process, send_data_tag);
		MPI::COMM_WORLD.Send (&this->normChange, 1 , MPI_DOUBLE, root_process, send_data_tag);
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
		for (int dk = -1; dk <= 1; dk++) if (di | dj | dk) {
				int haloLinearID = (1+di)*9 + (1+dj)*3 + (1+dk);
				int iHaloWidth = di ? halo_width : iBlockSize;
				int jHaloWidth = dj ? halo_width : jBlockSize;
				int kHaloWidth = dk ? halo_width : kBlockSize;
				int haloBlockSize = haloWidth_i * haloWidth_j * haloWidth_k;
				vacantHalo[haloLinearID] = vector<double> (haloBlockSize, 0);
			}
		return vacantHalo;
	}

};

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
	string dataSheetFileName = string(argv[5]);
	int iterationSteps = atoi(argv[6]);

	int num_procs_used = 1 + numCore_i * numCore_j * numCore_k;  // one for monitor and root process
	if (num_procs < num_procs_used) {
		cerr << "Not Enough Processes. " << endl;
		return 108;
	}
	int root_process = num_procs_used - 1;

	fstream dataSheet (dataSheetFileName);
	dataSheet >> Capsule::i_bound >> Capsule::j_bound >> Capsule::k_bound;

	Capsule::initCommonStruct();
	Capsule::prepareIndexMapping();

	int totalSize = Capsule::totalSize;

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
		printf("*#06#, %f\t%f\t%d\n", max_elapsed_secs, normChange, iterationSteps);
//		double err = 0;
//		for (double val : x) {
//			err += (val-1.0)*(val-1.0);
//		}
//		err = sqrt(err);
//		printf("Total Error: %f\n", err);
		dataSheet.close();
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

			capsule.prepareSendHalo();
			capsule.haloExchange(); 
			
//			err += capsule.GS_Z_block_update007_interior(i & 1 ?  - 1 : 1); // fwd/bwd
			err += capsule.GS_Z_block_update007_interior(1); // fwd
			
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
