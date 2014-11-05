// halo.cpp

// read regular
// distribute
// Z-order load
// Z-order send
// tags

#include <iostream>
#include <vector>
#include <string>
#include <cmath>
#include <fstream>
#include <mpi.h>

using namespace std;

#define send_data_tag 2001
#define return_data_tag 2002

#define max3(x, y, z) max(max((x), (y)), (z))

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


class Capsule {
public:
	static MPI_Status status;
	static int i_bound, j_bound, k_bound, totalSize;
	static int numCore_i, numCore_j, numCore_k, coresOnDuty, root_process;
	static int halo_width;

	int core_i, core_j, core_k, core_linear;
	int iStart, jStart, kStart; // as an extended

	int iBlockSize, jBlockSize, kBlockSize;

	int iExtBlockSize, jExtBlockSize, kExtBlockSize;
	int blockBitWidth, blockSize;

	vector<double> x_block;
	vector<double> b_block;

	vector<vector<double> > sendHalo; // x
	vector<vector<double> > recvHalo; // x

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

		x_block = vector<double> (blockSize, 0);
		b_block = vector<double> (blockSize, 0);

		sendHalo = haloInit ();
		recvHalo = haloInit ();

	}

	int zmix (int i, int j, int k) {
	    int zVal = 0;
	    for (int mask = (1 << (blockBitWidth-1)); mask; mask >>= 1) {
	        zVal = (i & mask) ? (zVal<<1) + 1 : zVal<<1;
	        zVal = (j & mask) ? (zVal<<1) + 1 : zVal<<1;
	        zVal = (k & mask) ? (zVal<<1) + 1 : zVal<<1;
	    }
		return zVal;
	}

	vector<int> zDeMix(int zVal) {
		int i = 0, j = 0, k = 0;
		for (int mask = (1 << (3*blockBitWidth-1)); mask; mask >>= 3) {
	        i = (zVal &  mask)     ? (i<<1) + 1 : i<<1;
	        j = (zVal & (mask>>1)) ? (j<<1) + 1 : j<<1;
	        k = (zVal & (mask>>2)) ? (k<<1) + 1 : k<<1;
	    }
		return {i, j, k};
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

	double GS_Z_block_update () {
		double normChange = 0;
		for (int zID = 0; zID < x_block.size(); zID++) {
			vector<int> regID = zDeMix(zID);
			int i_sub = regID[0], j_sub = regID[1], k_sub = regID[2];
			if (i_sub >= iBlockSize + halo_width || j_sub >= jBlockSize + halo_width || k_sub >= kBlockSize + halo_width|| 
				iStart + i_sub < 0         || jStart + j_sub <  0        || kStart + k_sub <  0       || 
				iStart + i_sub >= i_bound  || jStart + j_sub >= j_bound  || kStart + k_sub >= k_bound )
				continue;

			double val = 0;
			for (int di = -1; di <= 1; di++) 
			for (int dj = -1; dj <= 1; dj++)
			for (int dk = -1; dk <= 1; dk++) 
				if ((di | dj | dk) && 
					i_sub+di >= 0 && i_sub+di < iBlockSize + halo_width &&  
					j_sub+dj >= 0 && j_sub+dj < jBlockSize + halo_width &&  
					k_sub+dk >= 0 && k_sub+dk < kBlockSize + halo_width ){
					int zAdjID = zmix(i_sub+di, j_sub+dj, k_sub+dk);
					val += x_block[zAdjID]; // add coefficent later
				}
			val = (b_block[zID] - val) / (-26.0);
			// only count inner part change total.
			if (i_sub >= halo_width && j_sub >= halo_width && k_sub >= halo_width && 
				i_sub < iBlockSize && j_sub < jBlockSize && k_sub < kBlockSize)
				normChange += (val - x_block[zID]) * (val - x_block[zID]);
			x_block[zID] = val;
		}
		normChange /= (iBlockSize * jBlockSize * kBlockSize);
		return normChange;
	}

	int haloSend () {
		int ierr;
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
			if (haloSize) {
				ierr = MPI_Send (&sendHalo[srcHaloID][0], haloSize , MPI_DOUBLE, destCoreLinear, send_data_tag, MPI_COMM_WORLD);
			}
			
		}
		return ierr;
	}

	int haloRecv() {
		int ierr;
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
				ierr = MPI_Recv(&recvHalo[destHaloID][0], haloSize, MPI_DOUBLE, srcCoreLinear, send_data_tag, MPI_COMM_WORLD, &status);
			}
		}
		return ierr;
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

	int sendBackToRoot () {
		return MPI_Send (&x_block[0], blockSize , MPI_DOUBLE, root_process, send_data_tag, MPI_COMM_WORLD);
	}

	void syncBackAtRoot(vector<double> &x) {
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
			x[srcPos] = x_block[targetPos];
		}
	}

	void betaSend() {
		int n = 1, my_id;
		int cores = numCore_i * numCore_j * numCore_k;
		// MPI_Send (&n, 1 , MPI_INT, (core_linear+1) % cores, send_data_tag, MPI_COMM_WORLD);
		MPI_Send (&n, 1 , MPI_INT, (core_linear + 1) % cores, send_data_tag, MPI_COMM_WORLD);

	}


} ;

MPI_Status Capsule::status;
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
	if (argc < 6) {
		cerr << "Please input : numCore_i, numCore_j, numCore_k, halo_width, inputFileName" << endl;
		return 107;
	}

	MPI_Status status;
	int my_id, num_procs, an_id, sender, ierr;

	ierr = MPI_Init(&argc, &argv);
	ierr = MPI_Comm_rank(MPI_COMM_WORLD, &my_id);
	ierr = MPI_Comm_size(MPI_COMM_WORLD, &num_procs);

	int numCore_i  = Capsule::numCore_i  = atoi(argv[1]);
	int numCore_j  = Capsule::numCore_j  = atoi(argv[2]);
	int numCore_k  = Capsule::numCore_k  = atoi(argv[3]);
	int halo_width = Capsule::halo_width = atoi(argv[4]);
	string dataSheetFileName = string(argv[5]);

	int num_procs_used = 1 + numCore_i * numCore_j * numCore_k;  // one for monitor and root process
	if (num_procs < num_procs_used) {
		cerr << "Not Enough Proccesses. " << endl;
		return 108;
	}
	int root_process = num_procs_used - 1;

	fstream dataSheet (dataSheetFileName);
	dataSheet >> Capsule::i_bound >> Capsule::j_bound >> Capsule::k_bound;
	int totalSize = Capsule::totalSize = Capsule::i_bound * Capsule::j_bound * Capsule::k_bound;
	if (my_id == root_process) {		
		vector<double> x (totalSize, 1.0);
		vector<double> b (totalSize);
		
		for (int i = 0; i < x.size(); i++) {
			x[i] = 2 + 0.001 * i + 0.00007;
		}

		for (double &entry : b) {
			dataSheet >> entry;
		}
		
		for (int core_i = 0; core_i < numCore_i; core_i++)
		for (int core_j = 0; core_j < numCore_j; core_j++)
		for (int core_k = 0; core_k < numCore_k; core_k++) {
			vector<int> core_pack = {core_i, core_j, core_k};
			Capsule mainCull (core_i, core_j, core_k);
			vector<double> x_block = mainCull.extractFrom(x);
			vector<double> b_block = mainCull.extractFrom(b);

			ierr = MPI_Send (&core_pack[0], 3 , MPI_INT, mainCull.core_linear, send_data_tag, MPI_COMM_WORLD);
			ierr = MPI_Send (&x_block[0], x_block.size(), MPI_DOUBLE, mainCull.core_linear, send_data_tag, MPI_COMM_WORLD);
			ierr = MPI_Send (&b_block[0], b_block.size(), MPI_DOUBLE, mainCull.core_linear, send_data_tag, MPI_COMM_WORLD);

		}

		for (int core_i = 0; core_i < numCore_i; core_i++)
		for (int core_j = 0; core_j < numCore_j; core_j++)
		for (int core_k = 0; core_k < numCore_k; core_k++) {
			vector<int> core_pack = {core_i, core_j, core_k};
			Capsule mainCull (core_i, core_j, core_k);
			ierr = MPI_Recv(&mainCull.x_block[0], mainCull.blockSize, MPI_DOUBLE, mainCull.core_linear, send_data_tag, MPI_COMM_WORLD, &status);
			mainCull.syncBackAtRoot(x);
		}
		/*
		for (double val : x) {
			printf("%f\n", val);
		}
		*/


	}
	else if (my_id < root_process) { 
		vector<int> core_pack (3);
		ierr = MPI_Recv(&core_pack[0], 3, MPI_INT, root_process, send_data_tag, MPI_COMM_WORLD, &status);
		Capsule capsule (core_pack[0], core_pack[1], core_pack[2]);
		
		ierr = MPI_Recv(&capsule.x_block[0], capsule.blockSize, MPI_DOUBLE, root_process, send_data_tag, MPI_COMM_WORLD, &status);
		ierr = MPI_Recv(&capsule.b_block[0], capsule.blockSize, MPI_DOUBLE, root_process, send_data_tag, MPI_COMM_WORLD, &status);
		
		int n = 10;
		clock_t t0 = clock();
		while (n--) {
			capsule.GS_Z_block_update();
			capsule.prepareSendHalo();
			capsule.haloSend();
			capsule.haloRecv();
			capsule.updateFromRecvHalo();
		}
		clock_t t1 = clock();
    	double elapsed_secs = double(t1 - t0) / CLOCKS_PER_SEC;
		printf("Core: %02i, seconds_elapsed : %f\n", my_id, elapsed_secs);
		capsule.sendBackToRoot();
		
	}
	ierr = MPI_Finalize();
	return 0;
}