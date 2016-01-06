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

class Index {
public:
	int i, j, k;
};

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

double GS_Z_update (
	vector<double> &x, const vector<double> &b, 
	int iBlockSize, int jBlockSize, int kBlockSize, int blockBitWidth ) {
	
	double avgNormChange = 0;
	for (int zID = 0; zID < x.size(); zID++) {
		vector<int> regID = zDeMix(zID, blockBitWidth);
		int i_sub = regID[0], j_sub = regID[1], k_sub = regID[2];
		if (i_sub >= iBlockSize || j_sub >= jBlockSize || k_sub >= kBlockSize)
			continue;

		double val = 0;
		for (int di = -1; di <= 1; di++) 
		for (int dj = -1; dj <= 1; dj++)
		for (int dk = -1; dk <= 1; dk++) 
			if ((di | dj | dk) && 
				i_sub+di >= 0 && i_sub+di < iBlockSize &&  
				j_sub+dj >= 0 && j_sub+dj < jBlockSize &&  
				k_sub+dk >= 0 && k_sub+dk < kBlockSize ){
				int zAdjID = zmix(i_sub+di, j_sub+dj, k_sub+dk, blockBitWidth);
				val += x[zAdjID]; // add coefficent later
			}
		val = (b[zID] - val) / (-26.0);
		avgNormChange += (val - x[zID]) * (val - x[zID]);
		x[zID] = val;
	}
	avgNormChange = sqrt(avgNormChange) / (iBlockSize * jBlockSize * kBlockSize); 
	return avgNormChange;
}

double GS_Z_local_update (
	vector<double> &x, const vector<double> &b, 
	int i_bound, int j_bound, int k_bound, 
	int iStart, int jStart, int kStart, 
	int iBlockSize, int jBlockSize, int kBlockSize, 
	int halo_width, int blockBitWidth) {
	
	double normChange = 0;
	for (int zID = 0; zID < x.size(); zID++) {
		vector<int> regID = zDeMix(zID, blockBitWidth);
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
				int zAdjID = zmix(i_sub+di, j_sub+dj, k_sub+dk, blockBitWidth);
				val += x[zAdjID]; // add coefficent later
			}
		val = (b[zID] - val) / (-26.0);
		// only count inner part change total.
		if (i_sub >= halo_width && j_sub >= halo_width && k_sub >= halo_width && 
			i_sub < iBlockSize && j_sub < jBlockSize && k_sub < kBlockSize)
			normChange += (val - x[zID]) * (val - x[zID]);
		x[zID] = val;
	}
	normChange /= (iBlockSize * jBlockSize * kBlockSize);
	return normChange;
}


// inner 
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

//outer 
void updateFromHalo (
	vector<vector<double> > &recvHalo, vector<double> &x, 
	int iBlockSize, int jBlockSize, int kBlockSize, int halo_width, int blockBitWidth) {
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
				int zID = zmix(i_sub, j_sub, k_sub, blockBitWidth);
				int subHaloID = (i_sub - iHaloStart) * jHaloWidth * kHaloWidth + (j_sub - jHaloStart) * kHaloWidth + (k_sub - kHaloStart);
				x[zID] = recvHalo[haloLinearID][subHaloID];
			}
	}
}


int main (int argc, char *argv[]) {
	if (argc < 6) {
		cerr << "Please input : numCore_i, numCore_j, numCore_k, halo_width, inputFileName" << endl;
		return 107;
	}

	MPI_Status status;
	int my_id, num_procs, an_id, sender;

	int ierr = MPI_Init(&argc, &argv);
	ierr = MPI_Comm_rank(MPI_COMM_WORLD, &my_id);
	ierr = MPI_Comm_size(MPI_COMM_WORLD, &num_procs);

	int numCore_i  = atoi(argv[1]);
	int numCore_j  = atoi(argv[2]);
	int numCore_k  = atoi(argv[3]);
	int halo_width = atoi(argv[4]);
	string dataSheetFileName = string(argv[5]);

	int num_procs_used = 1 + numCore_i * numCore_j * numCore_k;  // one for monitor and root process

	if (num_procs < num_procs_used) {
		cerr << "Not Enough Proccesses. " << endl;
		return 108;
	}

	int root_process = num_procs_used - 1;

	int i_bound = -1, j_bound = -1, k_bound = -1, totalSize = -1;
	if (my_id == root_process) {
		fstream dataSheet (dataSheetFileName);
		dataSheet >> i_bound >> j_bound >> k_bound;

		totalSize = i_bound * j_bound * k_bound;
		for (int pid = 0; pid < root_process; pid++) {
			ierr = MPI_Send(&i_bound, 1 , MPI_INT, pid, send_data_tag, MPI_COMM_WORLD);
			ierr = MPI_Send(&j_bound, 1 , MPI_INT, pid, send_data_tag, MPI_COMM_WORLD);
			ierr = MPI_Send(&k_bound, 1 , MPI_INT, pid, send_data_tag, MPI_COMM_WORLD);
			ierr = MPI_Send(&totalSize, 1 , MPI_INT, pid, send_data_tag, MPI_COMM_WORLD);
		}

		vector<double> x (totalSize, 1.0);
		for (int i = 0; i < x.size(); i++) {
			x[i] = 2 + 0.001 * i + 0.00007;
		}

		vector<double> b (totalSize);
		for (double &entry : b) {
			dataSheet >> entry;
		}


		int iBlockSize = (i_bound % numCore_i) ? i_bound / numCore_i + 1: i_bound / numCore_i;
		int jBlockSize = (j_bound % numCore_j) ? j_bound / numCore_j + 1: j_bound / numCore_j;
		int kBlockSize = (k_bound % numCore_k) ? k_bound / numCore_k + 1: k_bound / numCore_k;

		int iExtBlockSize = iBlockSize + 2 * halo_width;
		int jExtBlockSize = jBlockSize + 2 * halo_width;
		int kExtBlockSize = kBlockSize + 2 * halo_width;
		
		for (int core_i = 0; core_i < numCore_i; core_i++)
		for (int core_j = 0; core_j < numCore_j; core_j++)
		for (int core_k = 0; core_k < numCore_k; core_k++) {
			int iStart = core_i * iBlockSize - halo_width;
			int jStart = core_j * jBlockSize - halo_width;
			int kStart = core_k * kBlockSize - halo_width;

			int blockBitWidth = ceil(max3(log2(iExtBlockSize), log2(jExtBlockSize), log2(kExtBlockSize)));
			int blockSize = 1 << (3 * blockBitWidth); 
			
			vector<double> x_block (blockSize);
			vector<double> b_block (blockSize);

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
				int targetPos = zmix(i_sub, j_sub, k_sub, blockBitWidth);
				x_block[targetPos] = x[srcPos];
				b_block[targetPos] = b[srcPos];
			}

			// send
			int core_linear = core_i * numCore_j * numCore_k + core_j * numCore_k + core_k;
			ierr = MPI_Send(&core_i, 1 , MPI_INT, core_linear, send_data_tag, MPI_COMM_WORLD);
			ierr = MPI_Send(&core_j, 1 , MPI_INT, core_linear, send_data_tag, MPI_COMM_WORLD);
			ierr = MPI_Send(&core_k, 1 , MPI_INT, core_linear, send_data_tag, MPI_COMM_WORLD);

			ierr = MPI_Send(&iStart, 1 , MPI_INT, core_linear, send_data_tag, MPI_COMM_WORLD);
			ierr = MPI_Send(&jStart, 1 , MPI_INT, core_linear, send_data_tag, MPI_COMM_WORLD);
			ierr = MPI_Send(&kStart, 1 , MPI_INT, core_linear, send_data_tag, MPI_COMM_WORLD);

			ierr = MPI_Send(&iExtBlockSize, 1 , MPI_INT, core_linear, send_data_tag, MPI_COMM_WORLD);
			ierr = MPI_Send(&jExtBlockSize, 1 , MPI_INT, core_linear, send_data_tag, MPI_COMM_WORLD);
			ierr = MPI_Send(&kExtBlockSize, 1 , MPI_INT, core_linear, send_data_tag, MPI_COMM_WORLD);

			ierr = MPI_Send(&blockBitWidth, 1 , MPI_INT, core_linear, send_data_tag, MPI_COMM_WORLD);
			ierr = MPI_Send(&blockSize, 1 , MPI_INT, core_linear, send_data_tag, MPI_COMM_WORLD);

			ierr = MPI_Send(&x_block[0], blockSize , MPI_DOUBLE, core_linear, send_data_tag, MPI_COMM_WORLD);
			ierr = MPI_Send(&b_block[0], blockSize , MPI_DOUBLE, core_linear, send_data_tag, MPI_COMM_WORLD);

		}

		// collect 
		
		for (int core_i = 0; core_i < numCore_i; core_i++)
		for (int core_j = 0; core_j < numCore_j; core_j++)
		for (int core_k = 0; core_k < numCore_k; core_k++) {
			int iStart = core_i * iBlockSize - halo_width;
			int jStart = core_j * jBlockSize - halo_width;
			int kStart = core_k * kBlockSize - halo_width;

			int blockBitWidth = ceil(max3(log2(iExtBlockSize), log2(jExtBlockSize), log2(kExtBlockSize)));
			int blockSize = 1 << (3 * blockBitWidth); 
			
			vector<double> x_block (blockSize);
			int core_linear = core_i * numCore_j * numCore_k + core_j * numCore_k + core_k;
			// ierr = MPI_Recv(&x_block[0], blockSize, MPI_DOUBLE, core_linear, send_data_tag, MPI_COMM_WORLD, &status);

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
				int targetPos = zmix(i_sub, j_sub, k_sub, blockBitWidth);
				x[srcPos] = x_block[targetPos];
			}

		}
		for (double val : x) {
			printf("%f\n", val);
		}
	}
	else if (my_id < root_process){
		ierr = MPI_Recv(&i_bound, 1, MPI_INT, root_process, send_data_tag, MPI_COMM_WORLD, &status);
		ierr = MPI_Recv(&j_bound, 1, MPI_INT, root_process, send_data_tag, MPI_COMM_WORLD, &status);
		ierr = MPI_Recv(&k_bound, 1, MPI_INT, root_process, send_data_tag, MPI_COMM_WORLD, &status);
		ierr = MPI_Recv(&totalSize, 1, MPI_INT, root_process, send_data_tag, MPI_COMM_WORLD, &status);

		int core_i, core_j, core_k; 
		ierr = MPI_Recv(&core_i, 1, MPI_INT, root_process, send_data_tag, MPI_COMM_WORLD, &status);
		ierr = MPI_Recv(&core_j, 1, MPI_INT, root_process, send_data_tag, MPI_COMM_WORLD, &status);
		ierr = MPI_Recv(&core_k, 1, MPI_INT, root_process, send_data_tag, MPI_COMM_WORLD, &status);

		int iStart, jStart, kStart; 
		ierr = MPI_Recv(&iStart, 1, MPI_INT, root_process, send_data_tag, MPI_COMM_WORLD, &status);
		ierr = MPI_Recv(&jStart, 1, MPI_INT, root_process, send_data_tag, MPI_COMM_WORLD, &status);
		ierr = MPI_Recv(&kStart, 1, MPI_INT, root_process, send_data_tag, MPI_COMM_WORLD, &status);

		int iExtBlockSize, jExtBlockSize, kExtBlockSize;
		ierr = MPI_Recv(&iExtBlockSize, 1, MPI_INT, root_process, send_data_tag, MPI_COMM_WORLD, &status);
		ierr = MPI_Recv(&jExtBlockSize, 1, MPI_INT, root_process, send_data_tag, MPI_COMM_WORLD, &status);
		ierr = MPI_Recv(&kExtBlockSize, 1, MPI_INT, root_process, send_data_tag, MPI_COMM_WORLD, &status);
		// printf("core: %i:%i %i \n",iStart, jStart, kStart);

		int iBlockSize = iExtBlockSize - halo_width;
		int jBlockSize = iExtBlockSize - halo_width;
		int kBlockSize = iExtBlockSize - halo_width;

		int blockBitWidth, blockSize;
		ierr = MPI_Recv(&blockBitWidth, 1, MPI_INT, root_process, send_data_tag, MPI_COMM_WORLD, &status);
		ierr = MPI_Recv(&blockSize, 1, MPI_INT, root_process, send_data_tag, MPI_COMM_WORLD, &status);

		vector<double> x_block (blockSize);
		vector<double> b_block (blockSize);
		ierr = MPI_Recv(&x_block[0], blockSize, MPI_DOUBLE, root_process, send_data_tag, MPI_COMM_WORLD, &status);
		ierr = MPI_Recv(&b_block[0], blockSize, MPI_DOUBLE, root_process, send_data_tag, MPI_COMM_WORLD, &status);

		vector<vector<double> > sendHalo = haloInit(iBlockSize, jBlockSize, kBlockSize, halo_width); // resize
		vector<vector<double> > recvHalo = haloInit(iBlockSize, jBlockSize, kBlockSize, halo_width); 
		
		double err = 1;
		int n = 10;
		while (n--) {
			double err = GS_Z_local_update (
				x_block, b_block, 
				i_bound, j_bound, k_bound, 
				iStart, jStart, kStart,
				 iBlockSize, jBlockSize, kBlockSize, 
				 halo_width, blockBitWidth);

			prepareSendHalo (sendHalo, x_block,  iBlockSize, jBlockSize, kBlockSize, halo_width, blockBitWidth);
			
			// send
			for (int di = -1; di <= 1; di++) 
			for (int dj = -1; dj <= 1; dj++)
			for (int dk = -1; dk <= 1; dk++) 
				if ((di | dj | dk) && 
					core_i + di >= 0 && core_i + di < numCore_i &&
					core_j + dj >= 0 && core_j + dj < numCore_j &&
					core_k + dk >= 0 && core_k + dk < numCore_k ) {
				int destCoreLinear = (core_i + di) * numCore_j * numCore_k + (core_j + dj) * numCore_k + (core_k + dk);
				int srcHaloID =  (1+di)*9 + (1+dj)*3 + (1+dk);
				int destHaloID = (1-di)*9 + (1-dj)*3 + (1-dk); // subject to verification
				int haloSize = sendHalo[srcHaloID].size();
				
				ierr = MPI_Send (&destHaloID, 1 , MPI_INT, destCoreLinear, send_data_tag, MPI_COMM_WORLD);
				ierr = MPI_Send (&haloSize, 1 ,   MPI_INT, destCoreLinear, send_data_tag, MPI_COMM_WORLD);
				ierr = MPI_Send (&sendHalo[srcHaloID][0], haloSize , MPI_DOUBLE, destCoreLinear, send_data_tag, MPI_COMM_WORLD);	
			}
			
			//receive 

			for (int di = -1; di <= 1; di++) 
			for (int dj = -1; dj <= 1; dj++)
			for (int dk = -1; dk <= 1; dk++) 
				if ((di | dj | dk) && 
					core_i + di >= 0 && core_i + di < numCore_i &&
					core_j + dj >= 0 && core_j + dj < numCore_j &&
					core_k + dk >= 0 && core_k + dk < numCore_k ) {
				int srcCoreLinear = (core_i + di) * numCore_j * numCore_k + (core_j + dj) * numCore_k + (core_k + dk);
				int srcHaloID;// = (1-di)*9 + (1-dj)*3 + (1-dk); 
				int destHaloID;//  =  (1+di)*9 + (1+dj)*3 + (1+dk); 
				int haloSize;

				ierr = MPI_Recv(&destHaloID, 1, MPI_INT, srcCoreLinear, send_data_tag, MPI_COMM_WORLD, &status);
				ierr = MPI_Recv(&haloSize,   1, MPI_INT, srcCoreLinear, send_data_tag, MPI_COMM_WORLD, &status);
				ierr = MPI_Recv(&recvHalo[destHaloID][0], haloSize, MPI_DOUBLE, srcCoreLinear, send_data_tag, MPI_COMM_WORLD, &status);
				
			}
			updateFromHalo (recvHalo, x_block, iBlockSize, jBlockSize, kBlockSize, halo_width, blockBitWidth);
		}

		// send back to root process 
		ierr = MPI_Send (&x_block[0], blockSize , MPI_DOUBLE, root_process, send_data_tag, MPI_COMM_WORLD);


	}

	ierr = MPI_Finalize();

	return 0;
}