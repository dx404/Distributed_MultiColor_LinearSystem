// main.cpp
#include <iostream>
#include <vector>
#include <mpi.h>
using namespace std;

#define send_data_tag 2014
#define return_data_tag 2015

int zmix(int i, int j, int k, int bitWidth) {
    int zVal = 0;
    for (int mask = (1 << (bitWidth-1)); mask; mask >>= 1) {
        zVal = (i & mask) ? (zVal<<1) + 1 : zVal<<1;
        zVal = (j & mask) ? (zVal<<1) + 1 : zVal<<1;
        zVal = (k & mask) ? (zVal<<1) + 1 : zVal<<1;
    }
	return zVal;
}

class HPCG_Vector {
public:
    int total_depth;  // steps recursively down to a point, e.g. 3
    int node_depth;   // steps recursively down to a point with in a node, e.g. 2

 	int depth_to_node; // depth_to_node = total_depth - node_depth;

    int total_width;  // number of sample points in each dim
    int total_size;
    int branch_width; // branch factor to the next level of tree

    double convergence_threshold; // 0.001
    
    int node_array_width;
    int node_array_size;

    int node_data_width;
    int node_data_size;

    int halo_width;   // 3

    //vector<double> data; // data to load

    vector<double> coef; // 27
    vector<vector<double> > node_data;  // data passing into each node

    HPCG_Vector(int total_depth, int node_depth, int halo_width = 2, double convergence_threshold = 0.001) {
    	this->total_depth = total_depth; // 4
    	this->node_depth = node_depth;   // 2
        this->branch_width = 2; // 2^3

    	// set up the default coefficients
    	for (int i = 0; i < 27; i++) {
    		coef.push_back(1.0);
    	}
    	coef[13] = -26;

    	// inter-node information
    	depth_to_node = total_depth - node_depth;
    	node_array_width = (1 << depth_to_node);
    	node_array_size = node_array_width * node_array_width * node_array_width; // number of nodes needed

    	// inner-node information
    	node_data_width =  (1 << node_depth) ;
    	node_data_size = node_data_width * node_data_width * node_data_width;
    	node_data =  vector<vector<double> >(node_array_size, vector<double> (node_data_size, 0));

    	total_width = (1 << total_depth) ;
    	total_size = total_width * total_width * total_width;
    }

    inline int getColorID (int i, int j, int k) const {
        return ((i & 1) << 2) + ((j & 1) << 1) + (k & 1);
    }

    inline int getNodeID (int i, int j, int k) const {
        return    ((i >> node_depth) << (depth_to_node * 2))
        		+ ((j >> node_depth) << depth_to_node)
				+  (k >> node_depth);
    }

    inline int getInnerID(int i, int j, int k) const {
    	return zmix(i, j, k, node_depth);
    }

    inline double getElement (int i, int j, int k) const {
    	 return node_data[getNodeID(i, j, k)][getInnerID(i, j, k)];
    }

    inline double setElement (int i, int j, int k, double newVal) {
        return node_data[getNodeID(i, j, k)][getInnerID(i, j, k)] = newVal;
    }

    void read_from_std_to_nodeData() {
        for (int i = 0; i < total_width; i++) {
        	 for (int j = 0; j < total_width; j++) {
        		 for (int k = 0; k < total_width; k++) {
        			 int nodeID = getNodeID(i, j, k);
        			 int innerID = getInnerID(i, j, k);
        			 cin >> node_data[nodeID][innerID];
        		 }
        	 }
        }
    }
};

int main(int argc, char *argv[]) {
	MPI_Status status;
	int ierr = MPI_Init(&argc, &argv);

	/* find out MY process ID, and how many processes were started. */
	int my_id, num_procs;
	int err;
	ierr = MPI_Comm_rank(MPI_COMM_WORLD, &my_id);
	ierr = MPI_Comm_size(MPI_COMM_WORLD, &num_procs);

	int root_process = num_procs - 1; // one more standby process

	int node_data_size;
	int sender;
	if (my_id == root_process) {
		cout << " hello, world" << endl;
		HPCG_Vector *x = new HPCG_Vector(4, 2);
		HPCG_Vector *b = new HPCG_Vector(4, 2);
		// read from IO
		// TO DO

		// Send Message
		for (int nodeID = 0; nodeID < x->node_array_size; nodeID++) {
			ierr = MPI_Send(&x->node_data_size, 1,
					MPI_INT, nodeID, send_data_tag, MPI_COMM_WORLD);
			ierr = MPI_Send(&x->node_data[0][0], x->total_size,
					MPI_DOUBLE, nodeID, send_data_tag, MPI_COMM_WORLD);
			ierr = MPI_Send(&b->node_data[0][0], b->total_size,
					MPI_DOUBLE, nodeID, send_data_tag, MPI_COMM_WORLD);
		}

		for (int nodeID = 0; nodeID < x->node_array_size; nodeID++) {
			ierr = MPI_Recv(&x->node_data[nodeID][0], x->node_data_size,
					MPI_DOUBLE, MPI_ANY_SOURCE, return_data_tag, MPI_COMM_WORLD, &status);
			sender = status.MPI_SOURCE;
		}

		delete x;
		delete b;
	}
	else {
		HPCG_Vector *x = new HPCG_Vector(4, 2);
		HPCG_Vector *b = new HPCG_Vector(4, 2);

		ierr = MPI_Recv(&node_data_size, 1, MPI_INT,
				root_process, send_data_tag, MPI_COMM_WORLD, &status);

		ierr = MPI_Recv(&x->node_data[my_id][0], node_data_size, MPI_DOUBLE,
				root_process, send_data_tag, MPI_COMM_WORLD, &status);

		ierr = MPI_Recv(&b->node_data[my_id][0], node_data_size, MPI_DOUBLE,
				root_process, send_data_tag, MPI_COMM_WORLD, &status);

		// recursively

	}

	ierr = MPI_Finalize();
	return 0;
}
