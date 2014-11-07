#include <iostream>
#include <cstdio>
#include <vector>
#include <mpi.h>

using namespace std;

int main(int argc, char *argv[]) {
  int  len, rc;
  char hostname[MPI_MAX_PROCESSOR_NAME];

  MPI::Status status;

  MPI::Init(argc, argv);
  if (rc != MPI_SUCCESS) {
    printf ("Error starting MPI program. Terminating.\n");
    MPI::COMM_WORLD.Abort (rc);
  }

  int numtasks = MPI::COMM_WORLD.Get_size();
  int rank     = MPI::COMM_WORLD.Get_rank();

  MPI_Get_processor_name(hostname, &len);

  vector<int> x = {2 * rank, 2 * rank + 1}; 
  vector<int> y = vector<int> (10);

  vector<int> a (1000); 
  if (rank == 0) {
      a = vector<int> ({1,2,3,4,5,6,7,8,9});   
  }
  
  MPI::COMM_WORLD.Bcast (&a[0], 3, MPI_INT, 0);

  MPI::COMM_WORLD.Gather(&x[0], 2, MPI::INT, &y[0], 2, MPI::INT, 4);

  printf ("Number of tasks= %d My rank= %d Running on %s\n", numtasks,rank, hostname);
  printf("y : %i, %i %i\n", rank, y[0], y[1]);
  /*******  do some work *******/
  if (rank == 4) {
    for (auto val : y)
      printf(":: %i :: %i\n", rank, val);

  }

  
  if (rank == 1){
    int x = 10, y = 20;
    MPI::COMM_WORLD.Send(&x, 1, MPI_INT, 2, 127);
    MPI::COMM_WORLD.Send(&y, 1, MPI_INT, 2, 128);  
  }
  else if (rank == 2) {
    int x = -10, y = -20;
    MPI::COMM_WORLD.Recv(&y, 1, MPI_INT, 1, 128);
    MPI::COMM_WORLD.Recv(&x, 1, MPI_INT, 1, 127);
    printf("---->2:::::%i %i\n", x, y );
  }

  MPI::Finalize();
  return 0;
}











