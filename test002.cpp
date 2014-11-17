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
  
  MPI::COMM_WORLD.Barrier();
  if (rank == 0) {
    printf("=========== %i \n", numtasks);
  }


  if (rank >= 3) {
    int t;
    MPI::COMM_WORLD.Reduce(&rank, &t, 1, MPI::DOUBLE, MPI::SUM, 4);
    printf("ttttt : %i \n", t); 
  }
  else {
    int gg = 0;
    // MPI::COMM_WORLD.Reduce(&gg, &gg, 1, MPI::DOUBLE, MPI::MAX, 4);
    // MPI::COMM_WORLD.Reduce(&gg, &gg, 1, MPI::DOUBLE, MPI::MAX, 2);
  }



  MPI::Finalize();
  return 0;
}











