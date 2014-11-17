#include <iostream>
#include <cstdio>
#include <unistd.h>
#include <vector>
#include <thread>
#include <chrono>
#include <pthread.h>

using namespace std;

#define NUM_THREADS 8
static pthread_mutex_t cs_mutex =  PTHREAD_RECURSIVE_MUTEX_INITIALIZER;


using namespace std;

void *wait (void *t) {
   int i;
   long tid;

   tid = ((long) t);

   sleep(1);
   // this_thread::sleep_for(chrono::milliseconds(10000));
   pthread_mutex_lock( &cs_mutex );
   cout << "Sleeping in thread " << endl;
   cout << "Thread with id : " << tid << "  ...exiting " << endl;
   // printf("%s \n %s %li %s \n", "Sleeping in thread ", "Thread with id : ", tid, "  ...exiting ");
   pthread_mutex_unlock( &cs_mutex );
   pthread_exit(NULL);
}

int main (int argc, char *argv[]) {
   vector<pthread_t> threads (NUM_THREADS);

    // Initialize and set thread joinable
   pthread_attr_t attr;
   pthread_attr_init (&attr);
   pthread_attr_setdetachstate (&attr, PTHREAD_CREATE_JOINABLE);

   void *status;
   for(int i = 0; i < NUM_THREADS; i++ ){
      pthread_mutex_lock( &cs_mutex );
      cout << "main() : creating thread, " << i << endl;
       pthread_mutex_unlock( &cs_mutex );
      int rc = pthread_create(&threads[i], &attr, wait, (void *) i );
      if (rc){
         cout << "Error:unable to create thread," << rc << endl;
         exit(-1);
      }
   }

   // free attribute and wait for the other threads
   pthread_attr_destroy(&attr);
   for(int i=0; i < NUM_THREADS; i++ ){
      int rc = pthread_join(threads[i], &status);
      if (rc){
         cout << "Error:unable to join," << rc << endl;
         exit(-1);
      }
      cout << "Main: completed thread id :" << i ;
      cout << "  exiting with status :" << status << endl;
   }

   cout << "Main: program exiting." << endl;
   pthread_exit(NULL);
}