#include <mpi.h>
#include <stdio.h>

int main(int argc, char** argv) {
  // Initialize the MPI environment
  MPI_Init(&argc, &argv);
  int ierr, thread_id, num_threads;

  // Get the number of processes
  ierr = MPI_Comm_size(MPI_COMM_WORLD, &num_threads);

  if (num_threads != 2) {
    printf("ERROR np nicht durch 2 teilbar\n");
    MPI_Abort(MPI_COMM_WORLD, 1);
  }

  // Get the rank of the process
  ierr = MPI_Comm_rank(MPI_COMM_WORLD, &thread_id);

  // Get the name of the processor
  char processor_name[MPI_MAX_PROCESSOR_NAME];
  int name_len;
  MPI_Get_processor_name(processor_name, &name_len);

  // Print off a hello world message
  printf("Hello world from processor %s, rank %d out of %d processors\n",
         processor_name, thread_id, num_threads);

  char message = 'm';

  if (thread_id == 0) {
    int destination = 1;
    ierr = MPI_Send(&message, 1, MPI_CHAR, destination, 0, MPI_COMM_WORLD);

    if (ierr == MPI_SUCCESS) {
      printf("thread nr %d send message %c to thread nr %d\n", thread_id,
             message, destination);
    } else {
      printf("error beim senden");
    }
  } else {
    char recvmessage;
    MPI_Status status;
    ierr = MPI_Recv(&recvmessage, 1, MPI_CHAR, 0, 0, MPI_COMM_WORLD, &status);

    if (ierr == MPI_SUCCESS) {
      printf("thread nr %d hat message %c empfangen von thread nr 0\n",
             thread_id, recvmessage);
    }
  }
  // Finalize the MPI environment.
  ierr = MPI_Finalize();
}