#include <mpi.h>
#include <stdio.h>

int main(int argc, char *argv[])
{
    int npes, myrank;
    MPI_Status status;
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &npes);
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

    printf("From process %d out of %d, Hello world!\n", myrank, npes);

    MPI_Finalize();
}
