#include <iostream>
#include <mpi.h>

int main(int argc, char** argv)
{
    int npes, myrank;
    MPI_Status status;
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &npes);
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
    std::cout << "Hello World!" << std::endl;

    MPI_Finalize();
    return 0;
}