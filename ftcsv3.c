#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

int main(int argc, char *argv[])
{
    // MPI things
    int npes, myrank;
    MPI_Status status;
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &npes);
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

    // Initialization
    double D = 0.1;
    double dx = 0.15;
    double dt = 0.1;
    double L = 1.0;
    double T = 1.0;
    double r = D * dt / pow(dx, 2);
    int Text = 300;
    int Tint = 100;
    int ntime = T / dt;
    int nspace = L / dx + 1;
    int root;
    if (myrank == 0)
    {
        root = 1;
    }
    else
        root = 0;

    if (root)
    {
        printf("nspace : %d, ntime : %d\n", nspace, ntime);
    }

    double results[ntime][nspace];

    // Initialize grid
    for (int n = 0; n < ntime; n++)
        for (int i = 0; i < nspace; i++)
            results[n][i] = 0.0;
    for (int i = 0; i < nspace; i++)
        results[0][i] = Tint;
    for (int n = 0; n < ntime; n++)
    {
        results[n][0] = Text;
        results[n][nspace - 1] = Text;
    }
    results[1][0] = Text;
    results[1][nspace - 1] = Text;
    //- Initialisation finished

    if (npes > nspace)
    {
        if (root)
            printf("Too much processors for considered problem.\n");
    }
    else if (npes == 1)
    {
        if (root)
            for (int n = 1; n < ntime; n++)
                for (int i = 1; i < nspace - 1; i++)
                    results[n][i] = results[n - 1][i] + r * (results[n - 1][i + 1] - 2 * results[n - 1][i] + results[n - 1][i - 1]);
        // Prints solution matrix
        if (root)
        {
            for (int i = 0; i < ntime; i++)
            {
                printf("rank : %d, line %d : ", myrank, i + 1);
                for (int j = 0; j < nspace; j++)
                {
                    printf("%3.2f ", results[i][j]);
                }
                printf("\n");
            }
        }
    }
    // case multi processors
    else
    {
        if (nspace % npes != 0)
        {
            if (root)
            {
                printf("Wrong number of processors (%d) for the problem size (%d)\n", npes, nspace);
                printf("Try a number of processors which is able to divide the problem size.\n");
            }
        }
        else
        {
            int newspace = nspace / npes;
        }
    }

    if (root)
        for (int i = 0; i < ntime; i++)
        {
            printf("rank : %d, line %d : ", myrank, i + 1);
            for (int j = 0; j < nspace; j++)
            {
                printf("%3.2f ", results[i][j]);
            }
            printf("\n");
        }
    MPI_Finalize();
}