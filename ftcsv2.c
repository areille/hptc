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
    double dx = 0.5;
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
    else
    {
        // case 3 processors 3 columns
        if (nspace % npes == 0)
        {
            int newspace = nspace / npes;
            double p_results[ntime][newspace];
            double boundary_mine, boundary_rcv;
            printf("newspace : %d\n", newspace);
            // INITIALISATION : FILL MATRIX WITH 0 AND BOUNDARY CONDITIONS
            for (int i = 0; i < ntime; i++)
            {
                for (int j = 0; j < newspace; j++)
                {
                    p_results[i][j] = results[i][j + (myrank * newspace)];
                }
            }
            // FILL UP LINE 1 to NTIME
            for (int n = 1; n < ntime; n++)
            {
                // FILL UP COLUMN 1 FIRST THEN THE OTHERS UNTIL MISSING VALUE
                for (int i = 1; i < newspace - 1; i++)
                {
                    p_results[n][i] = p_results[n - 1][i] + r * (p_results[n - 1][i + 1] - 2 * p_results[n - 1][i] + p_results[n - 1][i - 1]);
                }
                for (int prl = 0; prl < (npes - 1) / 2; prl++)
                {
                    if (myrank == prl)
                    {
                        printf("prl : %d\n", prl);
                        boundary_mine = p_results[n - 1][newspace - 1];
                        if (myrank < npes - 1)
                            MPI_Send(&boundary_mine, 1, MPI_DOUBLE, myrank + 1, 0, MPI_COMM_WORLD);
                        if (myrank > 1)
                        {
                            MPI_Recv(&boundary_rcv, 1, MPI_DOUBLE, myrank - 1, 0, MPI_COMM_WORLD, &status);
                            p_results[n][newspace - 1] = p_results[n - 1][newspace - 1] + r * (boundary_rcv - 2 * p_results[n - 1][newspace - 1] + p_results[n - 1][newspace - 2]);
                        }
                    }
                    if (myrank == prl + 1)
                    {
                        MPI_Recv(&boundary_rcv, 1, MPI_DOUBLE, myrank - 1, 0, MPI_COMM_WORLD, &status);
                        p_results[n][newspace - 1] = p_results[n - 1][newspace - 1] + r * (boundary_rcv - 2 * p_results[n - 1][newspace - 1] + p_results[n - 1][newspace - 2]);
                    }
                }
                for (int prr = npes - 1; prr > (npes - 1) / 2; prr--)
                {
                    if (myrank == prr)
                    {
                        printf("prr : %d\n", prr);
                        boundary_mine = p_results[n - 1][0];
                        if (myrank > 1)
                            MPI_Send(&boundary_mine, 1, MPI_DOUBLE, myrank - 1, 1, MPI_COMM_WORLD);
                        if (myrank < npes - 1)
                        {
                            MPI_Recv(&boundary_rcv, 1, MPI_DOUBLE, myrank + 1, 1, MPI_COMM_WORLD, &status);
                            p_results[n][0] = p_results[n - 1][0] + r * (p_results[n - 1][1] - 2 * p_results[n - 1][0] + boundary_rcv);
                        }
                    }
                    if (myrank == prr - 1)
                    {
                        MPI_Recv(&boundary_rcv, 1, MPI_DOUBLE, myrank + 1, 1, MPI_COMM_WORLD, &status);
                        p_results[n][0] = p_results[n - 1][0] + r * (p_results[n - 1][1] - 2 * p_results[n - 1][0] + boundary_rcv);
                    }
                }
            }
            // OUTPUT: PRINTS PARTICULAR RESULTS
            for (int i = 0; i < ntime; i++)
            {
                printf("rank : %d, line %d : ", myrank, i + 1);
                for (int j = 0; j < newspace; j++)
                {
                    printf("%3.2f ", p_results[i][j]);
                }
                printf("\n");
            }
        }
        else
        {
        }
    }

    MPI_Finalize();
}