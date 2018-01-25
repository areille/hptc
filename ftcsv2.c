#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

// TRY TO DO IT WRITING DIRECTLY IN A FILE

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
    double dx = 0.05;
    double dt = 0.01;
    double L = 1.0;
    double T = 1.0;
    double r = D * dt / pow(dx, 2);
    int Text = 300;
    int Tint = 100;
    int tag;
    int ntime = T / dt;
    int nspace = L / dx + 1;
    int root;
    if (myrank == 0)
        root = 1;
    else
        root = 0;

    if (root)
    {
        printf("\n");
        printf("nspace : %d, ntime : %d\n", nspace, ntime);
    }
    // SECURITY : CHECKS THE STABILITY CRITERION VALUE
    if (root)
    {
        printf("\n");
        printf("stability criterion r = %f\n", r);
    }
    if (0.5 <= r)
    {
        if (root)
        {
            printf("\n");
            printf("Exit : r >= 0.5 \n");
            printf("Try to change ∆x and ∆t\n");
        }
        exit(-1);
    }
    else
    {
        if (root)
        {
            printf("\n");
            printf("Continuing : stability criterion ok.\n");
        }
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
        {
            printf("\n");
            printf("Too much processors for considered problem.\n");
            exit(-1);
        }

    }
    else if (npes == 1 ||npes == nspace)
    {
        if (root)
            for (int n = 1; n < ntime; n++)
                for (int i = 1; i < nspace - 1; i++)
                    results[n][i] = results[n - 1][i] + r * (results[n - 1][i + 1] - 2 * results[n - 1][i] + results[n - 1][i - 1]);
        // Prints solution matrix
        if (root)
        {
            printf("\n");
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
        // case 1 (good problem division)
        if (nspace % npes == 0)
        {
            int newspace = nspace / npes;
            double p_results[ntime + 1][newspace + 2];
            if (root)
                printf("newspace : %d\n", newspace);
            // INITIALISATION : FILL MATRIX WITH 0 AND BOUNDARY CONDITIONS
            for (int i = 0; i <= ntime; i++)
            {
                for (int j = 0; j <= newspace; j++)
                {
                    if (i == 0)
                        p_results[i][j] = Tint;
                    else
                        p_results[i][j] = 0.0;
                }

                if (root)
                    p_results[i][0] = Text;
                if (myrank == npes - 1)
                    p_results[i][newspace - 1] = Text;

                // for (int i = 0; i < ntime; i++)
                // {
                //     printf("rank : %d, line %d : ", myrank, i + 1);
                //     for (int j = 0; j < newspace; j++)
                //     {
                //         printf("%3.2f ", p_results[i][j]);
                //     }
                //     printf("\n");
                // }
            }

            // FILL UP LINE 1 to NTIME
            for (int n = 1; n <= ntime; n++)
            {
                if (myrank > 0)
                {
                    tag = 1;
                    MPI_Send(&p_results[n - 1][1], 1, MPI_DOUBLE, myrank - 1, tag, MPI_COMM_WORLD);
                }
                if (myrank < npes - 1)
                {
                    tag = 1;
                    MPI_Recv(&p_results[n - 1][newspace + 1], 1, MPI_DOUBLE, myrank + 1, tag, MPI_COMM_WORLD, &status);
                }
                if (myrank < npes - 1)
                {
                    tag = 2;
                    MPI_Send(&p_results[n - 1][newspace], 1, MPI_DOUBLE, myrank + 1, tag, MPI_COMM_WORLD);
                }
                if (myrank > 0)
                {
                    tag = 2;
                    MPI_Recv(&p_results[n - 1][0], 1, MPI_DOUBLE, myrank - 1, tag, MPI_COMM_WORLD, &status);
                }
                for (int i = 1; i <= newspace; i++)
                {
                    p_results[n][i] = p_results[n - 1][i] + r * (p_results[n - 1][i + 1] - 2 * p_results[n - 1][i] + p_results[n - 1][i - 1]);
                }

                // MPI_Gather(p_results[n], newspace, MPI_DOUBLE, results[n], newspace, MPI_DOUBLE, 0, MPI_COMM_WORLD);
                if (root)
                    p_results[n][0] = Text;
                if (myrank == npes - 1)
                    p_results[n][newspace - 1] = Text;
                
                // MPI_Gather(p_results[n], newspace, MPI_DOUBLE, results[n], newspace, MPI_DOUBLE, 0, MPI_COMM_WORLD);
            }
            // ALL RESULTS ARE GATHERED ON PROCESSOR 0
            for (int i = 0; i < ntime; i++)
            {
                MPI_Gather(p_results[i], newspace, MPI_DOUBLE, results[i], newspace, MPI_DOUBLE, 0, MPI_COMM_WORLD);
            }
            // OUTPUT: PRINTS PARTICULAR RESULTS
            // for (int i = 0; i < ntime; i++)
            // {
            //     printf("rank : %d, line %d : ", myrank, i + 1);
            //     for (int j = 0; j < newspace; j++)
            //     {
            //         printf("%3.2f ", p_results[i][j]);
            //     }
            //     printf("\n");
            // }
        }
        else
        {
            if (root)
            {
                printf("\n");
                printf("Wrong number of processors (%d) for the problem size (%d)\n", npes, nspace);
                printf("Please try a number of processors which is able to divide the problem size.\n");
                exit(-1);
            }
        }
    }
    if (root)
    {
        printf("\n\n");
        printf("Results matrix : \n\n");
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
    MPI_Finalize();
}