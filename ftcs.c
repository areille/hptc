#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

void print_matrix(double **matrix, int width, int height)
{
    printf("Iterations over space : %d \n", width);
    printf("Iterations over time : %d \n", height);
    for (int i = 0; i < height; i++)
    {
        for (int j = 0; j < width; j++)
        {
            printf("%f", matrix[i][j]);
        }
        printf("\n");
    }
}

int main(int argc, char *argv[])
{
    // MPI things
    int npes, myrank;
    MPI_Status status;
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &npes);
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

    // printf("From process %d out of %d, Hello world!\n", myrank, npes);

    // Initialization
    double D = 0.1;
    double dx = 0.2;
    double dt = 0.2;
    double L = 1.0;
    double T = 1.0;
    double r = D * dt / pow(dx, 2);
    int Text = 300;
    int Tint = 100;
    int ntime = T / dt;
    int nspace = L / dx + 1;
    printf("nspace : %d, ntime : %d\n", nspace, ntime);

    double results[ntime][nspace];

    // Initialize grid
    for (int n = 0; n < ntime; n++)
    {
        for (int i = 0; i < nspace; i++)
        {
            results[n][i] = 0.0;
        }
    }
    for (int i = 0; i < nspace; i++)
    {
        results[0][i] = Tint;
    }
    for (int n = 0; n < ntime; n++)
    {
        results[n][0] = Text;
        results[n][nspace - 1] = Text;
    }
    results[1][0] = Text;
    results[1][nspace - 1] = Text;
    if (npes == 1)
    {
        // Fill the matrix according to equation
        for (int n = 1; n < ntime; n++)
        {
            for (int i = 1; i < nspace - 1; i++)
            {
                results[n][i] = results[n - 1][i] + r * (results[n - 1][i + 1] - 2 * results[n - 1][i] + results[n - 1][i - 1]);
            }
        }
    }
    if (npes == 2)
    {
        int newnspace = nspace / 2;
        printf("new nspace : %d\n", newnspace);
        double boundary_value_p1;
        double boundary_value_p2;

        // MPI_Type_vector(ntime, 1, newnspace, MPI_DOUBLE, &stype);
        // MPI_Type_commit(&stype);

        double p_results[ntime][newnspace];
        // INITIALISATION : FILL MATRIX WITH 0 AND BOUNDARY CONDITIONS
        for (int i = 0; i < ntime; i++)
        {
            for (int j = 0; j < newnspace; j++)
            {
                if (myrank == 0)
                    p_results[i][j] = results[i][j];
                else if (myrank == 1)
                    p_results[i][j] = results[i][j + newnspace];
            }
        }
        for (int n = 1; n < ntime; n++)
        {
            for (int i = 1; i < newnspace - 1; i++)
            {
                // cout << "n : " << n << ", i : " << i << endl;
                p_results[n][i] = p_results[n - 1][i] + r * (p_results[n - 1][i + 1] - 2 * p_results[n - 1][i] + p_results[n - 1][i - 1]);
                // cout << r << endl;
                // cout << "Results[n][i] : " << results1[n - 1][i] + r * (results1[n - 1][i + 1] - 2 * results1[n - 1][i] + results1[n - 1][i - 1]) << endl;
                // cout << "Results[n-1][i] : " << results1[n - 1][i] << endl;
                // cout << "Results[n-1][i-1] : " << results1[n - 1][i - 1] << endl;
                // cout << "Results[n-1][i+1] : " << results1[n - 1][i + 1] << endl;
            }
            if (myrank == 0)
            {
                boundary_value_p1 = p_results[n - 1][newnspace - 1];
                MPI_Send(&boundary_value_p1, 1, MPI_DOUBLE, 1, 0, MPI_COMM_WORLD);
                MPI_Recv(&boundary_value_p2, 1, MPI_DOUBLE, 1, 1, MPI_COMM_WORLD, &status);
                p_results[n][newnspace - 1] = p_results[n - 1][newnspace - 1] + r * (boundary_value_p2 - 2 * p_results[n - 1][newnspace - 1] + p_results[n - 1][newnspace - 2]);
            }
            else if (myrank == 1)
            {
                boundary_value_p2 = p_results[n - 1][0];
                MPI_Send(&boundary_value_p2, 1, MPI_DOUBLE, 0, 1, MPI_COMM_WORLD);
                MPI_Recv(&boundary_value_p1, 1, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD, &status);
                p_results[n][0] = p_results[n - 1][0] + r * (p_results[n - 1][1] - 2 * p_results[n - 1][0] + boundary_value_p1);
            }
        }

        // ALL RESULTS ARE GATHERED ON PROCESSOR 0
        for (int i = 0; i < ntime; i++)
        {
            MPI_Gather(p_results[i], newnspace, MPI_DOUBLE, results[i], newnspace, MPI_DOUBLE, 0, MPI_COMM_WORLD);
        }
    }
    if (myrank == 0)
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
    MPI_Finalize();
}
