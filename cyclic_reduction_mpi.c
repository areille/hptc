#include <mpi.h>
#include <stdio.h>
#include <math.h>

int main(int argc, char *argv[])
{

    /* PART 1 : INITIALIZATION */
    int i, j, k, size, index;
    int index1, index2;
    int npes, myrank;
    double alpha, gamma;
    const int numrows = 8; // = npes * 3 -1
    MPI_Status status;

    double D = 0.1;
    double dx = 0.2;
    double dt = 0.1;
    double L = 1.0;
    double T = 1.0;
    double r = D * dt / pow(dx, 2);

    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &npes);
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

    size = (int)pow(2, log2(npes + 1) + 1) - 1;
    // size = 15;
    double A[numrows][numrows];
    for (i = 0; i < numrows; i++)
    {
        for (j = 0; j < numrows; j++)
        {
            A[i][j] = 0.0;
        }
    }
    if (myrank == 0)
    {
        A[0][0] = 1 + 2 * r;
        A[0][1] = -r;
        A[1][0] = -r;
        A[1][1] = 1 + 2 * r;
        A[1][2] = -r;
        A[2][1] = -r;
        A[2][2] = 1 + 2 * r;
        A[2][3] = -r;
    }
    else if (myrank == (npes - 1))
    {
        index = 2 * myrank;
        A[0][index - 1] = -r;
        A[0][index] = 1 + 2 * r;
        A[0][index + 1] = -r;
        index = 2 * myrank + 1;
        A[1][index - 1] = -r;
        A[1][index] = 1 + 2 * r;
        A[1][index + 1] = -r;
        A[2][size - 2] = -r;
        A[2][size - 1] = 1 + 2 * r;
    }
    else
    {
        for (i = 0; i < 3; i++)
        {
            index = i + 2 * myrank;
            A[i][index - 1] = -r;
            A[i][index] = 1 + 2 * r;
            A[i][index + 1] = -r;
        }
    }
    // // RESULTS
    // if (myrank == 0)
    // {
    //     printf("\n");
    //     for (i = 0; i < size; i++)
    //     {
    //         for (j = 0; j < size; j++)
    //         {
    //             printf("%3.2f ", A[i][j]);
    //         }
    //         printf("\n");
    //     }
    // }
    for (i = 0; i < 3; i++)
    {
        A[i][size] = (double)(2 * myrank + i);
    }

    int numactivep = npes;
    int activep[npes];
    for (j = 0; j < numactivep; j++)
    {
        activep[j] = j;
    }
    for (j = 0; j < size + 1; j++)
    {
        A[3][j] = A[0][j];
        A[4][j] = A[2][j];
    }

    /* PART 2 : CYCLIC REDUCTION */

    for (i = 0; i < log2(size + 1) - 1; i++)
    {
        for (j = 0; j < numactivep; j++)
        {
            if (myrank == activep[j])
            {
                index1 = 2 * myrank + 1 - pow(2, i);
                index2 = 2 * myrank + 1 + pow(2, i);
                alpha = A[1][index1] / A[3][index1];
                gamma = A[1][index2] / A[4][index2];
                for (k = 0; k < size + 1; k++)
                    A[1][k] -= (alpha * A[3][k] + gamma * A[4][k]);
                if (numactivep > 1)
                {
                    if (j == 0)
                    {
                        MPI_Send(A[1], size + 1, MPI_DOUBLE, activep[1], 0,
                                 MPI_COMM_WORLD);
                    }
                    else if (j == numactivep - 1)
                    {
                        MPI_Send(A[1], size + 1, MPI_DOUBLE, activep[numactivep - 2],
                                 1, MPI_COMM_WORLD);
                    }
                    else if (j % 2 == 0)
                    {
                        MPI_Send(A[1], size + 1, MPI_DOUBLE, activep[j - 1],
                                 1, MPI_COMM_WORLD);
                        MPI_Send(A[1], size + 1, MPI_DOUBLE, activep[j + 1],
                                 0, MPI_COMM_WORLD);
                    }
                    else
                    {
                        MPI_Recv(A[3], size + 1, MPI_DOUBLE, activep[j - 1], 0,
                                 MPI_COMM_WORLD, &status);
                        MPI_Recv(A[4], size + 1, MPI_DOUBLE, activep[j + 1], 1,
                                 MPI_COMM_WORLD, &status);
                    }
                }
            }
        }
        numactivep = 0;
        for (j = activep[1]; j < npes; j = j + pow(2, i + 1))
        {
            activep[numactivep++] = j;
        }
    }

    /* PART 3 : BACK SUBSTITUTION */

    double x[npes];
    for (j = 0; j < npes; j++)
        x[j] = 0.0;
    if (myrank == activep[0])
    {
        x[myrank] = A[1][size] / A[1][(size - 1) / 2];
    }
    double tmp;
    for (i = log2(size + 1) - 3; i >= 0; i--)
    {
        tmp = x[myrank];
        MPI_Allgather(&tmp, 1, MPI_DOUBLE, x, 1, MPI_DOUBLE,
                      MPI_COMM_WORLD);
        numactivep = 0;
        for (j = activep[0] - pow(2, i); j < npes; j = j + pow(2, i + 1))
        {
            activep[numactivep++] = j;
        }
        for (j = 0; j < numactivep; j++)
        {
            if (myrank == activep[j])
            {
                x[myrank] = A[1][size];
                for (k = 0; k < npes; k++)
                {
                    if (k != myrank)
                        x[myrank] -= A[1][2 * k + 1] * x[k];
                }
                x[myrank] = x[myrank] / A[1][2 * myrank + 1];
            }
        }
    }
    tmp = x[myrank];
    // printf("rank : %d, x = %3.2f\n", myrank, tmp);
    MPI_Allgather(&tmp, 1, MPI_DOUBLE, x, 1, MPI_DOUBLE,
                  MPI_COMM_WORLD);

    // printf("x = %3.2f\n", x[0]);
    // printf("x = %3.2f\n", x[1]);
    // printf("x = %3.2f\n", x[2]);
    /* PART 4 : SOLVING FOR ODD ROWS */

    for (k = 0; k < npes; k++)
    {
        A[0][size] -= A[0][2 * k + 1] * x[k];
        A[2][size] -= A[2][2 * k + 1] * x[k];
    }
    A[0][size] = A[0][size] / A[0][2 * myrank];
    A[1][size] = x[myrank];
    A[2][size] = A[2][size] / A[2][2 * myrank + 2];

    // RESULTS
    if (myrank == 1)
    {
        printf("Rank 1 :\n");
        printf("\n");
        for (i = 0; i < size; i++)
        {
            for (j = 0; j < size; j++)
            {
                printf("%3.2f ", A[i][j]);
            }
            printf("\n");
        }
    }
    MPI_Finalize();
}
