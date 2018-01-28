#include <mpi.h>
#include <stdio.h>
#include <math.h>

const int size = 15;

int main(int argc, char *argv[])
{
    int i, j, k;
    int index1, index2, offset;
    double alpha, gamma;
    double D = 0.1;
    double dx = 0.2;
    double dt = 0.1;
    double L = 1.0;
    double T = 1.0;
    double r = D * dt / pow(dx, 2);

    /* PART 1 : MEMORY ALLOCATION AND MATRIX GENERATION */
    /* TO IMPROVE : USE AS MUCH MEMORY AS NECESSARY, AND NOT TOO MUCH */
    double x[size];
    for (i = 0; i < size; i++)
        x[i] = 0.0;

    double F[size];
    double A[size][size];

    for (int i = 0; i < size; i++)
    {
        for (j = 0; j < size; j++)
        {
            A[i][j] = 0.0;
        }
        F[i] = 100.0;
    }
    F[0] = 300.0;
    F[size - 1] = 300.0;

    A[0][0] = 1 + 2 * r;
    A[0][1] = -r;
    A[size - 1][size - 2] = -r;
    A[size - 1][size - 1] = 1 + 2 * r;

    for (i = 0; i < size - 1; i++)
    {
        A[i][i] = 1 + 2 * r;
        A[i][i - 1] = -r;
        A[i][i + 1] = -r;
    }

    /* PART 2 : CYCLIC REDUCTION */

    for (i = 0; i < log2(size + 1) - 1; i++)
    {
        for (j = pow(2, i + 1) - 1; j < size; j = j + pow(2, i + 1))
        {
            offset = pow(2, i);
            index1 = j - offset;
            index2 = j + offset;

            alpha = A[j][index1] / A[index1][index1];
            gamma = A[j][index2] / A[index2][index2];

            for (k = 0; k < size; k++)
            {
                A[j][k] -= (alpha * A[index1][k] + gamma * A[index2][k]);
            }
            F[j] -= (alpha * F[index1] + gamma * F[index2]);
        }
    }

    /* PART 3 : CYCLIC REDUCTION */

    int index = (size - 1) / 2;
    x[index] = F[index] / A[index][index];
    for (i = log2(size + 1) - 2; i >= 0; i--)
    {
        for (j = pow(2, i + 1) - 1; j < size; j = j + pow(2, i + 1))
        {
            offset = pow(2, i);
            index1 = j - offset;
            index2 = j + offset;
            x[index1] = F[index1];
            x[index2] = F[index2];
            for (k = 0; k < size; k++)
            {
                if (k != index1)
                    x[index1] -= A[index1][k] * x[k];
                if (k != index2)
                    x[index2] -= A[index2][k] * x[k];
            }
            x[index1] = x[index1] / A[index1][index1];
            x[index2] = x[index2] / A[index2][index2];
        }
    }
    for (i = 0; i < size; i++)
    {
        printf("F[%d] = %3.2f\n", i, F[i]);
    }
    for (i = 0; i < size; i++)
    {
        printf("x[%d] = %3.2f\n", i, x[i]);
    }

    return 0;
}