#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <mpi.h>

void thomasAlg(const double *a, const double *b, double *c, double *d, double *x, unsigned int n)
{
    int i;

    if (b[0] == 0)
    {
        printf("\n");
        printf("Warning : division by 0. Aborted.\n");
        return;
    }
    else
    {
        c[0] = c[0] / b[0];
        d[0] = d[0] / b[0];
    }

    double val;

    for (i = 1; i < n; i++)
    {
        if ((b[i] - c[i - 1] * a[i]) == 0)
        {
            printf("\n");
            printf("Warning : division by 0. Aborted.\n");
            exit(-1);
        }
        else
        {
            val = 1.0 / (b[i] - c[i - 1] * a[i]);
            c[i] = c[i] * val;
            d[i] = (d[i] - a[i] * d[i - 1]) * val;
        }
    }
    x[n - 1] = d[n - 1];
    for (i = n - 2; i > -1; i--)
        x[i] = d[i] - c[i] * x[i + 1];
}

void thomasAlg_parallel(int myrank, int npes, int N, double *b, double *a, double *c, double *x, double *q)
{
    int i, j, k, i_global;
    int newspace, p_offset;
    double S[2][2], T[2][2], s1tmp, s2tmp;
    MPI_Status status;
    double *l = (double *)malloc(N * sizeof(double));
    double *d = (double *)malloc(N * sizeof(double));
    double *y = (double *)malloc(N * sizeof(double));
    for (i = 0; i < N; i++)
        l[i] = d[i] = y[i] = 0.0;
    S[0][0] = S[1][1] = 1.0;
    S[1][0] = S[0][1] = 0.0;
    newspace = (int)floor(N / npes);
    p_offset = myrank * newspace;

    // Form local products of R_k matrices
    if (myrank == 0)
    {
        s1tmp = a[p_offset] * S[0][0];
        S[1][0] = S[0][0];
        S[1][1] = S[0][1];
        S[0][1] = a[p_offset] * S[0][1];
        S[0][0] = s1tmp;
        for (i = 1; i < newspace; i++)
        {
            s1tmp = a[i + p_offset] * S[0][0] - b[i + p_offset - 1] * c[i + p_offset - 1] * S[1][0];
            s2tmp = a[i + p_offset] * S[0][1] - b[i + p_offset - 1] * c[i + p_offset - 1] * S[1][1];
            S[1][0] = S[0][0];
            S[1][1] = S[0][1];
            S[0][0] = s1tmp;
        }
        S[0][1] = s2tmp;
    }
    else
    {
        for (i = 0; i < newspace; i++)
        {
            s1tmp = a[i + p_offset] * S[0][0] - b[i + p_offset - 1] * c[i + p_offset - 1] * S[1][0];
            s2tmp = a[i + p_offset] * S[0][1] - b[i + p_offset - 1] * c[i + p_offset - 1] * S[1][1];
            S[1][0] = S[0][0];
            S[1][1] = S[0][1];
            S[0][0] = s1tmp;
            S[0][1] = s2tmp;
        }
    }
    // Full-recursive doubling algorithm for distribution
    for (i = 0; i <= log2(npes); i++)
    {
        if (myrank + pow(2, i) < npes)
            MPI_Send(S, 4, MPI_DOUBLE, (int)(myrank + pow(2, i)), 0, MPI_COMM_WORLD);
        if (myrank - pow(2, i) >= 0)
        {
            MPI_Recv(T, 4, MPI_DOUBLE, (int)(myrank - pow(2, i)), 0, MPI_COMM_WORLD, &status);
            s1tmp = S[0][0] * T[0][0] + S[0][1] * T[1][0];
            S[0][1] = S[0][0] * T[0][1] + S[0][1] * T[1][1];
            S[0][0] = s1tmp;
            s1tmp = S[1][0] * T[0][0] + S[1][1] * T[1][0];
            S[1][1] = S[1][0] * T[0][1] + S[1][1] * T[1][1];
            S[1][0] = s1tmp;
        }
    }
    //Calculate last d_k first so that it can be distributed, //and then do the distribution.
    d[p_offset + newspace - 1] = (S[0][0] + S[0][1]) / (S[1][0] + S[1][1]);
    if (myrank == 0)
    {
        MPI_Send(&d[p_offset + newspace - 1], 1, MPI_DOUBLE,
                 1, 0, MPI_COMM_WORLD);
    }
    else
    {
        MPI_Recv(&d[p_offset - 1], 1, MPI_DOUBLE, myrank - 1, 0, MPI_COMM_WORLD, &status);
        if (myrank != npes - 1)
            MPI_Send(&d[p_offset + newspace - 1], 1, MPI_DOUBLE,
                     myrank + 1, 0, MPI_COMM_WORLD);
    }
    // Compute in parallel the local values of d_k and l_k
    if (myrank == 0)
    {
        l[0] = 0;
        d[0] = a[0];
        for (i = 1; i < newspace - 1; i++)
        {
            l[p_offset + i] = b[p_offset + i - 1] / d[p_offset + i - 1];
            d[p_offset + i] = a[p_offset + i] - l[p_offset + i] * c[p_offset + i - 1];
        }
        l[p_offset + newspace - 1] = b[p_offset + newspace - 2] /
                                     d[p_offset + newspace - 2];
    }
    else
    {
        for (i = 0; i < newspace - 1; i++)
        {
            l[p_offset + i] = b[p_offset + i - 1] /
                              d[p_offset + i - 1];
            d[p_offset + i] = a[p_offset + i] -
                              l[p_offset + i] * c[p_offset + i - 1];
        }
        l[p_offset + newspace - 1] = b[p_offset + newspace - 2] /
                                     d[p_offset + newspace - 2];
    }
    /***************************************************************/
    if (myrank > 0)
        d[p_offset - 1] = 0;
    // Distribute d_k and l_k to all processes
    double *tmp = (double *)malloc(N * sizeof(double));
    for (i = 0; i < N; i++)
        tmp[i] = d[i];
    MPI_Allreduce(tmp, d, N, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    for (i = 0; i < N; i++)
        tmp[i] = l[i];
    MPI_Allreduce(tmp, l, N, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    free(tmp);
    if (myrank == 0)
    {
        /* Forward Substitution [L][y] = [q] */ y[0] = q[0];
        for (i = 1; i < N; i++)
            y[i] = q[i] - l[i] * y[i - 1];
        /* Backward Substitution [U][x] = [y] */ x[N - 1] = y[N - 1] / d[N - 1];
        for (i = N - 2; i >= 0; i--)
            x[i] = (y[i] - c[i] * x[i + 1]) / d[i];
    }
    free(l);
    free(y);
    free(d);
    return;
}
int main(int argc, char *argv[])
{
    // MPI THINGS
    int npes, myrank;
    MPI_Status status;
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &npes);
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
    MPI_Request request;

    // Initialization
    double D = 0.1;
    double dx = 0.04;
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

    double results[ntime][nspace];
    for (int n = 0; n < ntime; n++)
        for (int i = 0; i < nspace; i++)
            results[n][i] = 0.0;
    // Initialisation with boudary conditions
    for (int i = 0; i < nspace; i++)
    {
        results[0][i] = Tint;
    }
    for (int n = 0; n < ntime; n++)
    {
        results[n][0] = Text;
        results[n][nspace - 1] = Text;
    }

    if (npes > nspace)
    {
        if (root)
        {
            printf("\n");
            printf("Too much processors for considered problem.\n");
            exit(-1);
        }
    }
    else if (npes == 1 || npes == nspace)
    {
        double a[nspace - 2];
        double b[nspace - 2];
        double c[nspace - 2];
        double c_copy[nspace - 2];
        double x[nspace - 2];
        double d[nspace];
        double d_reduced[nspace - 2];
        // Creating matrix C and vectors.
        for (int i = 0; i < nspace - 2; i++)
        {
            a[i] = -r;
            c[i] = -r;
            b[i] = 1 + 2 * r;
            x[i] = 0.0;
        }
        a[0] = 0.0;
        c[nspace - 3] = 0.0;

        for (int n = 1; n < ntime; n++)
        {
            // COPY THE PREVIOUS T° LINE IN d
            for (int i = 0; i < nspace; i++)
                d[i] = results[n - 1][i];
            // CHANGES THE BOUNDARIES VALUES
            d[1] += r * d[0];
            d[nspace - 2] += r * d[nspace - 1];
            for (int i = 0; i < nspace - 2; i++)
            {
                // COPY THE REDUCED VECTOR IN d_reduced
                d_reduced[i] = d[i + 1];
                // COPY GOOD VALUE OF C
                c_copy[i] = c[i];
            }
            // RESOLVES a[i]x[i-1]+b[i]x[i]+c[i]x[i+1] = d[i] // WARNING : C AND X ARE MODIFIED
            thomasAlg(a, b, c_copy, d_reduced, x, nspace - 2);

            // COPY THE RESULT x INTO RESULTS MATRIX
            for (int i = 0; i < nspace - 2; i++)
                results[n][i + 1] = x[i];
            results[n][0] = Text;
            results[n][nspace - 1] = Text;
        }
        if (root)
        {
            printf("\n\n");
            printf("Results matrix : \n\n");
            for (int i = 0; i < ntime; i++)
            {
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
        // if (nspace % npes == 0)
        // {

            double a[nspace - 2];
            double b[nspace - 2];
            double c[nspace - 2];
            double c_copy[nspace - 2];
            double x[nspace - 2];
            double d[nspace];
            double d_reduced[nspace - 2];
            int offset;
            // Creating matrix C and vectors.
            for (int i = 0; i < nspace - 2; i++)
            {
                a[i] = -r;
                c[i] = -r;
                b[i] = 1 + 2 * r;
                x[i] = 0.0;
            }

            for (int n = 1; n < ntime; n++)
            {
                // COPY THE PREVIOUS T° LINE IN d
                for (int i = 0; i < nspace; i++)
                    d[i] = results[n - 1][i];
                // CHANGES THE BOUNDARIES VALUES
                d[1] += r * d[0];
                d[nspace - 2] += r * d[nspace - 1];
                for (int i = 0; i < nspace - 2; i++)
                {
                    // COPY THE REDUCED VECTOR IN d_reduced
                    d_reduced[i] = d[i + 1];
                    // COPY GOOD VALUE OF C
                    c_copy[i] = c[i];
                }
                // RESOLVES a[i]x[i-1]+b[i]x[i]+c[i]x[i+1] = d[i] // WARNING : C AND X ARE MODIFIED
                thomasAlg_parallel(myrank, npes, nspace - 2, a, b, c_copy, x, d_reduced);

                // COPY THE RESULT x INTO RESULTS MATRIX
                for (int i = 0; i < nspace - 2; i++)
                    results[n][i + 1] = x[i];
                results[n][0] = Text;
                results[n][nspace - 1] = Text;
            }
            if (root)
            {
                printf("\n\n");
                printf("Results matrix ahem : \n\n");
                for (int i = 0; i < ntime; i++)
                {
                    for (int j = 0; j < nspace; j++)
                    {
                        printf("%3.2f ", results[i][j]);
                    }
                    printf("\n");
                }
            }
        }

        // if (root)
        // {
        //     printf("\n");
        //     printf("-----TODO------\n");
        // }
    // }
    MPI_Finalize();
}
