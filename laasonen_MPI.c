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
    double dx = 0.2;
    double dt = 0.1;
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
        c[nspace - 3] = 0;

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
        if (nspace % npes == 0)
        {
            int newspace = nspace / npes;
            double a[newspace];
            double b[newspace];
            double c[newspace];
            double d[newspace];
            double d_reduced[newspace];
            double x[newspace];
            double x_upper_h[newspace];
            double x_lower_h[newspace];
            double x_part[newspace];
            for (int i = 0; i < newspace; i++)
            {
                a[i] = -r;
                c[i] = -r;
                b[i] = 1 + 2 * r;
                d[i] = Tint;
            }
            if (root)
                d[0] = Text;
            if (myrank == npes - 1)
                d[newspace - 1] = Text;
            for (int n = 1; n < ntime; n++)
            {
                // COPY THE PREVIOUS T° LINE IN d
                for (int i = 0; i < newspace; i++)
                {
                    d[i] = results[n][myrank * newspace + i];
                    if (root)
                        d[0] += r * Text;
                    if (myrank == npes - 1)
                        d[newspace - 1] += r * Text;
                }

                // STEP 1 .........
                // FORWARD ELIMINATION
                x_upper_h[0] = c[0] / b[0];
                x_lower_h[0] = d[0] / b[0];
                for (int i = 1; i < newspace; i++)
                {
                    double denom = b[i] - c[i] * x_upper_h[i - 1];
                    if (!denom)
                    {
                        printf("Division by 0. Aborted.\n");
                        exit(-1);
                    }
                    else
                    {

                        x_upper_h[i] = c[i] / denom;
                        x_lower_h[i] = (d[i] - a[i] * x_lower_h[i - 1]) / denom;
                    }
                }
                // BACK SUBSTITUTION
                x_part[newspace - 1] = x_lower_h[newspace - 1];
                x_lower_h[newspace - 1] = -x_upper_h[newspace - 1];
                x_upper_h[newspace - 1] = a[newspace - 1] / b[newspace - 1];
                for (int i = newspace - 2; i > 0; i--)
                {
                    x_part[i] = x_lower_h[i] - x_upper_h[i] * x_part[i + 1];
                    x_lower_h[i] = -x_upper_h[i] - x_lower_h[i + 1];
                    double denom = b[i] - c[i] * x_upper_h[i - 1];
                    if (!denom)
                    {
                        printf("Division by 0. Aborted.\n");
                        exit(-1);
                    }
                    else
                    {
                        x_upper_h[i] = -a[i] / denom;
                    }
                }
                // FORWARD SUBSTITUTION
                x_upper_h[0] = -x_upper_h[0];
                for (int i = 1; i < newspace; i++)
                {
                    x_upper_h[i] = -x_upper_h[i] * x_upper_h[i - 1];
                }

                // STEP 2..........

                double OutData[8*npes];
                OutData[0] = -1.0;
                OutData[1] = x_upper_h[0];
                OutData[2] = x_lower_h[0];
                OutData[3] = -x_part[0];
                OutData[4] = x_upper_h[newspace - 1];
                OutData[5] = x_lower_h[newspace - 1];
                OutData[6] = -1.0;
                OutData[7] = -x_part[newspace - 1];

                // CONCATENATE OUTDATA ARRAYS
                double log2p = log(npes) / log(2);
                if (pow(2, log2p) == npes)
                    log2p++;
                for (int i = 0; i < log2p; i++)
                {
                    int nxfer = 8 * pow(2, i);
                    double toP = (myrank - (int)pow(2, i) + 2 * npes) % npes;
                    double frP = (myrank + (int)pow(2, i)) % npes;
                    MPI_Isend(OutData, nxfer, MPI_DOUBLE, toP, 0, MPI_COMM_WORLD, &request);
                    MPI_Recv(&OutData[nxfer + 1], nxfer, MPI_DOUBLE, frP, 0, MPI_COMM_WORLD, &status);
                }
                double reduca[2 * npes - 2];
                double reducb[2 * npes - 2];
                double reducc[2 * npes - 2];
                double reducr[2 * npes - 2];
                double coeff[2 * npes - 2];
                // PUT OUTDATA INTO REDUCED TRIDIAGONAL FORM
                int nsig = 8 * npes;
                int ifirst = 8 * (npes - myrank) + 5;
                for (int i = 0; i < 2 * npes - 2; i++)
                {
                    int ibase = (ifirst + 4 * (i - 1)) % nsig;
                    printf("ibase = %d\n", ibase);
                    reduca[i] = OutData[ibase - 1];
                    reducb[i] = OutData[ibase];
                    reducc[i] = OutData[ibase + 1];
                    reducr[i] = OutData[ibase + 2];
                }
                thomasAlg(reduca, reducb, reducc, reducr, coeff, 2 * npes - 2);
                int upper_coeff, lower_coeff;
                if (myrank != 0)
                    upper_coeff = coeff[2 * myrank - 2];
                else
                    upper_coeff = 0;
                if (myrank != npes)
                    lower_coeff = coeff[2 * myrank - 1];
                else
                    lower_coeff = 0;

                // printf("From rank %d \n", myrank);
                for (int i = 0; i < newspace; i++)
                {
                    x[i] = x_part[i] + upper_coeff * x_upper_h[i] + lower_coeff * x_lower_h[i];
                    // printf("%3.2f ", x[i]);
                }
                // printf("\n");

                MPI_Gather(x, newspace, MPI_DOUBLE, results[n], newspace, MPI_DOUBLE, 0, MPI_COMM_WORLD);
                results[n][0] = Text;
                results[n][nspace - 1] = Text;
            }
            // ALL RESULTS ARE GATHERED ON PROCESSOR 0
            // for (int i = 0; i < ntime; i++)
            // {
            // }
            if (root)
            {
                for (int n = 0; n < ntime; n++)
                {
                    printf("\n");
                    for (int i = 0; i < nspace; i++)
                    {
                        printf("%3.2f ", results[n][i]);
                    }
                }
            }
        }

        // if (root)
        // {
        //     printf("\n");
        //     printf("-----TODO------\n");
        // }
    }
    MPI_Finalize();
}
