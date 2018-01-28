#include <stdio.h>
#include <stdlib.h>
#include <math.h>

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

    printf("\n");
    printf("nspace : %d, ntime : %d\n", nspace, ntime);

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

    printf("\n\n");
    printf("Initialized matrix : \n\n");
    for (int i = 0; i < ntime; i++)
    {
        for (int j = 0; j < nspace; j++)
        {
            printf("%3.2f ", results[i][j]);
        }
        printf("\n");
    }
    // for (int i = 0; i < nspace; i++)
    // {
    //     // using an order one scheme to find n = 1
    //     results[1][i] = r * (results[0][i + 1] + results[0][i - 1]) + (1 - 2 * r) * results[0][i];
    // }
    // results[1][0] = Text;
    // results[1][nspace - 1] = Text;

    double a[nspace - 2];
    double b[nspace - 2];
    double c[nspace - 2];
    double e[nspace - 2];
    double f[nspace - 2];
    double g[nspace - 2];
    double c_copy[nspace - 2];
    double g_copy[nspace - 2];
    double x[nspace - 2];
    double x_new[nspace - 2];
    double d[nspace];
    double d_reduced[nspace - 2];

    // Creating matrices A and E and vectors.
    for (int i = 0; i < nspace - 2; i++)
    {
        a[i] = -r / 2;
        c[i] = -r / 2;
        b[i] = 1 + r;
        e[i] = r / 2;
        g[i] = r / 2;
        f[i] = 1 - r;
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
            g_copy[i] = g[i];
        }

        // CALCULATING NEW VECTOR X

        // x_new = E * d_reduced;
        x_new[0] = f[0] * d_reduced[0] + g[0] * d_reduced[1];
        for (int i = 1; i < nspace - 3; i++)
        {
            x_new[i] = e[0] * d_reduced[i - 1] + f[0] * d_reduced[i] + g[0] * d_reduced[i + 1];
        }
        x_new[nspace - 3] = e[0] * d_reduced[nspace - 4] + f[0] * d_reduced[nspace - 3];

        // RESOLVES a[i]x[i-1]+b[i]x[i]+c[i]x[i+1] = d[i] // WARNING : C AND X ARE MODIFIEd
        thomasAlg(a, b, c_copy, x_new, x, nspace - 2);

        // COPY THE RESULT x INTO RESULTS MATRIX
        for (int i = 0; i < nspace - 2; i++)
            results[n][i + 1] = x[i];
        results[n][0] = Text;
        results[n][nspace - 1] = Text;
    }

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
    return 0;
}
