#include <math.h>
#include <stdio.h>

#define PI 3.141592

double calculateSum(int max_m, double D, double L, double t, double x)
{
    double sum = 0.0;
    for (int m = 1; m < max_m; m++)
    {
        double term_1 = -D * (m * PI / L) * (m * PI / L) * t;
        double term_2 = (1 - pow(-1, m)) / (m * PI);
        double term_3 = m * PI * x / L;
        sum += exp(term_1) * term_2 * sin(term_3);
    }
    return sum;
}

int main(int argc, char *argv[])
{
    // Initialization
    double D = 0.1;
    double dx = 0.005;
    double dt = 0.001;
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

    double analytical[ntime][nspace];

    for (int n = 0; n < ntime; n++)
    {
        for (int i = 0; i < nspace; i++)
        {
            double x = i * dx;
            double t = n * dt;
            double sum = 0.0;

            sum = calculateSum(1000, D, L, t, x);

            analytical[n][i] = Text + 2 * (Tint - Text) * sum;
        }
    }

    printf("\n\n");
    printf("Results matrix : \n\n");
    for (int i = 0; i < ntime; i++)
    {
        if (i == ntime / 2)
        for (int j = 0; j < nspace; j++)
        {
                printf("%3.2f, ", analytical[i][j]);
        }
        // printf("\n");
    }
    return 0;
}
