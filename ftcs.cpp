#include <iostream>
#include <mpi.h>
#include <math.h>
#include <vector>

using namespace std;

void print_matrix(vector<vector<double> > matrix, int width, int height)
{
    cout << "iterations over space : " << width << endl;
    cout << "iterations over time : " << height << endl;
    for (int i = 0; i < height; i++)
    {
        for (int j = 0; j < width; j++)
        {
            cout << matrix[i][j] << " ";
        }
        cout << endl;
    }
}
int main(int argc, char **argv)
{
    // MPI things
    int npes, myrank;
    MPI_Status status;
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &npes);
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
    MPI_Comm comm;
    MPI_Datatype stype;

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
    cout << nspace << " " << ntime;

    vector<vector<double> > results(ntime, vector<double>(nspace));
    // Initialize grid
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
        cout << "New : " << newnspace << endl;
        double boundary_value_p1;
        double boundary_value_p2;

        MPI_Type_vector(ntime, 1, newnspace, MPI_DOUBLE, &stype);
        MPI_Type_commit(&stype);

        vector<vector<double> > p_results(ntime, vector<double>(newnspace));
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
                cout << "Out of process " << myrank << ", sent." << endl;
                MPI_Recv(&boundary_value_p2, 1, MPI_DOUBLE, 1, 1, MPI_COMM_WORLD, &status);
                cout << "Out of process " << myrank << ", received." << endl;
                p_results[n][newnspace - 1] = p_results[n - 1][newnspace - 1] + r * (boundary_value_p2 - 2 * p_results[n - 1][newnspace - 1] + p_results[n - 1][newnspace - 2]);
            } else if (myrank == 1) {
                boundary_value_p2 = p_results[n - 1][0];
                MPI_Send(&boundary_value_p2, 1, MPI_DOUBLE, 0, 1, MPI_COMM_WORLD);
                cout << "Out of process " << myrank << ", sent." << endl;
                MPI_Recv(&boundary_value_p1, 1, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD, &status);
                cout << "Out of process " << myrank << ", received." << endl;
                // results2[n][0] = boundary_value_p1;
                p_results[n][0] = p_results[n - 1][0] + r * (p_results[n - 1][1] - 2 * p_results[n - 1][0] + boundary_value_p1);
            }
        }
        print_matrix(p_results, newnspace, ntime);
        // PUTTING BACK RESULTS MATRIX INSIDE MAIN RESULTS
        for (int i = 0; i < ntime; i++)
        {
            for (int j = 0; j < newnspace; j++)
            {
                if (myrank == 0)
                    results[i][j] = p_results[i][j];
                else if (myrank == 1)
                    results[i][j + newnspace] = p_results[i][j];
            }
        }
    }

    // Stop MPI
    MPI_Finalize();
    print_matrix(results, nspace, ntime);

    return 0;
}
