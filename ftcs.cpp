#include <iostream>
#include <mpi.h>
#include <math.h>
#include <vector>

using namespace std;

void print_matrix(vector< vector<double> > matrix, int width, int height){
    cout << "iterations over space : " << width << endl;
    cout << "iterations over time : " << height << endl;
    for (int i = 0; i < height; i++) {
        for (int j = 0; j<width; j++){
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

    vector< vector<double> > results(ntime, vector<double>(nspace));
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


    // Fill the matrix according to equation
    for (int n = 1; n<ntime; n++) {
        for (int i = 1; i<nspace-1; i++) {
            results[n][i] = results[n-1][i] + r * (results[n-1][i+1] - 2* results[n-1][i] + results[n-1][i-1]);
        }
    }

    cout << "Hello World! I'm processor " << myrank+1 << " of " << npes << endl;

    // Stop MPI
    MPI_Finalize();

    print_matrix(results, nspace, ntime);

    return 0;
}
