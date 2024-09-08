#include <iostream>
#include <iomanip>
#include <chrono>
#include <cmath>

#include "helpers.hpp"

using namespace std;

int main(int argc, char **argv)
{
    int max_iterations = 1000;
    double epsilon = 1.0e-3;
    bool verify = true, print_config = false;

    // default values for M rows and N columns
    int N = 12;
    int M = 12;

    process_input(argc, argv, N, M, max_iterations, epsilon, verify, print_config);

    if ( print_config )
        std::cout << "Configuration: m: " << M << ", n: " << N << ", max-iterations: " << max_iterations << ", epsilon: " << epsilon << std::endl;

    auto time_1 = chrono::high_resolution_clock::now(); 

    int i, j;
    double diffnorm;
    int iteration_count = 0;

    Mat U(M, N);
    Mat W(M, N); 

    // Init & Boundary
    for (i = 0; i < M; ++i) {
        for (j = 0; j < N; ++j) {
            W[i][j] = U[i][j] = 0.0;
        }

        W[i][0] = U[i][0] = 0.05; // left side
        W[i][N-1] = U[i][N-1] = 0.1; // right side
    }

    for (j = 0; j < N; ++j) {
        W[0][j] = U[0][j] = 0.02; // top 
        W[M - 1][j] = U[M - 1][j] = 0.2; // bottom 
    }
    // End init

    iteration_count = 0;
    do
    {
        iteration_count++;
        diffnorm = 0.0;

        // Compute new values (but not on boundary) 
        for (i = 1; i < M - 1; ++i)
        {
            for (j = 1; j < N - 1; ++j)
            {
                W[i][j] = (U[i][j + 1] + U[i][j - 1] + U[i + 1][j] + U[i - 1][j]) * 0.25;
                diffnorm += (W[i][j] - U[i][j]) * (W[i][j] - U[i][j]);
            }
        }

        // Only transfer the interior points
        for (i = 1; i < M - 1; ++i)
            for (j = 1; j < N - 1; ++j)
                U[i][j] = W[i][j];

        diffnorm = sqrt(diffnorm); // all processes need to know when to stop
        
    } while (epsilon <= diffnorm && iteration_count < max_iterations);
    
    auto time_2 = chrono::high_resolution_clock::now();

    // Print time measurements 
    cout << "Elapsed time: "; 
    cout << std::fixed << std::setprecision(4) << chrono::duration<double>(time_2 - time_1).count(); 
    cout << " seconds, iterations: " << iteration_count << endl; 
 
    // verification     
    if ( verify ) {
        Mat U_sequential(M, N); // init another matrix for the verification

        int iteration_count_seq = 0;
        heat2d_sequential(U_sequential, max_iterations, epsilon, iteration_count_seq); 

        cout << "Verification: " << ( U.compare(U_sequential) && iteration_count == iteration_count_seq ? "OK" : "NOT OK") << std::endl;
    }

    return 0;
}