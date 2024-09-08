#include <iostream>
#include <iomanip>
#include <chrono>
#include <cmath>
#include <mpi.h>
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

    Mat U_total(M, N);
    int numprocs, rank;

    MPI_Init(&argc, &argv);
    auto start_time = MPI_Wtime(); 
    MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    int i, j;
    double diffnorm;
    int iteration_count = 0;

    int local_M = (M / numprocs) +2 + (rank==numprocs-1? M % numprocs:0);
    
    Mat U(local_M, N); 
    Mat W(local_M, N); 

    // Init & Boundary
    for (i = 0; i < local_M; ++i) {
        for (j = 0; j < N; ++j) {
            W[i][j] = U[i][j] = 0.0;
        }

        W[i][0] = U[i][0] = 0.05; // left side
        W[i][N-1] = U[i][N-1] = 0.1; // right side
    }

    if(rank==0){
        for (j = 0; j < N; ++j) {
            W[1][j] = U[1][j] = 0.02; // top 
        }
    }
    if(rank==numprocs-1){
        for (j = 0; j < N; ++j) {
            W[local_M - 2][j] = U[local_M - 2][j] = 0.2; // bottom 
        }
    }

    // End init
    int req_count = 0;
    iteration_count = 0;

    MPI_Status status;
    MPI_Request reqsend[2], reqrecv[2];

    const int start = (rank == 0) ? 2 : 1;
    const int end = (rank == numprocs-1) ? local_M - 3 : local_M-2;
    do
    {
        iteration_count++;
        diffnorm = 0.0;

        if(rank > 0){
            MPI_Isend(&U[1][1], N-2, MPI_DOUBLE, rank-1, 0, MPI_COMM_WORLD, &reqsend[0]); 
            MPI_Irecv(&U[0][1], N-2, MPI_DOUBLE, rank-1, 1, MPI_COMM_WORLD, &reqrecv[0]);

        }
        if(rank < numprocs-1){
            MPI_Isend(&U[local_M-2][1], N-2, MPI_DOUBLE, rank+1, 1, MPI_COMM_WORLD, &reqsend[1]);
            MPI_Irecv(&U[local_M-1][1], N-2, MPI_DOUBLE, rank+1, 0, MPI_COMM_WORLD, &reqrecv[1]);
        }
        for (i = start+1; i < end; ++i)
        {
            for (j = 1; j < N - 1; ++j)
            {
                W[i][j] = (U[i][j + 1] + U[i][j - 1] + U[i + 1][j] + U[i - 1][j]) * 0.25;
                diffnorm += (W[i][j] - U[i][j]) * (W[i][j] - U[i][j]);
            }
        }
   
       if (rank > 0){
            MPI_Wait(&reqrecv[0], &status);
            MPI_Wait(&reqsend[0], &status);
        }
        if (rank < numprocs - 1){
            MPI_Wait(&reqrecv[1], &status);
            MPI_Wait(&reqsend[1], &status);
        }
        for (j = 1; j < N - 1; ++j)
        {
            W[start][j] = (U[start][j + 1] + U[start][j - 1] + U[start + 1][j] + U[start - 1][j]) * 0.25;
            diffnorm += (W[start][j] - U[start][j]) * (W[start][j] - U[start][j]);
        }
        for (j = 1; j < N - 1; ++j)
        {
            W[end][j] = (U[end][j + 1] + U[end][j - 1] + U[end + 1][j] + U[end - 1][j]) * 0.25;
            diffnorm += (W[end][j] - U[end][j]) * (W[end][j] - U[end][j]);
        }
        // Only transfer the interior points
        for (i = start; i <= end; ++i)
            for (j = 1; j < N - 1; ++j)
                U[i][j] = W[i][j];
        

        MPI_Allreduce(MPI_IN_PLACE, &diffnorm, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
        diffnorm = sqrt(diffnorm);

    } while (epsilon <= diffnorm && iteration_count < max_iterations);
    
    
    auto end_time = MPI_Wtime();

    auto gather_start_time = MPI_Wtime();

    int sendcount = (local_M-2) * N; 
    int receive_counts[numprocs]; 
    int receive_displs[numprocs];

    for (int i = 0; i < numprocs; ++i) {
        receive_counts[i] = (local_M-2) * N;
        receive_displs[i] = (local_M-2) * N * i;
    }
    receive_counts[numprocs-1]+= (M % numprocs) * N;

    MPI_Gatherv(&U(1,0), sendcount, MPI_DOUBLE, &U_total(0,0), receive_counts, receive_displs, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    auto gather_end_time = MPI_Wtime();
    
    // verification     
    if (rank==0 && verify ) {
        // Print time measurements 
        auto gather_time = gather_end_time - gather_start_time;
        auto computation_time = end_time - start_time;
        cout << "Total Elapsed time: "; 
        cout << std::fixed << std::setprecision(4) <<computation_time+gather_time; 
        cout << " seconds, iterations: " << iteration_count << endl; 
        cout << "Gathering time: "; 
        cout << std::fixed << std::setprecision(4) << gather_time << endl; 
        cout << "Calculations time: "; 
        cout << std::fixed << std::setprecision(4) << computation_time << endl; 

        Mat U_sequential(M, N); // init another matrix for the verification
        int iteration_count_seq = 0;
        heat2d_sequential(U_sequential, max_iterations, epsilon, iteration_count_seq); 
        cout << "Verification: " << ( U_total.compare(U_sequential) && iteration_count == iteration_count_seq ? "OK" : "NOT OK") << std::endl;
    }
    MPI_Finalize();
    
    return 0;
}