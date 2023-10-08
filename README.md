# 2DHeatEquationSolver
Implemented parallel versions of the 2D heat equation solver with OpenMP and MPI.

OpenMP version: parallelizable opportunities were identified and the corresponding pragmas were applied.

MPI version: the work is distributed among processes. 
As each process works on a sub-part of the initial large matrix, and needs rows from other processes for its own calculations communication between those processes were established. 
