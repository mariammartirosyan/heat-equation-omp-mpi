# Heat Equation (2D) with OpenMP and MPI

The project aims to compare the sequential, OpenMP, and MPI versions of the 2D heat equation solver, measure their performance and determine which approach is best suited for different use cases.

### Performance Testing Parameters
The measurements were done on ALMA cluster with 6 nodes with the following configurations:
1) `--m 2688 --n 4096 --epsilon 0.01 --max-iterations 1000`
2) `--m 2688 --n 4096 --epsilon 0.01 --max-iterations 2000`
3) `--m 1152 --n 1152 --epsilon 0.01 --max-iterations 1000`

## OpenMP

### Parallelized Loops
The following loops were parallelized using  OpenMP `for` construct as the computations within iterations are independent:
- Initialization loops for U and W matrices
- Loop for computing new values
- Loop for transferring computed values from W to U

### Static vs Dynamic Scheduling
- **Static scheduling:** Distributes the workload evenly across the available threads and is more efficient for this problem, where the workload remains consistent across iterations.
- **Dynamic scheduling:** Assigns iterations to threads based on their availability, but it was less efficient for this task.

### Results

#### Configuration 1 

##### Static Scheduling
| Threads | Execution Time (Sec) | Speedup |
|---------|----------------|---------|
| 1       | 51.1           | 0.9     |
| 2       | 25.5           | 1.9     |
| 4       | 14.5           | 3.4     |
| 8       | 11.0           | 4.5     |
| 16      | 10.0           | 5.0     |
| 32      | 11.8           | 4.2     |

##### Dynamic Scheduling
| Threads | Execution Time (Sec) | Speedup |
|---------|----------------|---------|
| 1       | 51.8           | 0.9     |
| 2       | 31.6           | 1.5     |
| 4       | 17.9           | 2.7     |
| 8       | 13.3           | 3.7     |
| 16      | 12.5           | 4.0     |
| 32      | 11.6           | 4.3     |

#### Configuration 2 

##### Static Scheduling
| Threads | Execution Time (Sec) | Speedup |
|---------|----------------|---------|
| 1       | 112.6          | 0.9     |
| 2       | 55.9           | 2.0     |
| 4       | 34.8           | 3.2     |
| 8       | 24.8           | 4.5     |
| 16      | 21.3           | 5.2     |
| 32      | 22.1           | 5.0     |

##### Dynamic Scheduling
| Threads | Execution Time (Sec) | Speedup |
|---------|----------------|---------|
| 1       | 113.8           | 0.9     |
| 2       | 67.9           | 1.6     |
| 4       | 38.1           | 2.9     |
| 8       | 27.6           | 4.0     |
| 16      | 26.6           | 4.2     |
| 32      | 25.4           | 4.4     |

#### Configuration 3 

##### Static Scheduling
| Threads | Execution Time (Sec) | Speedup |
|---------|----------------|---------|
| 1       | 4.9          | 1.0    |
| 2       | 2.2           | 2.2     |
| 4       | 1.3           | 3.8     |
| 8       | 0.8           | 6.2     |
| 16      | 0.5           | 10.0     |
| 32      | 0.5           | 10.0     |

##### Dynamic Scheduling
| Threads | Execution Time (Sec) | Speedup |
|---------|----------------|---------|
| 1       | 4.9           | 1.0     |
| 2       | 4.1           | 1.2     |
| 4       | 2.1           | 2.3     |
| 8       | 1.2           | 4.1     |
| 16      | 0.8           | 6.2     |
| 32      | 0.8           | 6.2     |

---

## MPI

### Data Distribution
Each MPI process works on `M / numprocs + 2` rows of the matrix, where the extra 2 rows are for the halo/ghost rows. If `M % numprocs != 0`, leftover rows are assigned to the last process.

### Communication Between Processes
Every process except for rank 0 sends its second row to the process above and receives the second-to-last row from the process above. Also, every process except for the last one sends its second-to-last row to the process below and receives the second row from the process above.

The following MPI routines are used:
- **MPI_Isend** and **MPI_Irecv**: Used for non-blocking communication of halo rows.
- **MPI_Wait**: Waits for the processes to send/receive messages. In the meantime, the rows that don’t need the halo/ghost rows for their computations are being processed.
- **MPI_Allreduce**: Aggregates `diffnorm` across all processes to allow them to determine when to stop.

### Data Collection
After calculations are complete, data from all processes is gathered by rank 0 using `MPI_Gatherv`. This routine collects all computed rows, excluding the halo rows, and assembles them into the final matrix.

### Results

#### Configuration 1 

| Nodes | Processes | Speedup | Execution Time (Sec) | Gathering Time (Sec) |
|-------|-----------|---------|----------------|----------------|
| 1     | 2         | 1.5     | 32.4           | 0.01           |
| 1     | 4         | 2.3     | 21.3           | 0.01           |
| 1     | 8         | 3.1     | 15.7           | 0.01           |
| 1     | 16        | 6.3     | 7.9            | 0.02           |
| 1     | 32        | 6.4     | 7.8            | 0.04           |
| 2     | 32        | 9.4     | 5.3            | 0.3            |
| 4     | 64        | 12.8    | 3.9            | 0.5            |
| 6     | 96        | 14.2    | 3.5            | 0.6            |

#### Configuration 2 

| Nodes | Processes | Speedup | Execution Time (Sec) | Gathering Time (Sec) |
|-------|-----------|---------|----------------|----------------|
| 1     | 2         | 1.6     | 69.5           | 0.01           |
| 1     | 4         | 2.3     | 48.1           | 0.01           |
| 1     | 8         | 3.3     | 33.4           | 0.01           |
| 1     | 16        | 6.4     | 17.4            | 0.02           |
| 1     | 32        | 6.7     | 16.6           | 0.02           |
| 2     | 32        | 9.3     | 12.0            | 0.3            |
| 4     | 64        | 12.3    | 9.1            | 0.5            |
| 6     | 96        | 14.1    | 7.9            | 0.6            |

#### Configuration 3

| Nodes | Processes | Speedup | Execution Time (Sec) | Gathering Time (Sec) |
|-------|-----------|---------|----------------|----------------|
| 1     | 2         | 1.5     | 3.2           | 0.002           |
| 1     | 4         | 2.5     | 2.0           | 0.003           |
| 1     | 8         | 4.1     | 1.2           | 0.003           |
| 1     | 16        | 8.3     | 0.6            | 0.003          |
| 1     | 32        | 10.0     | 0.5            | 0.006          |
| 2     | 32        | 10.0     | 0.5            | 0.05            |
| 4     | 64        | 8.3    | 0.6            | 0.08            |
| 6     | 96        | 8.3    | 0.6            | 0.11           |
---
## OpenMP vs MPI

- **OpenMP:** A shared memory model, suitable for systems where all threads have access to the same memory. It requires minimal changes to the sequential code, mostly through pragmas. However, synchronization (e.g., for `diffnorm`) can introduce overhead.
  
- **MPI:** A distributed memory model, more appropriate for multi-node clusters. MPI is more complex due to the need for explicit communication between processes, which can introduce risks like deadlocks. 

### Matrix Size and Performance
Matrix size significantly affects performance due to synchronization/communication overhead:
- **OpenMP:** Synchronization overhead increases with matrix size.
- **MPI:** Larger matrices lead to higher inter-process communication overhead as data is transferred between processes.
---

## Conclusion

- **OpenMP**: Works well for smaller matrix sizes, but synchronization overhead increases with larger matrices.
- **MPI**: Outperforms OpenMP for large matrices in multi-node configurations due to better scalability and reduced synchronization overhead. 