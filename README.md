# Parallel-Programming
This repository contains two implementations of parallel computing using MPI (Message Passing Interface) in C++.
## Parallel Sorting Algorithms with MPI


---

## 1. Parallel Sample Sort (First Program) - f.cpp

### Description

This implementation performs a parallel sample sort on an array distributed across multiple MPI processes. The program:

* Generates a random integer array on process 0.
* Scatters chunks of the array to all processes.
* Locally sorts each chunk using quicksort.
* Performs regular sampling to choose pivots.
* Partitions each local chunk according to pivots.
* Exchanges partitions between processes.
* Merges received partitions and sorts them locally.
* Gathers the fully sorted partitions back to process 0.
* Outputs the final sorted array and execution time.

### Key Features

* Use of MPI collective communications such as `MPI_Scatterv`, `MPI_Gatherv`, `MPI_Bcast`.
* Regular sampling technique to select pivots for partitioning.
* Detailed partitioning and inter-process data exchange.
* Local sorting with quicksort.
* Execution time measurement using `MPI_Wtime()`.


## 2. Parallel Recursive Sorting with MPI Communicator Splitting (Second Program) - psrs.cpp

### Description

This implementation uses recursive sorting with MPI communicator splitting, based on a dimension `d = log2(num_processes)`. It works as follows:

* Generates a random array on process 0.
* Scatters parts of the array to all processes.
* Uses recursive pivot selection and partitioning.
* Communicates partitions between paired processes via `MPI_Sendrecv`.
* Splits MPI communicators recursively to progressively narrow sorting groups.
* Locally sorts final partitions.
* Gathers sorted partitions to process 0 for final output.
* Prints execution time.

### Key Features

* Recursive communicator splitting to handle sorting phases.
* Dynamic pivot selection with a median-of-three approach.
* MPI point-to-point communication (`MPI_Sendrecv`) for partition exchange.
* Local sorting via custom quicksort-like recursive function.
* Effective handling of powers-of-two number of processes.

