#include <iostream>
#include <cstdlib>
#include <cstdio>
#include <cmath>
#include <ctime>
#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>


void quickSort(int arr[], int left, int right) {
    int i = left, j = right;
    int pivot = arr[(left + right) / 2];

    // Partition
    while (i <= j) {
        while (arr[i] < pivot) i++;
        while (arr[j] > pivot) j--;
        if (i <= j) {
            int temp = arr[i];
            arr[i] = arr[j];
            arr[j] = temp;
            i++;
            j--;
        }
    }

    // Recursion
    if (left < j) quickSort(arr, left, j);
    if (i < right) quickSort(arr, i, right);
}

int main(int argc, char *argv[]) {
    int size = atoi(argv[1]); // size of array
    int *arr = new int[size]; // declare array
    int rank, p,i;
    int start_time, end_time;

    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &p);


      if (rank == 0) {
        start_time = MPI_Wtime();
    }

 
    // if (rank == 0) {
    //     /**Initialization checks*/
    //     if (p <= 0 || (p & (p- 1)) != 0) {
    //         std::cerr << "Number of processes must be a power of 2.";
    //         MPI_Abort(MPI_COMM_WORLD, 1);
    //     }

    //     if (p> size) {
    //         std::cerr << "Number of processes must be less than or equal to the number of elements.";
    //         MPI_Abort(MPI_COMM_WORLD, 1);
    //     }
    // }

    /**Generate random array
     * Prints the unsorted array
     */
    if (rank == 0) {
        printf("Number of processes: %d\n", p);

        for (i = 0; i < size; i++) {
            arr[i] = rand() % 100 + 1;
        }

        printf("\nThis is the unsorted array:\n");
        for (i = 0; i < size; i++) {
            printf("%d ", arr[i]);
        }

        printf("\n \n");
    }

    /**
     * Calculates sendcounts and displacements for each process
     * to Scatter the array accordingly
     */
    int *sendcounts = new int[p];
    int *displs = new int[p];
    int *chunk_array = new int[100]; // Initialize with appropriate size
    int chunk_size, rem = size % p, sum = 0;

int local_size = size / p;

for (i = 0; i < p; i++) {
    sendcounts[i] = (i < rem) ? local_size + 1 : local_size;
    displs[i] = sum;
    sum += sendcounts[i];
}
    if (rank == 0) {
        for (i = 0; i < p
; i++) {
            printf("sendcounts[%d] = %d\tdispls[%d] = %d\n", i, sendcounts[i], i, displs[i]);
        }
    }

    /** Scatter the array to each process from the root process */
    MPI_Scatterv(arr, sendcounts, displs, MPI_INT, chunk_array, 100, MPI_INT, 0, MPI_COMM_WORLD);
chunk_size = sendcounts[rank]; 
printf("Process %d received elements: ", rank);
    for (i = 0; i < chunk_size; i++) {
        printf("%d ", chunk_array[i]);
    }
    printf("\n");

    quickSort(chunk_array, 0, chunk_size - 1);

    printf("Process %d received and sorted elements: ", rank);
    for (i = 0; i < chunk_size; i++) {
        printf("%d ", chunk_array[i]);
    }
    printf("\n");

// ...
int regular_sample_count = size / (p * p); // Number of regular samples per process
int *regular_samples = new int[regular_sample_count];

// Calculate indices for regular samples
for (i = 0; i < regular_sample_count; i++) {
    regular_samples[i] = chunk_array[i * (local_size / regular_sample_count)];
}

// Print regular samples for each process
printf("Process %d regular samples: ", rank);
for (i = 0; i < regular_sample_count; i++) {
    printf("%d ", regular_samples[i]);
}
printf("\n");


 // Gather regular samples to process P0
    int *recvcounts = new int[p];
    int *displs_gather = new int[p];

    // Gather counts and displacements for MPI_Gatherv
    MPI_Gather(&regular_sample_count, 1, MPI_INT, recvcounts, 1, MPI_INT, 0, MPI_COMM_WORLD);

    int total_samples = 0;
    int *recv_buf = NULL;
    int *displs_buf = NULL;

    if (rank == 0) {
        for (i = 0; i < p; i++) {
            displs_gather[i] = total_samples;
            total_samples += recvcounts[i];
        }

        recv_buf = new int[total_samples];
        displs_buf = new int[p];
    }

    MPI_Gatherv(regular_samples, regular_sample_count, MPI_INT, recv_buf, recvcounts, displs_gather, MPI_INT, 0, MPI_COMM_WORLD);

 int pivot_count = p - 1; // Number of regular samples per process
 int *pivot_samples = new int[pivot_count];
   if (rank == 0) {
        printf("Regular samples collected by process P0 (before sorting): ");
        for (i = 0; i < total_samples; i++) {
            printf("%d ", recv_buf[i]);
        }
        printf("\n");

        // Sort the received regular samples in process P0
        quickSort(recv_buf, 0, total_samples - 1);

        printf("Regular samples collected by process P0 (after sorting): ");
        for (i = 0; i < total_samples; i++) {
            printf("%d ", recv_buf[i]);
        }
        printf("\n");

       

/**Regular Sampling*/
int j;
//int step = (i + 1) * local_size + (local_size / 2) - 1;
for (j = 0; j < pivot_count; j++) {
    pivot_samples[j] = recv_buf[(j + 1) * p + (p / 2) - 1];
}
    
     printf("Chosen pivot values by process P0: ");
        for (i=0; i < pivot_count; i++) {
            printf("%d ", pivot_samples[i]);
        }
        printf("\n");

   }
   MPI_Bcast(pivot_samples, pivot_count, MPI_INT, 0, MPI_COMM_WORLD);

     // Print received pivot values in all processes
    printf("Process %d received pivot values: ", rank);
    for (i = 0; i < pivot_count; i++) {
        printf("%d ", pivot_samples[i]);
    }
    printf("\n");




/**Partition the local sorted sub-list into P partitions using received pivot values**/

int *partition_sizes = new int[p]; // Array to store sizes of each partition
int *partition_indices = new int[p + 1]; // Array to store starting indices of each partition

// Initialize partition_sizes array to 0
for (i = 0; i < p; i++) {
    partition_sizes[i] = 0;
}

// Determine partition sizes based on pivot values
int current_partition = 0;
for (i = 0; i < local_size; i++) {
    if (current_partition < pivot_count && chunk_array[i] > pivot_samples[current_partition]) {
        current_partition++;
    }
    partition_sizes[current_partition]++;
}

// Calculate starting indices of partitions
partition_indices[0] = 0;
for (i = 1; i <= p; i++) {
    partition_indices[i] = partition_indices[i - 1] + partition_sizes[i - 1];
}

// Allocate memory for partitions
int **partitions = new int*[p];
for (i = 0; i < p; i++) {
    partitions[i] = new int[partition_sizes[i]];
}

// Assign elements to partitions based on pivot values
int *current_indices = new int[p];
for (i = 0; i < p; i++) {
    current_indices[i] = partition_indices[i];
}

current_partition = 0;
for (i = 0; i < local_size; i++) {
    if (current_partition < pivot_count && chunk_array[i] > pivot_samples[current_partition]) {
        current_partition++;
    }
    partitions[current_partition][current_indices[current_partition] - partition_indices[current_partition]] = chunk_array[i];
    current_indices[current_partition]++;
}

// Print the partitions (for demonstration purposes)
for (i = 0; i < p; i++) {
    printf("Process %d, Partition %d: ", rank, i);
    for (int j = 0; j < partition_sizes[i]; j++) {
        printf("%d ", partitions[i][j]);
    }
    printf("\n");
}


//COMMMUNICATION

int *partition_lengths = new int[p]; // Array to store lengths of each partition to be sent
int *recv_partition = NULL; // Buffer to receive partition

// Determine partition lengths to be sent to other processes
for (i = 0; i < p; i++) {
    partition_lengths[i] = partition_sizes[i];
}

// Pi keeps its ith partition
int *local_partition = partitions[rank];
int local_partition_size = partition_sizes[rank];

// Send jth partition to process Pj (excluding own partition)
for (int j = 0; j < p; j++) {
    if (j != rank) {
        // Send the length of the partition to be received
        MPI_Send(&partition_lengths[j], 1, MPI_INT, j, 0, MPI_COMM_WORLD);

        // Send the partition elements to process Pj
        MPI_Send(partitions[j], partition_lengths[j], MPI_INT, j, 0, MPI_COMM_WORLD);
    }
}


// Receive jth partition from process Pj
for (int j = 0; j < p; j++) {
    if (j != rank) {
        int recv_length;

        // Receive the length of the partition to be received
        MPI_Recv(&recv_length, 1, MPI_INT, j, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

        // Allocate memory to receive the partition elements
        recv_partition = new int[recv_length];

        // Receive the partition elements from process Pj
        MPI_Recv(recv_partition, recv_length, MPI_INT, j, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

        // Replace the jth partition with the received partition elements
        delete[] partitions[j];
        partitions[j] = new int[recv_length];
        for (i = 0; i < recv_length; i++) {
            partitions[j][i] = recv_partition[i];
        }

        // Clean up memory for received partition
        delete[] recv_partition;
    }
}

// Calculate the total size of merged partitions
int total_merged_size = 0;
for (i = 0; i < p; i++) {
    total_merged_size += partition_sizes[i];
}

// Merge the partitions into a single array for each process
int *merged_partition = new int[total_merged_size];

// Merge all partitions into 'merged_partition'
int current_index = 0;
for (i = 0; i < p; i++) {
    for (int j = 0; j < partition_sizes[i]; j++) {
        merged_partition[current_index++] = partitions[i][j];
    }
}

// Sort the merged_partition locally within each process
quickSort(merged_partition, 0, total_merged_size - 1);

// Print the locally sorted merged partition for each process (for demonstration purposes)
printf("Process %d, Locally Sorted Merged Partition: ", rank);
for (i = 0; i < total_merged_size; i++) {
    printf("%d ", merged_partition[i]);
}
printf("\n");



// Gather locally sorted merged partitions from all processes onto Process 0
int *recv_merged_partitions = NULL;
int *recv_counts = NULL;
int *displs_gatherv = NULL;
if (rank == 0) {
    recv_merged_partitions = new int[total_merged_size * p];
    recv_counts = new int[p];
    displs_gatherv = new int[p];
}

// Gather counts of merged partitions to be received by Process 0
MPI_Gather(&total_merged_size, 1, MPI_INT, recv_counts, 1, MPI_INT, 0, MPI_COMM_WORLD);

// Calculate displacement for MPI_Gatherv
if (rank == 0) {
    displs_gatherv[0] = 0;
    for (i = 1; i < p; i++) {
        displs_gatherv[i] = displs_gatherv[i - 1] + recv_counts[i - 1];
    }
}

// Gather locally sorted merged partitions from all processes onto Process 0
MPI_Gatherv(merged_partition, total_merged_size, MPI_INT, recv_merged_partitions,
            recv_counts, displs_gatherv, MPI_INT, 0, MPI_COMM_WORLD);

// Sort the received merged partitions on Process 0
if (rank == 0) {
    // Calculate the total size of the received merged partitions
    int total_received_size = 0;
    for (i = 0; i < p; i++) {
        total_received_size += recv_counts[i];
    }

    // Sort the received merged partitions on Process 0
    quickSort(recv_merged_partitions, 0, total_received_size - 1);

    // Print the final sorted list on Process 0
    printf("Process %d, Final Sorted List: ", rank);
    for (i = 0; i < total_received_size; i++) {
        printf("%d ", recv_merged_partitions[i]);
    }
    printf("\n");

}

double local_start_time = MPI_Wtime();
double local_end_time = MPI_Wtime();
double global_start_time, global_end_time;

MPI_Allreduce(&local_start_time, &global_start_time, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
MPI_Allreduce(&local_end_time, &global_end_time, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);

if (rank == 0) {
    double execution_time = global_end_time - global_start_time;
    printf("Execution time: %.6f seconds\n", execution_time);
}


// Clean up memory
if (rank == 0) {
    delete[] recv_merged_partitions;
    delete[] recv_counts;
    delete[] displs_gatherv;
}

delete[] merged_partition;

delete[] partition_lengths;

for (i = 0; i < p; i++) {
    delete[] partitions[i];
}
delete[] partitions;

delete[] partition_sizes;
delete[] partition_indices;
delete[] current_indices;
    delete[] pivot_samples;
    delete[] recvcounts;
    delete[]displs_gather;
    delete[] regular_samples;
    delete[] arr;
    delete[] sendcounts;
    delete[] displs;
    delete[] chunk_array;

    MPI_Finalize();
    return 0;
}