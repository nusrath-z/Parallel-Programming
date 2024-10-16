#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <vector>
#include <iostream>
#include <string.h>

using namespace std;


 int getNextGroup (int i, int id , int d);
 int choosePivot(int A[], int size);
 void split(int A[], int size, int pivot, int &midIndex, int &lowLen, int &highLen);
 bool shouldPassLargerList (int iteration, int comm_rank, int d);
  int getCommLink (int iteration, int comm_rank, int dimensions);
  void sort(int A[], int size) ;

void printBuffer(int A[], int size) {
    std::cout << "[";
    for (int i = 0; i < size; i++) {
        std::cout << A[i];
        if (i < size - 1) {
            std::cout << ", ";
        }
    }
    std::cout << "]";
}
  
 //void performPivotSelection(int id, int num_process, std::vector<int>& new_chunk_array,MPI_Comm currComm);


int main(int argc, char *argv[]) {
    
int size = atoi(argv[1]),                       // size of array
        arr[size],                                  // declare array
        sorted_array[size], 
        id, 
        num_process,
        d,                       // declare sorted array
        i;

  double local_start_time = MPI_Wtime();
    
    
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &id);
    MPI_Comm_size(MPI_COMM_WORLD, &num_process);

    

    d = static_cast<int>(log2(num_process));
    int gp_size = num_process/2;

    if (id == 0){

   /**Initialization checks*/
    if (num_process <= 0 || (num_process & (num_process - 1)) != 0) {
        std::cerr << "Number of processes must be a power of 2.";
        MPI_Abort(MPI_COMM_WORLD, 1);}

    if(num_process > size){
          std::cerr << "Number of processes must less than or equal to number of elements.";
          MPI_Abort(MPI_COMM_WORLD, 1);}    

    }
    
    /**Generate random array
     * Prints the unsorted array
     */
    if(id==0){ 
          
          printf("Number of processes:%d\n",num_process);  
          printf("It is %d dimension\n",d);                           
          
          for(i=0; i<size; i++){                  
            arr[i] = rand()%100 + 1;
          }

          printf("\nThis is the unsorted array:\n");
          for(i=0; i<size; i++){
            printf("%d ", arr[i]);
          }
       
          printf("\n \n");
             
    }

    /**
     * Calculates sendcounts and displacements for each process
     * to Scatter the array accordinly
     * 
    */
    int sendcounts[num_process],
        displs[num_process],
        chunk_array[100],
        chunk_size,
        rem =  size%num_process,
        sum = 0;
    
    
    for (i = 0; i < num_process; i++) {
        sendcounts[i] = size/num_process;
        if (rem > 0) {
            sendcounts[i]++;
            rem--;
        }

        displs[i] = sum;
        sum += sendcounts[i];
    }

    if (id == 0) {
        for (i = 0; i < num_process; i++) {
            printf("sendcounts[%d] = %d\tdispls[%d] = %d\n", i, sendcounts[i], i, displs[i]);
        }
                                                                                                                                                                                                                                                                                                                                                                                                                                   
    }

    /**Scatter the array to each process from the root process*/
    MPI_Scatterv(arr, sendcounts, displs, MPI_INT, &chunk_array, 100, MPI_INT, 0, MPI_COMM_WORLD);

    /**new_chunk_array is the chunk recieved by each processor**/
    vector<int> new_chunk_array; 

chunk_size = sendcounts[id]; 
    if(id == 0){
    
        
        for(int i=0; i<chunk_size; i++){
            new_chunk_array.push_back(chunk_array[i]);
        }
        printf("process:%d chunck:",id);
        for(int i=0; i<chunk_size; i++){
          // printf("%d ",chunk_array[i]);
        }
        printf("\n \n");
    }
    else{

        //chunk_size = sendcounts[id];
        for(int i=0; i<chunk_size; i++){
           new_chunk_array.push_back(chunk_array[i]);
        }
        printf("process:%d chunck:",id);
        for(int i=0; i<chunk_size; i++){
           // printf("%d ", chunk_array[i]);
        }
        printf("\n \n");
    
    }

/**COMMUNICATION AND EXCHANGE*/
    MPI_Comm currComm = MPI_COMM_WORLD;
    for(i=1;i<=d;i++){

    int currentRank;
    int pivot;
    int currentSize;
    MPI_Comm_rank(currComm, &currentRank);
    MPI_Comm_size(currComm, &currentSize);

    
    
    if(currentRank==0){

       //printf("\ndimension: %d\n",i);

       pivot = choosePivot(chunk_array, chunk_size);
       

    }   

     MPI_Bcast(&pivot, 1, MPI_INT, 0, currComm);

     if (currentRank == 0){
          
      // printf("[MPI process %d] I am the broadcast root, and send value %d.\n", currentRank, pivot);
    } 
    else {
       
        //printf("[MPI process %d] I am a broadcast receiver, and obtained value %d.\n", currentRank,pivot);
    }

    int midIndex, lowLen, highLen;
    split(chunk_array, chunk_size, pivot, midIndex, lowLen, highLen);

   printf("\npivot: %d\n",pivot);
    // std::cout << "Elements after partition with pivot (" << pivot << ") : ";
    // int j;
    // for (j = 0; j < chunk_size; ++j) {
    //     std::cout << chunk_array[j] << " ";
    // }
    // std::cout << "(Length: " << lowLen << ")" << std::endl;
    // std::cout << "(Length: " << highLen << ")" << std::endl;
    
    // determine details for this communication
    bool PassLargerList = shouldPassLargerList(i,id,d);
    int commLink = getCommLink(i,id,d);
    int recvLen;

    if (PassLargerList) {
      MPI_Sendrecv(&highLen, 1, MPI_INT, commLink, 0, &recvLen, 1, MPI_INT, commLink, 0, MPI_COMM_WORLD, NULL);
    } else {
      MPI_Sendrecv(&lowLen, 1, MPI_INT, commLink, 0, &recvLen, 1, MPI_INT, commLink, 0, MPI_COMM_WORLD, NULL);
    }

    // initialize new array
    int keepLen = PassLargerList ? lowLen : highLen;
    int recvBuffer[recvLen + keepLen];

    if (PassLargerList) {
      MPI_Sendrecv(&chunk_array[midIndex], highLen, MPI_INT, commLink, 1, recvBuffer, recvLen, MPI_INT, commLink, 1, MPI_COMM_WORLD, NULL);
    } else {
      MPI_Sendrecv(chunk_array, lowLen, MPI_INT, commLink, 1, recvBuffer, recvLen, MPI_INT, commLink, 1, MPI_COMM_WORLD, NULL);
    }

    if (PassLargerList) {
    memcpy(&recvBuffer[recvLen], chunk_array, lowLen * sizeof(int));
} else {
    memcpy(&recvBuffer[recvLen], &chunk_array[midIndex], highLen * sizeof(int));
}

int z;
for (z = 0; z < recvLen + keepLen; ++z) {
    chunk_array[z] = recvBuffer[z];
}

chunk_size = recvLen + keepLen;

    // split communicator
    MPI_Comm nextComm;
    int nextGroup = getNextGroup(i, id,d);
    MPI_Comm_split(currComm, nextGroup, id, &nextComm);
    currComm = nextComm;

 }
   
   sort(chunk_array, chunk_size);

int* sortMe;
int *gatherSizes, *dis;
  if (id==0) { 
    gatherSizes = (int*) malloc(num_process * sizeof(int)); 
    dis = (int*) malloc(num_process * sizeof(int));
  }
  MPI_Gather(&chunk_size, 1, MPI_INT, gatherSizes, 1, MPI_INT, 0, MPI_COMM_WORLD);
  if (id==0) {
    dis[0] = 0;
    for (int i = 1; i < num_process; i++) {
      dis[i] = dis[i-1] + gatherSizes[i-1];
    }
    sortMe = (int*) malloc(size * sizeof(int));
  }
  MPI_Gatherv(chunk_array, chunk_size, MPI_INT, sortMe, gatherSizes, dis, MPI_INT, 0, MPI_COMM_WORLD);

if (id==0) {
    // Printing the sorted array received after gathering
    std::cout << "Sorted Array: ";
    printBuffer(sortMe, size); // Assuming you have a printBuffer function defined

    // Free the allocated memory
    free(sortMe);
    free(gatherSizes);
    free(dis);
}
   

double local_end_time = MPI_Wtime();
double global_start_time, global_end_time;

MPI_Allreduce(&local_start_time, &global_start_time, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
MPI_Allreduce(&local_end_time, &global_end_time, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);

if (id == 0) {
    double execution_time = global_end_time - global_start_time;
    printf("Execution time: %.6f seconds\n", execution_time);
}
     MPI_Finalize();
     return 0;
 }




  int getNextGroup (int i, int id , int d) {
    int bit = d - i;
    return (id & (1 << bit)) == 0 ? 0 : 1;
  }

int choosePivot(int A[], int size) {
    int firstIndex = rand() % size;
    int secondIndex = rand() % size;
    int thirdIndex = rand() % size;
    int firstSample = A[firstIndex];
    int secondSample = A[secondIndex];
    int thirdSample = A[thirdIndex];
    if (firstSample >= secondSample && firstSample < thirdSample) {
        return firstSample;
    }
    if (secondSample >= firstSample && secondSample < thirdSample) {
        return secondSample;
    }
    return thirdSample;
}

void split(int A[], int size, int pivot, int &midIndex, int &lowLen, int &highLen) {
    int i = 0;
    lowLen = 0;
    while (i < size && A[i] <= pivot) {
        i += 1;
        lowLen += 1;
    }
    midIndex = i;
    while (i < size) {
        if (A[i] <= pivot) {
            int tmp = A[i];
            A[i] = A[midIndex];
            A[midIndex] = tmp;
            midIndex += 1;
            lowLen += 1;
        }
        i += 1;
    }
    highLen = size - lowLen;
}

void sort(int A[], int size) {
    if (size > 1) {
        bool same = true;
        for (int i = 0; i < size - 1; i++) {
            if (A[i] != A[i + 1]) {
                same = false;
                break;
            }
        }
        if (same) {
            return;
        }
        int pivot = choosePivot(A, size);
        int midIndex, lowLen, highLen;
        split(A, size, pivot, midIndex, lowLen, highLen);
        sort(A, lowLen);
        sort(&A[midIndex], highLen);
    }
}

 bool shouldPassLargerList (int iteration, int comm_rank, int d) {
    int bit = d - iteration;
    return (comm_rank & (1 << bit)) == 0 ? true : false;
  }

   int getCommLink (int iteration, int comm_rank, int dimensions) {
    int bit = dimensions - iteration;
    int mask = (int) pow(2, bit);
    return comm_rank ^ mask;
  }


