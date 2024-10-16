#include <iostream>
#include <cstdlib>
#include <ctime>
#include <sys/time.h> 

using namespace std;

void swap(int &a, int &b) {
    int temp = a;
    a = b;
    b = temp;
}

int partition(int arr[], int start, int end)
{
 
    int pivot = arr[start];
 
    int count = 0;
    for (int i = start + 1; i <= end; i++) {
        if (arr[i] <= pivot)
            count++;
    }
 
    // Giving pivot element its correct position
    int pivotIndex = start + count;
    swap(arr[pivotIndex], arr[start]);
 
    // Sorting left and right parts of the pivot element
    int i = start, j = end;
 
    while (i < pivotIndex && j > pivotIndex) {
 
        while (arr[i] <= pivot) {
            i++;
        }
 
        while (arr[j] > pivot) {
            j--;
        }
 
        if (i < pivotIndex && j > pivotIndex) {
            swap(arr[i++], arr[j--]);
        }
    }
 
    return pivotIndex;
}

void quickSort(int arr[], int start, int end)
{
 
    // base case
    if (start >= end)
        return;
 
    // partitioning the array
    int p = partition(arr, start, end);
 
    // Sorting the left part
    quickSort(arr, start, p - 1);
 
    // Sorting the right part
    quickSort(arr, p + 1, end);
}

int main(int argc, char* argv[]) {

    if (argc != 2) {
        cout << "Usage: " << argv[0] << " <array_size>" << endl;
        return 1;
    }

    int n = atoi(argv[1]); // Convert command line argument to integer

    if (n <= 0) {
        cout << "Array size must be a positive integer." << endl;
        return 1;
    }

struct timeval start_time, end_time;
    int* arr = new int[n]; // Dynamically allocate array of size 'n'
gettimeofday(&start_time, NULL);
    // Generating random numbers and populating the array
    srand(static_cast<unsigned int>(time(NULL)));
    for (int i = 0; i < n; i++) {
        arr[i] = rand() % 100; // Generates random numbers between 0 and 99
    }

 
    // Sorting the array using Quick Sort algorithm
    quickSort(arr, 0, n - 1);


  

    // Displaying the sorted array
    cout << "Sorted array:" << endl;
    for (int i = 0; i < n; i++) {
        cout << arr[i] << " ";
    }
    cout << endl;

     gettimeofday(&end_time, NULL);
     double execution_time = (end_time.tv_sec - start_time.tv_sec) + (end_time.tv_usec - start_time.tv_usec) / 1000000.0;
 cout << "\nExecution time: " << execution_time << " seconds" << endl;
    delete[] arr; // Free allocated memory

    return 0;
}