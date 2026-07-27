#include <stdio.h>
#include <stdlib.h>
#include <limits.h>
#include <time.h>


void swap(int *a, int *b) 
{
    int temp = *a;
    *a = *b;
    *b = temp;
}

// ============
// 1. MERGESORT 
// ============
void merge(int A[], int p, int q, int r) 
{
    int n1 = q - p + 1;
    int n2 = r - q;

    int *L = (int *)malloc((n1 + 1) * sizeof(int));
    int *R = (int *)malloc((n2 + 1) * sizeof(int));

    for (int i = 0; i < n1; i++) 
    {
        L[i] = A[p + i];
    }
    for (int j = 0; j < n2; j++) 
    {
        R[j] = A[q + 1 + j];
    }

    
    L[n1] = INT_MAX;
    R[n2] = INT_MAX;

    int i = 0;
    int j = 0;
    for (int k = p; k <= r; k++) 
    {
        if (L[i] <= R[j]) 
        {
            A[k] = L[i];
            i++;
        } 
        else 
        {
            A[k] = R[j];
            j++;
        }
    }

    free(L);
    free(R);
}

void merge_sort(int A[], int p, int r)
{
    if (p < r) 
    {
        int q = p + (r - p) / 2;
        merge_sort(A, p, q);
        merge_sort(A, q + 1, r);
        merge(A, p, q, r);
    }
}

// ===========
// 2. HEAPSORT 
// ===========
void max_heapify(int A[], int i, int heap_size) 
{
    int l = 2 * i + 1; 
    int r = 2 * i + 2; 
    int largest = i;

    if (l < heap_size && A[l] > A[largest])
    {
        largest = l;
    }
    if (r < heap_size && A[r] > A[largest]) 
    {
        largest = r;
    }
    if (largest != i)
    {
        swap(&A[i], &A[largest]);
        max_heapify(A, largest, heap_size);
    }
}

void build_max_heap(int A[], int n) 
{
    for (int i = n / 2 - 1; i >= 0; i--)
    {
        max_heapify(A, i, n);
    }
}

void heapsort(int A[], int n) 
{
    build_max_heap(A, n);
    int heap_size = n;
    for (int i = n - 1; i >= 1; i--) 
    {
        swap(&A[0], &A[i]);
        heap_size--;
        max_heapify(A, 0, heap_size);
    }
}

// ============
// 3. QUICKSORT 
// ============
int partition(int A[], int p, int r) 
{
    int x = A[r]; 
    int i = p - 1;

    for (int j = p; j <= r - 1; j++) 
    {
        if (A[j] <= x) 
        {
            i++;
            swap(&A[i], &A[j]);
        }
    }
    swap(&A[i + 1], &A[r]);
    return i + 1;
}

void quicksort(int A[], int p, int r) 
{
    if (p < r) 
    {
        int q = partition(A, p, r);
        quicksort(A, p, q - 1);
        quicksort(A, q + 1, r);
    }
}

//helper function
void copy_array(const int src[], int dest[], int n) 
{
    for (int i = 0; i < n; i++) {
        dest[i] = src[i];
    }
}

int main(void) 
{
    srand((unsigned int)time(NULL));

    int iterations = 100;
    
    printf("===============================================================\n");
    printf(" Size (N) | MergeSort (ms) | HeapSort (ms) | QuickSort (ms)    \n");
    printf("===============================================================\n");

    for (int size = 100; size <= 1000; size += 100) 
    {
        double total_merge_time = 0.0;
        double total_heap_time  = 0.0;
        double total_quick_time = 0.0;

        int *base_arr  = (int *)malloc(size * sizeof(int));
        int *work_arr  = (int *)malloc(size * sizeof(int));

        for (int iter = 0; iter < iterations; iter++) 
        {
            // Generate random data
            for (int k = 0; k < size; k++) 
            {
                base_arr[k] = rand() % 10000;
            }

            clock_t start, end;

            // MergeSort Time
            copy_array(base_arr, work_arr, size);
            start = clock();
            merge_sort(work_arr, 0, size - 1);
            end = clock();
            total_merge_time += ((double)(end - start)) / CLOCKS_PER_SEC;

            // HeapSort Time
            copy_array(base_arr, work_arr, size);
            start = clock();
            heapsort(work_arr, size);
            end = clock();
            total_heap_time += ((double)(end - start)) / CLOCKS_PER_SEC;

            // QuickSort Time
            copy_array(base_arr, work_arr, size);
            start = clock();
            quicksort(work_arr, 0, size - 1);
            end = clock();
            total_quick_time += ((double)(end - start)) / CLOCKS_PER_SEC;
        }

        // Calculate average time in milliseconds
        double avg_merge = (total_merge_time / iterations) * 1000.0;
        double avg_heap  = (total_heap_time / iterations)  * 1000.0;
        double avg_quick = (total_quick_time / iterations) * 1000.0;

        printf(" %8d | %14.4f | %13.4f | %14.4f \n", 
               size, avg_merge, avg_heap, avg_quick);

        free(base_arr);
        free(work_arr);
    }

    printf("===============================================================\n");
    return 0;
}
