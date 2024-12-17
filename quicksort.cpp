/*Written with reference to: 
https://www.geeksforgeeks.org/cpp-program-for-quicksort/
https://www.softwaretestinghelp.com/quick-sort/
https://opendsa-server.cs.vt.edu/embed/quicksortAV
Github user gauravjain2's quicksort implementation https://github.com/gauravjain2/quicksort/blob/master/quicksort.cpp
Cracking the Coding Interview: 6th Edition pg 146-148*/

#include "quicksort.h"
#include <iostream>

int divideSortPivot(double arr[], int left, int right, size_t indicies[]){
    int pivotUnsorted = ((right-left)/2)+left;
    double pivot = arr[pivotUnsorted];//pivot from middle of array; right - left gives us the range, divide than range in half, then add that to the start position(left)
    int i=-1;
    for (int j=left;j<=right;++j){//loop finds true index of pivot to use for comparison in next loop(left/right comparison swaps) - affects time complexity as entire data structure must be iterated through
        if(arr[j]<=pivot){
            ++i;
        }
    }
    int pivotIndex = i+left;//amount of elements less than/equal to pivot, plus index of where we began searching
    i=left-1;
    for(int j=left;j<=right;++j){
        if(arr[j]<=pivot){
            ++i;
            if(j==pivotUnsorted)
                pivotUnsorted=i;
            std::swap(arr[i],arr[j]);
            std::swap(indicies[i],indicies[j]);
        }
    }
    std::swap(arr[pivotIndex],arr[pivotUnsorted]);
    std::swap(indicies[pivotIndex],indicies[pivotUnsorted]);
    return (pivotIndex);
}

void quicksort(double arr[], int left, int right, size_t indicies[]){
    if(left >= right)//base case: if leftmost index of "sub" array we are recursively sorting is greater than rightmost, then return.
        return;

    int pivot = divideSortPivot(arr, left, right, indicies);

    quicksort(arr, left, pivot-1, indicies);
    quicksort(arr, pivot+1, right, indicies);
}

void quicksort(double arr[],int arrLength, size_t indicies[]){
    quicksort(arr, 0, arrLength-1, indicies);
}
