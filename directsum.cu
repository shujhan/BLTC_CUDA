#include<iostream>
using std::cout;
using std::endl;

#define TIMEFLAG 0 // For extra timing output

#if TIMEFLAG
#include<ctime>
#endif

__device__ void kernel(double *res, double x, double y){
    res[0] = x*y;
}

double kernel_serial(double x, double y){
    return x*y;
}

__global__ void direct_e_sum(double *d_efield, double *d_particles, double *d_target, double *d_weights,
        size_t source_size, size_t target_size){

   int idx = blockIdx.x*blockDim.x + threadIdx.x;

   if (idx >= target_size){return;}

   double target_loc = d_target[idx];
   double kernel_res;

   for (size_t k=0; k<source_size;k++){
       kernel(&kernel_res, target_loc, d_particles[k]);
       d_efield[idx] += kernel_res * d_weights[k];
   }
}

void directsum_serial(double *e_field, double *source_particles, double *target_particles, double *weights,
        size_t source_size, size_t target_size){

    for(size_t k=0;k<target_size;k++){
        e_field[k] = 0.0;
    }

#if TIMEFLAG
    const std::clock_t start = std::clock();
#endif

    for (size_t k=0;k<target_size;k++){
        for (size_t j=0;j<source_size;j++){
            e_field[k] += kernel_serial(target_particles[k], source_particles[j])*weights[j];
        }
    }

#if TIMEFLAG
    const std::clock_t end = std::clock();

    cout << "Direct sum serial time (ms): " << 1000 * double(end-start)/CLOCKS_PER_SEC << endl;
#endif

}

void directsum(double *e_field, double *source_particles, double *target_particles, double *weights,
        size_t source_size, size_t target_size){

    // Initilize output array to zero
    for(size_t k=0;k<target_size;k++){
        e_field[k] = 0.0;
    }


    // Allocate and transfer data to GPU
   double *d_particles;
   double *d_target;
   double *d_weights;
   double *d_efield;

   cudaError_t errcode;
    errcode = cudaMalloc(&d_particles, source_size*sizeof(double));
    if (errcode != cudaSuccess){
         cout << "Failed to allocate source particles on device with code " << errcode << " " << cudaGetErrorString(errcode) <<  endl;
    }
    errcode = cudaMalloc(&d_weights, source_size*sizeof(double));
    if (errcode != cudaSuccess){
         cout << "Failed to allocate weights on device with code " << errcode << " " << cudaGetErrorString(errcode) <<endl;
    }
    errcode = cudaMalloc(&d_target, target_size*sizeof(double));
    if (errcode != cudaSuccess){
         cout << "Failed to allocate target particles on device with code " << errcode << " " << cudaGetErrorString(errcode) <<  endl;
    }
    errcode = cudaMalloc(&d_efield, target_size*sizeof(double));
    if (errcode != cudaSuccess){
         cout << "Failed to allocate e_field on device with code " << errcode << " " << cudaGetErrorString(errcode) <<endl;
    }

    errcode = cudaMemcpy(d_particles, source_particles, source_size*sizeof(double), cudaMemcpyHostToDevice);
    if (errcode != cudaSuccess){
         cout << "Failed to transfer source particles to device with code " << errcode << " " << cudaGetErrorString(errcode) <<endl;
    }
    errcode = cudaMemcpy(d_target, target_particles, target_size*sizeof(double), cudaMemcpyHostToDevice);
    if (errcode != cudaSuccess){
         cout << "Failed to transfer target particles to device with code " << errcode << " " << cudaGetErrorString(errcode) <<endl;
    }
 
    errcode = cudaMemcpy(d_weights, weights, source_size*sizeof(double), cudaMemcpyHostToDevice);
    if (errcode != cudaSuccess){
         cout << "Failed to transfer particle weights to device with code " << errcode << " " << cudaGetErrorString(errcode) <<endl;
    }

    
    int blocksize = 1024;
    int gridlen = target_size / blocksize;
    if (target_size % blocksize != 0) {gridlen++;}

#if TIMEFLAG
    float elapsed_time = 0;
    cudaEvent_t start,stop;
    cudaEventCreate(&start);
    cudaEventCreate(&stop);

    cudaEventRecord(start);
#endif

    // Call the kernel
    direct_e_sum<<<gridlen,blocksize>>>(d_efield, d_particles, d_target, d_weights, source_size, target_size);
    errcode = cudaDeviceSynchronize();
    if (errcode != cudaSuccess){
        cout << "Kernel launch failed with code " << errcode << " " << cudaGetErrorString(errcode) << endl;
    }

#if TIMEFLAG
    cudaEventRecord(stop);
    cudaEventSynchronize(stop);
    cudaEventElapsedTime(&elapsed_time, start, stop);

    cudaEventDestroy(start);
    cudaEventDestroy(stop);

    cout << "Direct summation parallel time (ms): " << elapsed_time << endl;
#endif

    errcode = cudaMemcpy(e_field, d_efield, target_size*sizeof(double), cudaMemcpyDeviceToHost);
    if (errcode != cudaSuccess){
        cout << "Failed to transfer calculated e_field to host with code " << errcode << " " << cudaGetErrorString(errcode) << endl;
    }

    cudaFree(d_efield);
    cudaFree(d_target);
    cudaFree(d_weights);
    cudaFree(d_particles);

}
