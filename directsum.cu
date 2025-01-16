#include "quicksort.h"
#include<iostream>
#include<cmath>
using std::cout;
using std::endl;


__device__ double kernelp(double x, double y){
/*    double const eps = 1e-8;
    double z = x - y;
    z = z - round(z);
    return 0.5 * z * sqrt(1 + 4 * eps * eps) / sqrt( z*z + eps*eps  ) - z;
    */
    return x*y;
}

double kernels(double x, double y){
    /*
    double const eps = 1e-8;
    double z = x - y;
    z = z - round(z);
    return 0.5 * z * sqrt(1 + 4 * eps * eps) / sqrt( z*z + eps*eps  ) - z;
    */
    return x*y;
}

//////////////////////////////
// Dynamic parallel version //
//////////////////////////////
// TODO Replace the atomicAdd with a reduce sum over an array, should be a bit faster (less serialization)
__global__ void direct_e_particle(double *d_efield_p, double *d_particles, double target_loc, double *d_weights, size_t source_size){
    int idx = blockIdx.x*blockDim.x + threadIdx.x;

    if (idx >= source_size){return;}

    double local_eval = kernelp(target_loc, d_particles[idx]) * d_weights[idx];

    atomicAdd(d_efield_p, local_eval);
}

__global__ void direct_e_sum_dynamic(double *d_efield, double *d_particles, double *d_target, double *d_weights,
        size_t source_size, size_t target_size){

    int idx = blockIdx.x*blockDim.x + threadIdx.x;

   if (idx >= target_size){return;}
   // Initilize to zero
   d_efield[idx] = 0.0;

   double target_loc = d_target[idx];

   int blocksize = 1024;
   int gridlen = (source_size + blocksize - 1) / blocksize;
   direct_e_particle<<<gridlen,blocksize>>>(d_efield + idx, d_particles, target_loc, d_weights, source_size);
}

//////////////////////////////////
// Non-dynamic parallel version //
//////////////////////////////////
__global__ void direct_e_sum(double *d_efield, double *d_particles, double *d_target, double *d_weights,
        size_t source_size, size_t target_size){

   int idx = blockIdx.x*blockDim.x + threadIdx.x;

   if (idx >= target_size){return;}

   double target_loc = d_target[idx];
   double local_e = 0.0;

   for (size_t k=0; k<source_size;k++){
       local_e += kernelp(target_loc, d_particles[k]) * d_weights[k];
   }
   d_efield[idx] = local_e;
}

void directsum_serial(double *e_field, double *source_particles, double *target_particles, double *weights,
        size_t source_size, size_t target_size){

    for(size_t k=0;k<target_size;k++){
        e_field[k] = 0.0;
    }


    for (size_t k=0;k<target_size;k++){
        for (size_t j=0;j<source_size;j++){
            e_field[k] += kernels(target_particles[k], source_particles[j])*weights[j];
        }
    }


}

void directsum(double *e_field, double *source_particles, double *target_particles, double *weights,
        size_t source_size, size_t target_size, bool dynamic){


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

    
    int blocksize = 128;
    int gridlen = target_size / blocksize;
    if (target_size % blocksize != 0) {gridlen++;}


    // Call the kernel
    if(dynamic){
	cudaDeviceSetLimit(cudaLimitDevRuntimePendingLaunchCount, 32768);
        direct_e_sum_dynamic<<<gridlen,blocksize>>>(d_efield, d_particles, d_target, d_weights, source_size, target_size);
    }
    else{
        direct_e_sum<<<gridlen,blocksize>>>(d_efield, d_particles, d_target, d_weights, source_size, target_size);
    }
    errcode = cudaDeviceSynchronize();
    if (errcode != cudaSuccess){
        cout << "Kernel launch failed with code " << errcode << " " << cudaGetErrorString(errcode) << endl;
    }


    errcode = cudaMemcpy(e_field, d_efield, target_size*sizeof(double), cudaMemcpyDeviceToHost);
    if (errcode != cudaSuccess){
        cout << "Failed to transfer calculated e_field to host with code " << errcode << " " << cudaGetErrorString(errcode) << endl;
    }

    cudaFree(d_efield);
    cudaFree(d_target);
    cudaFree(d_weights);
    cudaFree(d_particles);

}
