#include<iostream>
#include<cmath>
using std::cout;
using std::endl;

static const double pi = 3.14159265358979323846;
const double L = 4 * pi;

#define TESTFLAG 0 

// As a test, made these and ceps macros
// did not change speed substantially
double const eps = 1e-1;
const double epsoverLsq = eps*eps/(L*L);
const double Linv = 1.0/L;

//const double ceps = sqrt(1.0);//0.5 * sqrt(1 + 4 * epsoverLsq);
__device__ inline double kernelp(double x, double y, double3 kernel_params){
    double z = (x - y) * kernel_params.x;
    z -= round(z);
//    z -= (z>0.5) - (z<-0.5); No speed difference
    return 0.5 * z * kernel_params.z * rsqrt( z*z + kernel_params.y  ) - z;
// Below two lines might be slightly faster/more accurate than the above one (compiler might already be doing this anyways)
//    return __fma_rz( 0.5*sqrt( __fma_rz(4.0, epsoverLsq, 1.0)  ), z*rsqrt( __fma_rz(z, z, epsoverLsq)  ), -z); 
//    return fma( 0.5*sqrt(fma(4.0, epsoverLsq, 1.0)  ), z*rsqrt(fma(z, z, epsoverLsq)), -z); 

}

double kernels(double x, double y){
    
//    double const eps = 1e-1;
    double z = (x - y)/L;
    z = z - round(z);
    return 0.5 * z * sqrt(1 + 4 * eps * eps/(L*L)) * rsqrt( z*z + eps*eps/(L*L)  ) - z;
    
//    return x*y;
}

//////////////////////////////
// Dynamic parallel version //
//////////////////////////////
// TODO Replace the atomicAdd with a reduce sum over an array, should be a bit faster (less serialization)
__global__ void direct_e_particle(double *d_efield_p, double *d_particles, double target_loc, double *d_weights, size_t source_size, double3 kernel_params){
    int idx = blockIdx.x*blockDim.x + threadIdx.x;

    if (idx >= source_size){return;}

    double local_eval = kernelp(target_loc, d_particles[idx], kernel_params) * d_weights[idx];

    atomicAdd(d_efield_p, local_eval);
}

__global__ void direct_e_sum_dynamic(double *d_efield, double *d_particles, double *d_target, double *d_weights,
        size_t source_size, size_t target_size, double3 kernel_params){

    int idx = blockIdx.x*blockDim.x + threadIdx.x;

   if (idx >= target_size){return;}
   // Initilize to zero
   d_efield[idx] = 0.0;

   double target_loc = d_target[idx];

   int blocksize = 1024;
   int gridlen = (source_size + blocksize - 1) / blocksize;
   direct_e_particle<<<gridlen,blocksize>>>(d_efield + idx, d_particles, target_loc, d_weights, source_size, kernel_params);
}

//////////////////////////////////
// Non-dynamic parallel version //
//////////////////////////////////

__global__ void direct_e_sum(double *d_efield, double *d_particles, double *d_target, double *d_weights,
        const size_t source_size, size_t target_size, double3 kernel_params){

   const int idx = blockIdx.x*blockDim.x + threadIdx.x;

   double local_e = 0.0;

   if (idx >= target_size){return;}

   const double target_loc = d_target[idx];

   for (size_t k=0; k<source_size;k++){
       local_e += kernelp(target_loc, d_particles[k], kernel_params) * d_weights[k];
   }

   d_efield[idx] = local_e;
}

/*
__global__ void direct_e_sum(double *d_efield, double *d_particles, double *d_target, double *d_weights,
        const size_t source_size, size_t target_size){

   const int idx = blockIdx.x*blockDim.x + threadIdx.x;

   double local_e = 0.0;

   const double target_loc = idx < target_size ? d_target[idx] : 0.0;

   __shared__ double source_points[128];
   __shared__ double weights[128];

   for(size_t j=0;j<source_size/128;j++){
       source_points[threadIdx.x] = d_particles[128*j + threadIdx.x];
       weights[threadIdx.x] = d_weights[128*j + threadIdx.x];
       __syncthreads();

       if (idx < target_size){
           for (size_t k=0; k<128;k++){
               local_e += kernelp(target_loc, source_points[k]) * weights[k];
               //local_e += kernelp(target_loc, source_points[k]) * d_weights[k];
           }
       }
   }

   if (idx < target_size){

       //for(size_t j=source_size/128;j<source_size;j++){
       //    local_e += kernelp(target_loc, source_points[j]) * weights[j];
       //}

       d_efield[idx] = local_e;
   }
}
*/
// No significant speedup over non-nested
__global__ void direct_e_sum_nested(double *d_efield, double *d_particles, double *d_target, double *d_weights,
        const size_t source_size, size_t target_size, double3 kernel_params){

    const int idx_x = blockIdx.x*blockDim.x + threadIdx.x;
    const int idx_y = blockIdx.y*blockDim.y + threadIdx.y;

    if ( (idx_x >= target_size) or (idx_y >= source_size)){return;}

    double local_target = d_particles[idx_x];
    double local_source = d_particles[idx_y];

    double local_eval = kernelp(local_target, local_source, kernel_params) * d_weights[idx_y];

    atomicAdd(d_efield + idx_x, local_eval);

   // printf("idx_x = %d\t idx_y = %d\t local_e = %.5f\t e = %.5f\n", idx_x, idx_y, local_eval, d_efield[idx_x]);
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
        size_t source_size, size_t target_size, bool nested){


    // Initilize output array to zero
    for(size_t k=0;k<target_size;k++){
        e_field[k] = 0.0;
    }

    double3 kernel_params;
    kernel_params.x = Linv;
    kernel_params.y = epsoverLsq;
    kernel_params.z = sqrt(1.0 + 4.0 * epsoverLsq);


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

    errcode = cudaMemset(d_efield, 0.0, target_size*sizeof(double));
    if (errcode != cudaSuccess){
         cout << "Failed to initilize e_field to 0.0 with code " << errcode << " " << cudaGetErrorString(errcode) <<endl;
    }
    
    int blocksize = 128;
    int gridlen = target_size / blocksize;
    if (target_size % blocksize != 0) {gridlen++;}


    // Call the kernel
    if(!nested){
        //cudaDeviceSetLimit(cudaLimitDevRuntimePendingLaunchCount, 32768);
        //direct_e_sum_dynamic<<<gridlen,blocksize>>>(d_efield, d_particles, d_target, d_weights, source_size, target_size);
        direct_e_sum<<<gridlen,blocksize>>>(d_efield, d_particles, d_target, d_weights, source_size, target_size, kernel_params);

    }
    else{
        int blocksize = 16;
        int gridlen = target_size / blocksize;
        if (target_size % blocksize != 0) {gridlen++;}

        dim3 grid(gridlen, gridlen, 1);
        dim3 blockdim(blocksize, blocksize, 1);
        direct_e_sum_nested<<<grid,blockdim>>>(d_efield, d_particles, d_target, d_weights, source_size, target_size, kernel_params);
    }
    errcode = cudaDeviceSynchronize();
    if (errcode != cudaSuccess){
        cout << "Kernel launch failed with code " << errcode << " " << cudaGetErrorString(errcode) << endl;
    }


    errcode = cudaMemcpy(e_field, d_efield, target_size*sizeof(double), cudaMemcpyDeviceToHost);
    if (errcode != cudaSuccess){
        cout << "Failed to transfer calculated e_field to host with code " << errcode << " " << cudaGetErrorString(errcode) << endl;
    }

#if TESTFLAG
    for (size_t k=0;k<source_size;k++){
        cout << "e_direct[" << k << "] = " << e_field[k] << endl;
    }
#endif


    cudaFree(d_efield);
    cudaFree(d_target);
    cudaFree(d_weights);
    cudaFree(d_particles);

}
