#include "FieldStructure.hpp"
#include <cmath>
#include <cassert>
#include <cstring>
#include <sys/times.h>
// #include <openacc.h>
#include <iostream>
#include <cfloat> // dbl_min
#include <cstddef>
using namespace std;
#if OPENACC_ENABLED
#include <accelmath.h>
#endif

ElectricField::~ElectricField() = default;


E_MQ_DirectSum::E_MQ_DirectSum() {}
E_MQ_DirectSum::E_MQ_DirectSum(double L, double epsilon) : L(L), epsilon(epsilon) {}
E_MQ_DirectSum::~E_MQ_DirectSum() = default;

void E_MQ_DirectSum::operator() (double *e_field, double *source_particles, double *target_particles, double *weights,
        size_t source_size, size_t target_size)
{    
    // Initilize output array to zero
    for(size_t k=0;k<target_size;k++){
        e_field[k] = 0.0;
    }

    double3 kernel_params;
    kernel_params.x = 1.0 / L;
    kernel_params.y = epsilon * epsilon / (L * L);
    kernel_params.z = sqrt(1.0 + 4.0 * kernel_params.y);


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


    
    direct_e_sum<<<gridlen,blocksize>>>(d_efield, d_particles, d_target, d_weights, source_size, target_size, kernel_params);
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
