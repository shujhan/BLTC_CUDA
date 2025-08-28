#include "FieldStructure.hpp"

#include<cmath>
#include <iostream>
#include <iomanip>
#include <assert.h>
#include <vector>

using namespace std;

ElectricField::~ElectricField() = default;

//////////////////////////////////
/////////// Direct Sum /////////////
//////////////////////////////////

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

void E_MQ_DirectSum::__global__ void direct_e_sum(double *d_efield, double *d_particles, double *d_target, double *d_weights,
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



void E_MQ_DirectSum::__device__ inline double kernelp(double x, double y, double3 kernel_params){
     double z = (x - y) * kernel_params.x;
     z -= round(z);
     return 0.5 * z * kernel_params.z * rsqrt( z*z + kernel_params.y  ) - z;
}


//////////////////////////////////
/////////// Treecode /////////////
//////////////////////////////////


E_MQ_Treecode::E_MQ_Treecode() {}
E_MQ_Treecode::E_MQ_Treecode(double L, double epsilon, double beta) : 
L(L), epsilon(epsilon), mac(mac), beta(beta) {}

// Definitions for E_MQ_Treecode
E_MQ_Treecode::E_MQ_Treecode(double L, double epsilon,
                             double mac, int degree, int max_source, int max_target,
                             int verbosity)
    : L(L), epsilon(epsilon), mac(mac), degree(degree),
      max_source(max_source), max_target(max_target), verbosity(verbosity),
      lambda(nullptr), particles_x(nullptr), keps_tc_reord(nullptr), keps_tc_noreord(nullptr), 
      tree_members{nullptr, nullptr}, leaf_members{nullptr, nullptr}, interaction_list_far(nullptr), 
      interaction_list_far_size(nullptr), max_far_size(0), interaction_list_near(nullptr), interaction_list_near_size(nullptr), max_near_size(0), 
      cluster_list_t1(nullptr), cluster_list_moments(nullptr) {
        #ifdef OPENACC_ENABLED
        #pragma acc enter data copyin(this)
        #endif
      }


    E_MQ_Treecode::~E_MQ_Treecode() 
    {
    #ifdef OPENACC_ENABLED
      #pragma acc exit data delete(this)
    #endif
    };


// use BLTC function 
void E_MQ_Treecode::operator()(double *e_field, double *source_particles, double *target_particles, double *weights, 
        size_t e_field_size, size_t source_size, size_t target_size){
    // Initilize device vars
    panel *d_tree_list;
    double *d_particles;
    double *d_targets;
    double *d_weights;
    double *d_efield;
    int *d_near_list;
    int *d_far_list;
    int *d_leaf_indicies;

    cudaError_t errcode;


    errcode = cudaMalloc(&d_particles, source_size*sizeof(double));
    if(checkcudaerr(errcode) !=0){cout << "Failed allocating source particles" << endl;}

    errcode = cudaMalloc(&d_targets, target_size*sizeof(double));
    if(checkcudaerr(errcode) !=0){cout << "Failed allocating target particles" << endl;}

    errcode = cudaMalloc(&d_weights, source_size*sizeof(double));
    if(checkcudaerr(errcode) != 0){cout << "Failed allocating source weights" << endl;}

    errcode = cudaMalloc(&d_efield, target_size*sizeof(double));
    if(checkcudaerr(errcode) != 0){cout << "Failed allocating e_field" << endl;}
    errcode = cudaMemset(d_efield, 0.0, target_size*sizeof(double));
    if(checkcudaerr(errcode) != 0){cout << "Failed initilizing e_field to 0" << endl;}


#if TESTFLAG
    cout << "Input particles:" << endl;
    for(size_t k=0;k<source_size;k++){
        cout << "x[" << k << "] = " << source_particles[k] << endl;
    }
#endif

    

    // Set up root panel
    panel root;
    double xmin = source_particles[0];
    double xmax = source_particles[0];
    root.members[0] = 0;
    root.members[1] = source_size-1;
    for (size_t k=0; k<source_size; k++){
        if(source_particles[k] < xmin){ xmin = source_particles[k]; }
        if(source_particles[k] > xmax){ xmax = source_particles[k]; }
    }
    root.xinterval[0] = xmin-0.001;
    root.xinterval[1] = xmax+0.001;
    root.xc = (root.xinterval[0] + root.xinterval[1])/2;
    root.level = 0;
    root.num_members = source_size;
    for (int k=0;k<PP;k++){
        root.s[k] = 0.5 * ( root.xinterval[0] + root.xinterval[1] + std::cos(k*pi/PP)*( root.xinterval[1] - root.xinterval[0]  ));
    }
    
    int tree_size = 1;
    int leaf_size = 0;

    // TODO Need to sort particles and also re-sort the corresponding weights
    size_t* source_indicies = (size_t*)malloc(sizeof(size_t)*source_size);
    for(size_t k=0;k<source_size;k++){
        source_indicies[k] = k;
    }

    if (source_size > N0){
        split_panel(&root, source_particles, &tree_size, &leaf_size, source_indicies);
    }
    else{ 
        // TODO Should just abort to direct sum here
        root.left_child = NULL;
        root.right_child = NULL;
        leaf_size=1;  
    }

    cout << "Set up root panel" << endl;


    double *sorted_particles = (double*)malloc(sizeof(double) * source_size);
    double *sorted_weights = (double*)malloc(sizeof(double) * source_size);
    for(size_t i=0;i<source_size;i++){
        sorted_particles[i] = source_particles[source_indicies[i]];
        sorted_weights[i] = weights[source_indicies[i]];
    }

#if TESTFLAG
    cout << "Sorted indicies:" << endl;
        for(size_t k=0;k<source_size;k++){
            cout << source_indicies[k] << endl;
        }
    cout << endl << "Sorted particles:" << endl << endl;;
    for(size_t k=0;k<source_size;k++){
        cout << "x[" << k << "] = " << sorted_particles[k] << endl;
    }
#endif


    errcode = cudaMemcpy(d_particles, sorted_particles, source_size*sizeof(double), cudaMemcpyHostToDevice);
    if(checkcudaerr(errcode) != 0){cout << "Failed copying source particles to device" << endl;}

    errcode = cudaMemcpy(d_targets, target_particles, target_size*sizeof(double), cudaMemcpyHostToDevice);
    if(checkcudaerr(errcode) != 0){cout << "Failed copying target particles to device" << endl;}

    errcode = cudaMemcpy(d_weights, sorted_weights, source_size*sizeof(double), cudaMemcpyHostToDevice);
    if(checkcudaerr(errcode) != 0){cout << "Failed copying source weights to device" << endl;}

    //free(sorted_particles);
    //free(sorted_weights);


    panel tree_list[tree_size];
    int leaf_indicies[leaf_size];

    int id = 0;
    int leaf_id = 0;

    cout << "Initilizing tree list" << endl;

    init_tree_list(&root, tree_list, &id, leaf_indicies, &leaf_id);

#if TESTFLAG
    cout << endl;
    cout << "Tree list is" << endl;
    for (int k=0;k<tree_size;k++){
        cout << "Range: [" << tree_list[k].xinterval[0] << "," << tree_list[k].xinterval[1] << "]\t Particles: [" << tree_list[k].members[0] << "," << tree_list[k].members[1] <<  "]\t ID: " << tree_list[k].id << endl;
    }
    cout << endl;
    cout << "Leaf indicies are" << endl;
    for (int k=0;k<leaf_size;k++){
        cout << "Range: [" << tree_list[leaf_indicies[k]].xinterval[0] << "," << tree_list[leaf_indicies[k]].xinterval[1] << "]\t ID: " << leaf_indicies[k] << endl;
    }
    cout << endl;

    cout << "Tree list length: " << tree_size << endl;
    cout << "Number of leafs: " << leaf_size << endl;
    cout << endl;
#endif

    cout << "Initilized tree list" << endl;

    int near_interactions[leaf_size*leaf_size];
    int far_interactions[leaf_size*leaf_size];
    for (size_t k=0;k<leaf_size*leaf_size;k++){
        near_interactions[k] = -1;
        far_interactions[k] = -1;
    }

    int total_nears = 0;
    for(int k=0;k<leaf_size;k++){
        int near_index = 0;
        int far_index = 0;
        init_interaction_lists(tree_list+leaf_indicies[k], &root, near_interactions, far_interactions, &near_index, &far_index, k, leaf_size, L, &total_nears); 
    }

    cout << "Initilized interaction lists" << endl;

#if TESTFLAG
    cout << endl;
    for(int k=0;k<leaf_size;k++){
        cout << "Leaf " << k << " [" << tree_list[leaf_indicies[k]].xinterval[0] << "," << tree_list[leaf_indicies[k]].xinterval[1] << "] has ID " << leaf_indicies[k] << " and near Ids" << endl;
        for(int i=0;i<leaf_size;i++){
            if (near_interactions[k*leaf_size + i] == -1){break;}
            cout << near_interactions[k*leaf_size + i] << "\t";
        }
        cout << endl;
        for(int i=0;i<leaf_size;i++){
            if (near_interactions[k*leaf_size + i] == -1){break;}
            cout << "[" << tree_list[near_interactions[k*leaf_size + i]].xinterval[0] << " , " << tree_list[near_interactions[k*leaf_size+i]].xinterval[1] << "]\t";
        }
        cout << endl;
        cout << "and far Ids" << endl;
        for(int i=0;i<leaf_size;i++){
            if (far_interactions[k*leaf_size + i] == -1){break;}
            cout << far_interactions[k*leaf_size + i] << "\t";
        }
        cout << endl;
        for(int i=0;i<leaf_size;i++){
            if (far_interactions[k*leaf_size + i] == -1){break;}
            cout << "[" << tree_list[far_interactions[k*leaf_size + i]].xinterval[0] << " , " << tree_list[far_interactions[k*leaf_size+i]].xinterval[1] << "]\t";
        }
        cout << endl;
        cout << endl;
    }
    cout << endl;
#endif

    errcode = cudaMalloc(&d_tree_list, tree_size*sizeof(panel));
    if(checkcudaerr(errcode) != 0){cout << "Failed allocating tree list" << endl;}

    
    errcode = cudaMalloc(&d_leaf_indicies, leaf_size*sizeof(int));
    if(checkcudaerr(errcode) != 0){cout << "Failed allocating leaf indicies" << endl;}

    errcode = cudaMemcpy(d_tree_list, tree_list, tree_size*sizeof(panel), cudaMemcpyHostToDevice);
    if(checkcudaerr(errcode) != 0){cout << "Failed copying tree list to device" << endl;}
 
    
    errcode = cudaMemcpy(d_leaf_indicies, leaf_indicies, leaf_size*sizeof(int), cudaMemcpyHostToDevice);
    if(checkcudaerr(errcode) != 0){cout << "Failed copying leaf indicies to device" << endl;}


    errcode = cudaMalloc(&d_near_list, leaf_size*leaf_size*sizeof(int));
    if(checkcudaerr(errcode) != 0){cout << "Failed allocating near interaction list" << endl;}
    
    errcode = cudaMalloc(&d_far_list, leaf_size*leaf_size*sizeof(int));
    if(checkcudaerr(errcode) != 0){cout << "Failed allocating far interaction list" << endl;}

    errcode = cudaMemcpy(d_near_list, near_interactions, leaf_size*leaf_size*sizeof(int), cudaMemcpyHostToDevice);
    if(checkcudaerr(errcode) != 0){cout << "Failed copying near interactions to device" << endl;}

    errcode = cudaMemcpy(d_far_list, far_interactions, leaf_size*leaf_size*sizeof(int), cudaMemcpyHostToDevice);
    if(checkcudaerr(errcode) != 0){cout << "Failed copying far interactions to device" << endl;}

    // TODO This might be better as the smallest divisor of 32 larger than PP (might need to be slightly careful with some warp-level primitives though)
    int blocksize = 32;
    //int gridlen = (PP*tree_size + blocksize - 1) / blocksize;
    int gridlen = (blocksize*tree_size + blocksize - 1) / blocksize;
    init_modified_weights<<<gridlen,blocksize>>>(d_tree_list, d_particles, d_weights, source_size, tree_size);

    cout << "Initilized modified weights" << endl;

#if TESTFLAG
    errcode = cudaMemcpy(tree_list, d_tree_list, tree_size*sizeof(panel), cudaMemcpyDeviceToHost);
    if(errcode != 0){cout << "Failed to copy tree list to host" << endl;}

    for(int k=0;k<tree_size;k++){
        cout << "Modified weights for panel " << k << " at [" << tree_list[k].xinterval[0] << "," << tree_list[k].xinterval[1] << "]" << endl;
        for (int j=0;j<PP;j++){
            cout << tree_list[k].modified_weights[j] << "\t";
        }
        cout << endl;
        cout << "Chebyshev Points for panel " << k << " at [" << tree_list[k].xinterval[0] << "," << tree_list[k].xinterval[1] << "]" << endl;
        for (int j=0;j<PP;j++){
            cout << tree_list[k].s[j] << "\t";
        }
        cout << endl << endl;
    }
    cout << endl;
#endif



    cudaDeviceSynchronize();

    double3 kernel_params;
    kernel_params.x = Linv;
    kernel_params.y = epsoverLsq;
    kernel_params.z = 0.5*ceps;


    // Might be some tuning to be done here
    blocksize = 256;
    
    cudaStream_t streams_far[leaf_size];
    for (int i=0; i<leaf_size;i++){
        cudaStreamCreate(&streams_far[i]);

        gridlen = (tree_list[leaf_indicies[i]].num_members * tree_list[leaf_indicies[i]].far_size + blocksize - 1) / blocksize;
        computepanelsum_far<<<gridlen,blocksize, 0, streams_far[i]>>>(d_efield, d_tree_list + leaf_indicies[i], d_tree_list, d_particles, d_particles, d_far_list, i, leaf_size, kernel_params);
    }

    cout << "Queued far interactions" << endl;

/*
    for (int i=0; i<leaf_size;i++){
            cudaStreamCreate(&streams_near[i]);

            gridlen = (tree_list[leaf_indicies[i]].num_members * tree_list[leaf_indicies[i]].near_size + blocksize - 1) / blocksize;
            computepanelsum_near<<<gridlen,blocksize, 0, streams_near[i]>>>(d_efield, d_tree_list + leaf_indicies[i], d_tree_list, d_targets, d_particles, d_weights, d_near_list, i, leaf_size);
        }
*/

    /*
    int current_idx = 0;
    for(int i=0; i<leaf_size;i++){
        int leaf_particles = tree_list[leaf_indicies[i]].near_size;
        cudaStreamCreate(&streams_near[current_idx]);
        for (int j=0; j<leaf_particles; j++){

            int near_id = near_interactions[i*leaf_size + j];
            size_t left_mem = tree_list[near_id].members[0];
            size_t right_mem = tree_list[near_id].members[1];

            gridlen = (tree_list[leaf_indicies[i]].num_members * (right_mem - left_mem + 1)/16 + blocksize - 1) / blocksize; 

            computepanelsum_near<<<gridlen,blocksize, 0, streams_near[current_idx]>>>(d_efield, d_tree_list + leaf_indicies[i], left_mem, right_mem, d_targets, d_particles, d_weights, leaf_size);

        }
        current_idx += 1;
    }
    */

    /*
    size_t *near_data = (size_t*)malloc(4 * total_nears * sizeof(size_t));
    int interaction_id;
    int running_nears = 0;
    for(size_t i=0;i<leaf_size;i++){
        for(size_t j=0;j<leaf_size;j++){
            interaction_id = near_interactions[i*leaf_size + j];
            if (interaction_id == -1){continue;} // break might be better here

            near_data[4*running_nears] = tree_list[leaf_indicies[i]].members[0];
            near_data[4*running_nears+1] = tree_list[leaf_indicies[j]].members[0];
            near_data[4*running_nears+2] = tree_list[leaf_indicies[i]].num_members;
            near_data[4*running_nears+3] = tree_list[leaf_indicies[j]].num_members;
            running_nears += 1;
        }
    }
    */

    // Interactions are symmetric, so total_nears is always even excluding the self-interaction
    // (just counting the (undirected) edges in a graph where each panel is a node and an edge is an interaction)
//    size_t num_near_interactions = (total_nears - leaf_size) / 2;
    uint4 near_data[total_nears];
    int interaction_id;
    //int leaf_id;
    int running_nears = 0;
    for(size_t i=0;i<leaf_size;i++){
        for(size_t j=0;j<leaf_size;j++){
            interaction_id = near_interactions[i*leaf_size + j];
            if (interaction_id == -1){continue;} // break might be better here
            leaf_id = leaf_indicies[i];
            //if (leaf_indicies[i] > leaf_indicies[j]) {continue;} // Handle symmetry
            //cout << leaf_id << "\t" << interaction_id << endl;
            near_data[running_nears].x = tree_list[leaf_id].members[0];
            near_data[running_nears].y = tree_list[interaction_id].members[0];
            near_data[running_nears].z = tree_list[leaf_id].num_members;
            near_data[running_nears].w = tree_list[interaction_id].num_members;

            running_nears += 1;
        }
    }

    uint4 *d_near_data;
    errcode = cudaMalloc(&d_near_data, total_nears*sizeof(uint4));
    if(checkcudaerr(errcode) != 0){cout << "Failed allocating near data" << endl;}

    errcode = cudaMemcpy(d_near_data, near_data, total_nears*sizeof(uint4), cudaMemcpyHostToDevice);
    if(checkcudaerr(errcode) != 0){cout << "Failed copying near data to device" << endl;}
/* 5/14/25
    blocksize = 512;
    gridlen = (total_nears*N0*N0 + blocksize - 1) / blocksize;
    computepanelsum_near<<<gridlen, blocksize>>>(d_efield, d_tree_list, d_near_list, d_targets, d_particles, d_weights, leaf_size, d_leaf_indicies, d_near_data);
*/

    int xblocklen = 16;
    int yblocklen = 32;
    int zblocklen = 1;
    dim3 blockdim(xblocklen, yblocklen, zblocklen);

    int xgridlen = N0 / xblocklen;
    if (N0 % xblocklen != 0){xgridlen += 1;}


    
    int ygridlen = N0 / yblocklen;
    if (N0 % yblocklen != 0){ygridlen += 1;}


    // Not 100% sure the spillover is handled correctly
//    int ygridlen = (N0/UNROLLNUM) / yblocklen;
//    if (N0 % (UNROLLNUM*yblocklen) != 0){ygridlen += 1;}

    int zgridlen = total_nears / zblocklen;
    if (total_nears % zblocklen != 0){zgridlen += 1;}

    dim3 griddim(xgridlen, ygridlen, zgridlen);
    computepanelsum_near<<<griddim, blockdim>>>(d_efield, d_particles, d_particles, d_weights, d_near_data, total_nears, kernel_params);

    /*
    dim3 blockdim(8,128,1);
    int xgridlen = leaf_size / 8;
    int ygridlen = N0 / 128;
    if (leaf_size %8 != 0){xgridlen += 1;}
    if (N0 %128 != 0){ygridlen += 1;}
    dim3 griddim(xgridlen, ygridlen, 1);
    computepanelsum_near<<<griddim, blockdim>>>(d_efield, d_tree_list, d_near_list, d_targets, d_particles, d_weights, leaf_size, d_leaf_indicies);
*/

    cout << "Queued near interactions" << endl;

        //computepanelsum<<<gridlen,blocksize, 0, streams[i]>>>(d_efield, d_tree_list + leaf_indicies[i], d_tree_list, d_targets, d_particles, d_weights, d_near_list, d_far_list, i, leaf_size);

    cudaDeviceSynchronize();
    for (int i=0; i<leaf_size;i++){
    //    cudaStreamDestroy(streams1[i]);
     //   cudaStreamDestroy(streams2[i]);
     cudaStreamDestroy(streams_far[i]);
//     cudaStreamDestroy(streams_near[i]);
    }

    cout << "Destroyed far streams" << endl;

    // This segfaults, but works if we use current_idx instead
//    for (int i=0; i<current_idx; i++){
//        cudaStreamDestroy(streams_near[i]);
//    }

    cout << "Computed BLTC sum" << endl;

    size_t* d_source_indicies;
    errcode = cudaMalloc(&d_source_indicies, source_size*sizeof(size_t));
    if(checkcudaerr(errcode) !=0){cout << "Failed allocating source indicies" << endl;}
    errcode = cudaMemcpy(d_source_indicies, source_indicies, source_size*sizeof(size_t), cudaMemcpyHostToDevice);
    if(checkcudaerr(errcode) !=0){cout << "Failed to copy source indicies" << endl;}
    
    double* d_ordered_e_field;
    errcode = cudaMalloc(&d_ordered_e_field, source_size*sizeof(double));
    if(checkcudaerr(errcode) !=0){cout << "Failed allocating ordered e_field" << endl;}

    gridlen = (source_size + blocksize - 1) / blocksize;
    re_order<<<gridlen, blocksize>>>(d_ordered_e_field, d_efield, d_source_indicies, source_size);

    errcode = cudaMemcpy(e_field, d_ordered_e_field, target_size*sizeof(double), cudaMemcpyDeviceToHost);
    if(checkcudaerr(errcode) != 0){cout << "Failed to copy e_field to host" << endl;}


#if TESTFLAG
    for (size_t k=0;k<source_size;k++){
        cout << "e[" << k << "] = " << e_field[k] << endl;
    }
#endif


    if (root.left_child){
        free_tree_list(root.left_child);
    }
    if (root.right_child){
        free_tree_list(root.right_child);
    }

    free(sorted_particles);
    free(sorted_weights);

    cudaFree(d_tree_list);
    cudaFree(d_particles);
    cudaFree(d_weights);
    cudaFree(d_leaf_indicies);
    cudaFree(d_efield);
    cudaFree(d_ordered_e_field);
    cudaFree(d_source_indicies);

}

void E_MQ_Treecode::split_panel(panel *p, double *source_particles, int *tree_size, int *leaf_size, size_t *indicies){
    panel *left_child = new panel();
    panel *right_child = new panel();

    // Set left child characteristics
    left_child->xinterval[0] = p->xinterval[0];
    left_child->xinterval[1] = p->xc;
    left_child->xc = 0.5 * (p->xinterval[0] + p->xc);
    left_child->level = p->level + 1;
    left_child->parent = p;
    for (int k=0;k<PP;k++){
        left_child->s[k] = 0.5 * ( left_child->xinterval[0] + left_child->xinterval[1] + std::cos(k*pi/P)*( left_child->xinterval[1] - left_child->xinterval[0]  ));
    }

    // Set right child characteristics
    right_child->xinterval[0] = p->xc;
    right_child->xinterval[1] = p->xinterval[1];
    right_child->xc = 0.5 * (p->xinterval[1] + p->xc);
    right_child->level = p->level + 1;
    right_child->parent = p;
    for (int k=0;k<PP;k++){
        right_child->s[k] = 0.5 * ( right_child->xinterval[0] + right_child->xinterval[1] + std::cos(k*pi/P)*( right_child->xinterval[1] - right_child->xinterval[0]  ));
    }

    std::vector<size_t> left_indicies;
    std::vector<size_t> right_indicies;

    // indicies maps between the original indicies and the (partially) sorted indicies
    // indicies[j] = k means that the jth particle in the new index scheme was the kth particle in the original
    for(size_t i=p->members[0];i<=p->members[1];i++){
        if(source_particles[indicies[i]] <= p->xc){
            left_indicies.push_back(indicies[i]);
        }
        else{
            right_indicies.push_back(indicies[i]);
        }
    }
    
    left_child->num_members = left_indicies.size();
    right_child->num_members = right_indicies.size();

    if(left_child->num_members > 0){
        p->left_child = left_child;
        *tree_size += 1;
    }
    else{
        p->left_child = NULL;
    }
    if (right_child->num_members > 0){
        p->right_child = right_child;
        *tree_size += 1;
    }
    else{
        p->right_child = NULL;
    }

    left_child->members[0] = p->members[0];
    left_child->members[1] = p->members[0] + left_indicies.size() - 1;
    right_child->members[0] = p->members[0] + left_indicies.size();
    right_child->members[1] = p->members[1];

#if TESTFLAG
    cout << " Split panels sucessfully, found " << left_child->num_members << " left particles and " << right_child->num_members << " right particles" << endl;

    cout << "Left panel has particles" << endl;
    for(size_t i=0;i<left_indicies.size();i++){
        cout << left_indicies[i] << endl;
    }

    cout << endl << "Right panel has particles" << endl;
    for(size_t i=0;i<right_indicies.size();i++){
        cout << right_indicies[i] << endl;
    }

#endif



    for(size_t i=0;i<left_indicies.size();i++){
        indicies[i + p->members[0]] = left_indicies[i];
    }

    for(size_t i=0;i<right_indicies.size();i++){
        indicies[i + p->members[0] + left_indicies.size()] = right_indicies[i];
    }

    if(left_child->num_members > N0){
        split_panel(left_child, source_particles, tree_size, leaf_size, indicies);
    }
    else if (left_child->num_members != 0){
        *leaf_size += 1;
        left_child->right_child = NULL;
        left_child->left_child = NULL;
    }

    if(right_child->num_members > N0){
        split_panel(right_child, source_particles, tree_size, leaf_size, indicies);
    }
    else if (right_child->num_members != 0){
        *leaf_size += 1;
        right_child->right_child = NULL;
        right_child->left_child = NULL;
    }

}

/* free_tree_list
 *
 * This function will recursivley free a panel and all its children.
 * Passing the root to this function will free the entire tree.
 *
*/
void free_tree_list(panel *panel){
    if(panel->left_child){free_tree_list(panel->left_child);}
    if(panel->right_child){free_tree_list(panel->right_child);}
    free(panel);
}


void E_MQ_Treecode::init_tree_list(panel *p, panel *tree_list, int *current_id, int *leaf_indicies, int *leaf_id){
    p->id = *current_id;
    tree_list[*current_id] = *p;
    *current_id += 1;

    
    if (p->left_child){
        init_tree_list(p->left_child, tree_list, current_id, leaf_indicies, leaf_id);
    }
    if (p->right_child){
        init_tree_list(p->right_child, tree_list, current_id, leaf_indicies, leaf_id);
    }
   // Handle leafs
    if (!(p->left_child) && !(p->right_child)){
#if TESTFLAG
        cout << "Setting leaf index to " << *current_id-1 << endl;
#endif
        leaf_indicies[*leaf_id] = *current_id - 1;
        *leaf_id += 1;
    }
}

void E_MQ_Treecode::free_tree_list(panel *panel){
    if(panel->left_child){free_tree_list(panel->left_child);}
    if(panel->right_child){free_tree_list(panel->right_child);}
    free(panel);
}

void E_MQ_Treecode::__global__ void init_modified_weights(panel *d_tree_list, double *d_particles, double *d_weights, int source_size, int tree_size){
    const int tree_idx = blockIdx.x;
    const int cheb_idx = threadIdx.x;
    double sum = 0.0; // Memory only cleared after warp is finished, so unused threads should still have this set (?)
    if (cheb_idx >= PP || tree_idx >= tree_size){return;}

    panel *p = d_tree_list + tree_idx;

    double w1;
    if (cheb_idx == 0 || cheb_idx == (PP-1)) {
        w1 = 0.5;
    }
    else {
        w1 = 1.0;
    }
    //if (cheb_idx % 2 == 1) {
    if ( cheb_idx & 0x1 == 1){
        w1 *= -1.0;
    }
    
    double a1; // the clartiy term in paper 
    double modified_weight = 0.0;
    double y;
    const double cheb_pt = p->s[cheb_idx];
    int flag;

    // set up modified weights 
    // Could unroll this loop and have the second half of a warp do half of these iterations if we are fine with limiting the interpolation degree to be <=14
    // Could also just have a switch to check that (and a few versions of this function with various levels of loop unrolling)
    for (int k = p->members[0]; k <= p->members[1]; k++) {
        y = d_particles[k]; // particles in cluster
        flag = -1;
        sum = 0.0;

        if(fabs(y - cheb_pt) <= DBL_MIN) {
            flag = cheb_idx;
        }
        else{
            a1 = w1 /(y - cheb_pt);
            sum += a1;
        }

        // This possibly avoids some thread divergence at the cost of using an extra register
        // and some extra arithmetic ops
        // Appears to produce the same result as above, but this hasn't been tested thouroughly
       //close = fabs(y - p->s[i]) <= DBL_MIN;
       //flag = i*close - (1-close); // = i if close=1, =-1 if close=0
       //a1[i] = w1[i] / (y - p->s[i] + close) * (1-close); // =a1[i] if close=0, =0 if close=1
       //sum += a1[i];
        

        if (flag > -1) {
            // shfl_xor_sync might be better here
            flag = __shfl_sync(FULL_MASK, flag, threadIdx.x);
            a1 = 1.0;
        }
        if (flag > -1 && fabs(y - cheb_pt) > DBL_MIN){
            a1 = 0.0;
        }


        // Use a shfl_down_sync here to reduce sum
        for (int offset = 16; offset>0; offset /= 2){
            sum += __shfl_down_sync(FULL_MASK, sum, offset);
        }
        sum = __shfl_sync(FULL_MASK, sum, 0);

        modified_weight += a1 * d_weights[k] / sum;

    }
    p->modified_weights[threadIdx.x] = modified_weight;
}


void E_MQ_Treecode::init_interaction_lists(panel* leaf, panel* source_panel, int *near_ids, int *far_ids, int *near_index, int *far_index, int leaf_id, int leaf_size, double period, int *total_nears){
    // Check if source panel is a leaf
    if (!source_panel->left_child && !source_panel->right_child){
        near_ids[leaf_id*leaf_size + *near_index] = source_panel->id;
        *near_index += 1;
        *total_nears += 1;
        leaf->near_size += 1;
    }
    else{
        double leaf_radius = leaf->xinterval[1] - leaf->xc;
        double source_radius = source_panel->xinterval[1] - source_panel->xc;
        double distance = std::fabs(leaf->xc - source_panel->xc);
        ///////////////////////////// THIS LINE FOR PERIODIC CONDITIONS //////////////////////////
        distance = std::fmin(distance, std::fabs(period-distance));
        //////////////////////////////////////////////////////////////////////////////////////////
        
#if TESTFLAG
        cout << "Leaf Id, Source ID, Leaf Interval, Source Interval, Leaf radius, source radius, distance, ratio: " << endl;
        cout << leaf->id << "\t" << source_panel->id << "\t [";
        cout << leaf->xinterval[0] << "," << leaf->xinterval[1] << "]\t [";
        cout <<  source_panel->xinterval[0] << "," << source_panel->xinterval[1] << "]\t";
        cout << leaf_radius << "\t" << source_radius << "\t" << distance << "\t";
        cout << (leaf_radius + source_radius) / distance << endl;
#endif
        if ( (leaf_radius + source_radius) / distance < MAC){
            far_ids[leaf_id*leaf_size + *far_index] = source_panel->id;
            *far_index += 1;
            leaf->far_size += 1;
        }
        // If not a far interaction, recursivley check source_panel's children
        else{
            if(source_panel->left_child){
                init_interaction_lists(leaf, source_panel->left_child, near_ids, far_ids, near_index, far_index, leaf_id, leaf_size, period, total_nears);
            }
            if(source_panel->right_child){
                init_interaction_lists(leaf, source_panel->right_child, near_ids, far_ids, near_index, far_index, leaf_id, leaf_size, period, total_nears);
            }
        }
    }
}



double inline E_MQ_Treecode::__device__ kernel(double x, double y, double3 kernel_params){
    double z = (x - y) * kernel_params.x;
    z -= round(z);
    return z * kernel_params.z * rsqrt( z*z + kernel_params.y ) - z;
}


void E_MQ_Treecode::__global__ computesum(double *e_field, panel *tree_list, int *leaf_indicies, double *target_particles, double *source_particles, double *weights, int *d_near_list, int *d_far_list, int leaf_size){
    int idx = blockIdx.x*blockDim.x + threadIdx.x;

    if(idx >= leaf_size) {return;}

    panel *leaf_panel = tree_list + leaf_indicies[idx];

    int blocksize = 128;
    int gridlen = (leaf_panel->num_members + blocksize - 1) / blocksize;
    computepanelsum<<<gridlen, blocksize>>>(e_field, leaf_panel, tree_list, target_particles, source_particles, weights, d_near_list, d_far_list, idx, leaf_size);

}



void E_MQ_Treecode::__global__ computepanelsum(double *e_field, panel *leaf_panel, panel *tree_list, double *target_particles, double *source_particles, double *weights, int *d_near_list, int *d_far_list, int leaf_id, int leaf_size){
    const int idx = blockIdx.x*blockDim.x + threadIdx.x;

    if(idx >= leaf_panel->num_members){return;}

    const int member_idx = leaf_panel->members[0]+idx;
    const double px = source_particles[member_idx];
    double local_e = 0.0;
    double z;
    double norm_epssq = sqrt(1.0 + 4.0 * epsoverLsq);
    panel far_panel;
    panel near_panel;
    //size_t near_members[2];

    // This is running for a single panel, so the far panels and near panels are the same for every thread
    // So, should put relevant source particles, weights, modified weights, and chebyshev points in shared memory
    // Will still need some global acesses across blocks, maybe

    // Should be able to further parallelize these loops
    // Double loop over leafs (outer) and far panels (inner) should be done outside kernel
    for(size_t k=0;k<leaf_panel->far_size;k++){
            far_panel = tree_list[d_far_list[leaf_size * leaf_id + k]];
            for (size_t j=0;j<PP;j++){
//                local_e += kernel(px, far_panel.s[j]) * far_panel.modified_weights[j];
                z = (px - far_panel.s[j])*Linv;
                z = z - round(z);
                local_e += (0.5 * z * norm_epssq * rsqrt( z*z + epsoverLsq  ) - z) * far_panel.modified_weights[j];
            }
       } 

    // Same comment for this loop
    for(size_t k=0;k<leaf_panel->near_size;k++){
        near_panel = tree_list[d_near_list[leaf_size * leaf_id + k]];
        // Might be slightly faster to not load entire panel into memory
        //near_members[0] = tree_list[d_near_list[leaf_size * leaf_id + k]].members[0];
        //near_members[1] = tree_list[d_near_list[leaf_size * leaf_id + k]].members[1];
        for (size_t j=near_panel.members[0];j<=near_panel.members[1];j++){
            // local_e += kernel(px, source_particles[j]) * weights[j];
            z = (px - source_particles[j])*Linv;
            z = z - round(z);
            local_e += (0.5 * z * norm_epssq * rsqrt( z*z + epsoverLsq  ) - z) * weights[j];
        }
    }
    e_field[member_idx] = local_e;
}



void E_MQ_Treecode::__global__ computepanelsum_far(double *e_field, panel *leaf_panel, panel *tree_list, double *target_particles, double *source_particles, int *d_far_list, int leaf_id, int leaf_size, double3 kernel_params){
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    if (idx >= leaf_panel->num_members * leaf_panel->far_size){return;}

    int member_idx = leaf_panel->members[0] + idx / leaf_panel->far_size;
    int far_idx = idx % leaf_panel->far_size;

    double px = source_particles[member_idx];
    double local_e = 0.0;
    panel far_panel = tree_list[d_far_list[leaf_size * leaf_id + far_idx]];
    for (size_t j=0;j<PP;j++){
        local_e += kernel(px, far_panel.s[j], kernel_params) * far_panel.modified_weights[j];
    }

    // Could make a bit faster using the same __shfl_reduce_sum from init_modified_weights
    atomicAdd(e_field + member_idx, local_e);
}



void E_MQ_Treecode::__global__ computepanelsum_near(double *e_field, double *target_particles, double *source_particles, double *weights, uint4* near_data, int total_nears, double3 kernel_params){
    int interaction_id = blockIdx.z * blockDim.z + threadIdx.z;
    if (interaction_id >= total_nears) {return;}

    int target_id = blockIdx.x * blockDim.x + threadIdx.x;
    int source_id = blockIdx.y * blockDim.y + threadIdx.y;

    // Not sure doing this as shared memory helps at all
    /*__shared__ size_t target_mem_0;
    __shared__ size_t source_mem_0;
    __shared__ size_t target_size;
    __shared__ size_t source_size;

    if ( (threadIdx.x == 0) && (threadIdx.y == 0)){
        target_mem_0 = near_data[4*interaction_id];
        source_mem_0 = near_data[4*interaction_id+1];
        target_size = near_data[4*interaction_id+2];
        source_size = near_data[4*interaction_id+3];
    }

    __syncthreads();
*/
/*
    size_t target_mem_0 = near_data[4*interaction_id];
    size_t source_mem_0 = near_data[4*interaction_id+1];
    size_t target_size = near_data[4*interaction_id+2];
    size_t source_size = near_data[4*interaction_id+3];
    */

    uint4 interaction = near_data[interaction_id];
    size_t target_mem_0 = interaction.x;
    size_t source_mem_0 = interaction.y;
    size_t target_size = interaction.z;
    size_t source_size = interaction.w;

    
    if(target_id >= target_size){return;}
    if(source_id >= source_size){return;}

    double target_x = target_particles[target_mem_0 + target_id];
    double source_x = target_particles[source_mem_0 + source_id];
    double source_weight = weights[source_mem_0 + source_id];

    double local_e = kernel(target_x, source_x, kernel_params) * source_weight;


    // IDK how much partially unrolling this loop helps, just trying to get the % of inactive threads down
/*    if(target_id >= target_size){return;}
    size_t reduced_source_size = source_size / UNROLLNUM;
    size_t spilled_source_size = source_size % UNROLLNUM;
    if(source_id >= reduced_source_size + 1){return;}

    double local_e = 0.0;
    double target_x = target_particles[target_mem_0 + target_id];
    if (source_id < reduced_source_size){
#pragma unroll
        for (size_t i=0;i<UNROLLNUM;i++){
            double source_x = source_particles[source_mem_0 + UNROLLNUM*source_id + i];
            double source_weights = weights[source_mem_0 + UNROLLNUM*source_id + i];
            local_e += kernel(target_x, source_x, kernel_params)*source_weights;
        }
    }
    else{
        for (size_t i=0;i<spilled_source_size;i++){
            double source_x = source_particles[source_mem_0 + UNROLLNUM*source_id + i];
            double source_weight = weights[source_mem_0 + UNROLLNUM*source_id + i];
            local_e += kernel(target_x, source_x, kernel_params) * source_weight;
        }
    }
*/
    atomicAdd(e_field + target_mem_0 + target_id, local_e);
}


void E_MQ_Treecode::__global__ re_order(double* e_field, double* ordered_efield, size_t* indicies, int source_size){
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    if (idx >= source_size){return;}

    e_field[indicies[idx]] = ordered_efield[idx];
}

int E_MQ_Treecode::checkcudaerr(cudaError_t err){
    if (err != cudaSuccess){
        cout << "Cuda Error " << err << ": " << cudaGetErrorString(err) << endl;
        return 1;
    }
    else{
        return 0;
    }
}




