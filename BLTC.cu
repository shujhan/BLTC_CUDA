# TEST IF IN FAST

#include "BLTC.hpp"
#include "quicksort.h"
#include<cmath>

#include <iostream>
#include <iomanip>
#include<assert.h>
using std::cout; 
using std::endl;
using namespace std;

const double L = 4*pi;

#define TESTFLAG 0 

#define cdpErrchk(ans) { cdpAssert((ans), __FILE__, __LINE__); }
__device__ void cdpAssert(cudaError_t code, const char *file, int line, bool abort=true)
{
    if (code != cudaSuccess)
    {
        printf("GPU kernel assert: %s %s %d\n", cudaGetErrorString(code), file, line);
	if (abort) assert (0);
    }
}

/* split panel
 *
 * Creates the left and right children of the passed panel p.  If these children are 
 * too large and must themselves be split again, the function is called recursivley.
 * The entire tree can be constructed by passing in the root panel.  In this case
 * tree_size will contain the total number of panels in the tree (minus the root) and
 * leaf_size will contain the number of leaves.
 *
 * This function assumes that xinterval, xc, level, members, and num_members are set in p.
 * It will set these values, set the parent of both children to be p and left_child and
 * right_child of p, initilize near_ids and far_ids to -1, and set the Chebyshev points.  
 * The id and modified_weights attributes are not set.
 *
 */
void split_panel(panel *p, double* source_particles, int *tree_size, int *leaf_size){

    panel *left_child = new panel();
    panel *right_child = new panel();

    left_child->xinterval[0] = p->xinterval[0];
    left_child->xinterval[1] = p->xc;
    left_child->xc = 0.5 * (p->xinterval[0] + p->xc);
    left_child->level = p->level + 1;
    left_child->parent = p;
    for (int k=0;k<PP;k++){
        left_child->s[k] = 0.5 * ( left_child->xinterval[0] + left_child->xinterval[1] + std::cos(k*pi/P)*( left_child->xinterval[1] - left_child->xinterval[0]  ));
    }

    right_child->xinterval[0] = p->xc;
    right_child->xinterval[1] = p->xinterval[1];
    right_child->xc = 0.5 * (p->xinterval[1] + p->xc);
    right_child->level = p->level + 1;
    right_child->parent = p;
    for (int k=0;k<PP;k++){
        right_child->s[k] = 0.5 * ( right_child->xinterval[0] + right_child->xinterval[1] + std::cos(k*pi/P)*( right_child->xinterval[1] - right_child->xinterval[0]  ));
    }

#if TESTFLAG
    cout << "Members are " << p->members[0] << " through " << p->members[1] << endl;
    cout << endl;
    cout << "Splitting particles with " << p->num_members << " Members" << endl;
#endif

    // Handle empty panel case
    if (source_particles[p->members[0]] > p->xc){
        left_child->num_members = 0;
        right_child->members[0] = p->members[0];
        right_child->members[1] = p->members[1];
        right_child->num_members = p->num_members;
    }
    else if (source_particles[p->members[1]] <= p->xc){
        right_child->num_members = 0;
        left_child->members[0] = p->members[0];
        left_child->members[1] = p->members[1];
        left_child->num_members = p->num_members;
    }
    else{
        left_child->members[0] = p->members[0];
        right_child->members[1] = p->members[1];

        for(size_t k=p->members[0];p->members[1];k++){
            if(source_particles[k] > p->xc){
                left_child->members[1] = k-1;
                right_child->members[0] = k;
                break;
            }
        }
        left_child->num_members = left_child->members[1] - left_child->members[0] + 1;
        right_child->num_members = right_child->members[1] - right_child->members[0] + 1;
    }

    
#if TESTFLAG
    cout << " Split panels sucessfully, found " << left_child->num_members << " left particles and " << right_child->num_members << " right particles" << endl;
#endif


    if (left_child->num_members > 0){
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

    if( left_child->num_members > N0  ){
        split_panel(p->left_child, source_particles, tree_size, leaf_size);
    }
    else if (left_child->num_members != 0){
        *leaf_size += 1;
        left_child->right_child = NULL;
        left_child->left_child = NULL;
    }
    if (right_child->num_members > N0 ){
        split_panel(p->right_child, source_particles, tree_size, leaf_size);
    }
    else if (right_child->num_members != 0){
        *leaf_size += 1;
        right_child->right_child = NULL;
        right_child->left_child = NULL;
    }

}

void free_tree_list(panel *panel){
    if(panel->left_child){free_tree_list(panel->left_child);}
    if(panel->right_child){free_tree_list(panel->right_child);}
    free(panel);
}

//TODO This kernel has very low occupancy -- try to optimize
//Might be able to do some loop unrolling / dynamic parallelism
//If nothing else, can identify 1 thread per PP loop iteration (should have next to no communication
//in each case)
//Currently has very low warp efficiency (lots of thread divergence)
// - Should launch with PP*tree_size threads, so that each thread gets one panel and one iteration of each PP loop
//   - Should also keep local variables for several vars (w1[i], p->modified_weights[k], etc) to reduce global memory acesses
__global__ void init_modified_weights(panel *d_tree_list, double *d_particles, double *d_weights, int source_size, int tree_size){
    int idx = blockIdx.x*blockDim.x + threadIdx.x;

    if (idx >= tree_size){return;}

    panel *p = d_tree_list + idx;

    double w1[PP]; 
    for (int i = 0; i < PP; i++) {
        if (i == 0 || i == (PP-1)) {
            w1[i] = 0.5;
        }
        else {
            w1[i] = 1.0;
        }
        if (i % 2 == 1) {
            w1[i] *= -1.0;
        }
    }
    
    double a1[PP]; // the clartiy term in paper 
   
    // Initilize modified weights to zero
    for (int k = 0; k<PP; k++){
        p->modified_weights[k] = 0.0;
    }

    //bool close;
    // set up modified weights 
    for (int k = p->members[0]; k <= p->members[1]; k++) {
        double y = d_particles[k]; // particles in cluster
        int flag = -1;
        double sum = 0.0;
        for (int i = 0; i < PP; i++) {
            if(fabs(y - p->s[i]) <= DBL_MIN) {
                flag = i;
            }
            else{
                a1[i] = w1[i] /(y - p->s[i]);
                sum += a1[i];
            }

            // This possibly avoids some thread divergence at the cost of using an extra register
            // and some extra arithmetic ops
            // Appears to produce the same result as above, but this hasn't been tested thouroughly
           //close = fabs(y - p->s[i]) <= DBL_MIN;
           //flag = i*close - (1-close); // = i if close=1, =-1 if close=0
           //a1[i] = w1[i] / (y - p->s[i] + close) * (1-close); // =a1[i] if close=0, =0 if close=1
           //sum += a1[i];
        }

        if (flag > -1) {
            sum = 1.0;
            for (int j = 0; j < PP; j++) {
                a1[j] = 0.0;
            }  
            a1[flag] = 1.0;
        }

        double D = 1.0 / sum;
        for (int i = 0; i < PP; i++) {
            p->modified_weights[i] += a1[i] * D * d_weights[k];
        }
    }
}

void init_interaction_lists(panel* leaf, panel* source_panel, int *near_ids, int *far_ids, int *near_index, int *far_index, int leaf_id, int leaf_size, double period){
    // Check if source panel is a leaf
    if (!source_panel->left_child && !source_panel->right_child){
        near_ids[leaf_id*leaf_size + *near_index] = source_panel->id;
        *near_index += 1;
        leaf->near_size += 1;
    }
    else{
        //double leaf_radius = leaf->xc - leaf->xinterval[0];
        double leaf_radius = leaf->xinterval[1] - leaf->xc;
        //double source_radius = source_panel->xc - source_panel->xinterval[0];
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
        else{
            if(source_panel->left_child){
                init_interaction_lists(leaf, source_panel->left_child, near_ids, far_ids, near_index, far_index, leaf_id, leaf_size, period);
            }
            if(source_panel->right_child){
                init_interaction_lists(leaf, source_panel->right_child, near_ids, far_ids, near_index, far_index, leaf_id, leaf_size, period);
            }
        }
    }
}


void init_tree_list(panel *p, panel *tree_list, int *current_id, int *leaf_indicies, int *leaf_id){
    
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

// This will eventually read L from the interface class
__device__ double kernel(double x, double y){
    const double eps = 1e-1;
    double z = (x - y)/L;
    z = z - round(z);
    return 0.5 * z * sqrt(1.0 + 4.0 * eps * eps / (L*L)) * rsqrt( z*z + eps*eps/(L*L)  ) - z;
    //return x*y;
}

__global__ void computepanelsum(double *e_field, panel *leaf_panel, panel *tree_list, double *target_particles, double *source_particles, double *weights, int *d_near_list, int *d_far_list, int leaf_id, int leaf_size){
    int idx = blockIdx.x*blockDim.x + threadIdx.x;

    if(idx >= leaf_panel->num_members){return;}

    int member_idx = leaf_panel->members[0]+idx;
    double px = source_particles[member_idx];
    double local_e = 0.0;

    for(size_t k=0;k<leaf_panel->far_size;k++){
            panel far_panel = tree_list[d_far_list[leaf_size * leaf_id + k]];
            for (size_t j=0;j<PP;j++){
                local_e += kernel(px, far_panel.s[j]) * far_panel.modified_weights[j];
            }
       } 

    for(size_t k=0;k<leaf_panel->near_size;k++){
        panel near_panel = tree_list[d_near_list[leaf_size * leaf_id + k]];
        for (size_t j=near_panel.members[0];j<=near_panel.members[1];j++){
            local_e += kernel(px, source_particles[j]) * weights[j];
        }
    }
    e_field[member_idx] = local_e;
}

__global__ void computesum(double *e_field, panel *tree_list, int *leaf_indicies, double *target_particles, double *source_particles, double *weights, int *d_near_list, int *d_far_list, int leaf_size){
    int idx = blockIdx.x*blockDim.x + threadIdx.x;

    if(idx >= leaf_size) {return;}

    panel *leaf_panel = tree_list + leaf_indicies[idx];

    int blocksize = 128;
    int gridlen = (leaf_panel->num_members + blocksize - 1) / blocksize;
    computepanelsum<<<gridlen, blocksize>>>(e_field, leaf_panel, tree_list, target_particles, source_particles, weights, d_near_list, d_far_list, idx, leaf_size);

}

int checkcudaerr(cudaError_t err){
    if (err != cudaSuccess){
        cout << "Cuda Error " << err << ": " << cudaGetErrorString(err) << endl;
        return 1;
    }
    else{
        return 0;
    }
}

// TODO
// Multi-GPU
// non-unity weights (just need to re-order properly)
// target particles differing from source particles (probably just need to swap source to target in a few places)

void BLTC(double *e_field, double *source_particles, double *target_particles, double *weights, 
        size_t e_field_size, size_t source_size, size_t target_size){

#if TESTFLAG
    cout << "Input particles:" << endl;
    for(size_t k=0;k<source_size;k++){
        cout << "x[" << k << "] = " << source_particles[k] << endl;
    }
#endif
    // TODO Need to sort particles and also re-sort the corresponding weights
    size_t* source_indicies = (size_t*)malloc(sizeof(size_t)*source_size);
    for(size_t k=0;k<source_size;k++){
        source_indicies[k] = k;
    }
    quicksort(source_particles, (int)source_size, source_indicies);
#if TESTFLAG
    cout << endl << "Sorted particles:" << endl << endl;;
    for(size_t k=0;k<source_size;k++){
        cout << "x[" << k << "] = " << source_particles[k] << endl;
    }
#endif



    panel root;
    double xmin = source_particles[0];
    double xmax = source_particles[0];
    root.members[0] = 0;
    root.members[1] = source_size-1;
    for (size_t k=0; k<source_size; k++){
        if(source_particles[k] < xmin){ xmin = source_particles[k]; }
        if(source_particles[k] > xmax){ xmax = source_particles[k]; }
    }
    //root.xinterval[0] = xmin - 0.001;
    //root.xinterval[1] = xmax + 0.001;
    //root.xc = (xmax + xmin) / 2;
    //root.xinterval[0] = 0.0;
    //root.xinterval[1] = L;
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


    if (source_size > N0){
        split_panel(&root, source_particles, &tree_size, &leaf_size);
    }
    else{ 
        root.left_child = NULL;
        root.right_child = NULL;
        leaf_size=1;  
    }

    cout << "Set up root panel" << endl;

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

    for(int k=0;k<leaf_size;k++){
        int near_index = 0;
        int far_index = 0;
        init_interaction_lists(tree_list+leaf_indicies[k], &root, near_interactions, far_interactions, &near_index, &far_index, k, leaf_size, L); 
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

    // TODO Move initilization/memory transfer of everything except tree items to be before the tree is
    // constructed (concurrency between data tranfer to GPU and CPU operations)

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

    errcode = cudaMalloc(&d_tree_list, tree_size*sizeof(panel));
    if(checkcudaerr(errcode) != 0){cout << "Failed allocating tree list" << endl;}

    errcode = cudaMalloc(&d_particles, source_size*sizeof(double));
    if(checkcudaerr(errcode) !=0){cout << "Failed allocating source particles" << endl;}

    errcode = cudaMalloc(&d_targets, target_size*sizeof(double));
    if(checkcudaerr(errcode) !=0){cout << "Failed allocating target particles" << endl;}

    errcode = cudaMalloc(&d_weights, source_size*sizeof(double));
    if(checkcudaerr(errcode) != 0){cout << "Failed allocating source weights" << endl;}

    errcode = cudaMalloc(&d_efield, target_size*sizeof(double));
    if(checkcudaerr(errcode) != 0){cout << "Failed allocating e_field" << endl;}
    
    errcode = cudaMalloc(&d_near_list, leaf_size*leaf_size*sizeof(int));
    if(checkcudaerr(errcode) != 0){cout << "Failed allocating near interaction list" << endl;}
    
    errcode = cudaMalloc(&d_far_list, leaf_size*leaf_size*sizeof(int));
    if(checkcudaerr(errcode) != 0){cout << "Failed allocating far interaction list" << endl;}

    errcode = cudaMalloc(&d_leaf_indicies, leaf_size*sizeof(int));
    if(checkcudaerr(errcode) != 0){cout << "Failed allocating leaf indicies" << endl;}

    errcode = cudaMemcpy(d_tree_list, tree_list, tree_size*sizeof(panel), cudaMemcpyHostToDevice);
    if(checkcudaerr(errcode) != 0){cout << "Failed copying tree list to device" << endl;}
 
    errcode = cudaMemcpy(d_particles, source_particles, source_size*sizeof(double), cudaMemcpyHostToDevice);
    if(checkcudaerr(errcode) != 0){cout << "Failed copying source particles to device" << endl;}

    errcode = cudaMemcpy(d_targets, target_particles, target_size*sizeof(double), cudaMemcpyHostToDevice);
    if(checkcudaerr(errcode) != 0){cout << "Failed copying target particles to device" << endl;}

    errcode = cudaMemcpy(d_weights, weights, source_size*sizeof(double), cudaMemcpyHostToDevice);
    if(checkcudaerr(errcode) != 0){cout << "Failed copying source weights to device" << endl;}

    errcode = cudaMemcpy(d_near_list, near_interactions, leaf_size*leaf_size*sizeof(int), cudaMemcpyHostToDevice);
    if(checkcudaerr(errcode) != 0){cout << "Failed copying near interactions to device" << endl;}

    errcode = cudaMemcpy(d_far_list, far_interactions, leaf_size*leaf_size*sizeof(int), cudaMemcpyHostToDevice);
    if(checkcudaerr(errcode) != 0){cout << "Failed copying far interactions to device" << endl;}

    errcode = cudaMemcpy(d_leaf_indicies, leaf_indicies, leaf_size*sizeof(int), cudaMemcpyHostToDevice);
    if(checkcudaerr(errcode) != 0){cout << "Failed copying leaf indicies to device" << endl;}

    int blocksize = 128;
    // TODO Should this be tree_size + blocksize + 1 ?
    int gridlen = (source_size + blocksize - 1) / blocksize;
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

    if (root.left_child){
        free_tree_list(root.left_child);
    }
    if (root.right_child){
        free_tree_list(root.right_child);
    }

    //Should set this dynamically eventually
    cudaDeviceSetLimit(cudaLimitDevRuntimePendingLaunchCount, 16384);

    gridlen = (leaf_size + blocksize - 1) / blocksize;
    computesum<<<gridlen,blocksize>>>(d_efield, d_tree_list, d_leaf_indicies, d_targets, d_particles, d_weights, d_near_list, d_far_list, leaf_size);

    cout << "Computed BLTC sum" << endl;

    double* ordered_e_field = (double*)malloc(sizeof(double)*source_size);
    errcode = cudaMemcpy(ordered_e_field, d_efield, target_size*sizeof(double), cudaMemcpyDeviceToHost);
    if(checkcudaerr(errcode) != 0){cout << "Failed to copy e_field to host" << endl;}

    // TODO This can be done more quickly with the GPU
    for(size_t k=0;k<source_size;k++){
        e_field[source_indicies[k]] = ordered_e_field[k];
    }
    free(source_indicies);

#if TESTFLAG
    for (size_t k=0;k<source_size;k++){
        cout << "e[" << k << "] = " << e_field[k] << endl;
    }
#endif

    cudaFree(d_tree_list);
    cudaFree(d_particles);
    cudaFree(d_weights);
    cudaFree(d_leaf_indicies);
    cudaFree(d_efield);

}
