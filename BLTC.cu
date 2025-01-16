#include "BLTC.hpp"
#include "quicksort.h"
#include<cmath>

#include <iostream>
#include <iomanip>
#include<assert.h>
using std::cout; 
using std::endl;
using namespace std;

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
    left_child->xc = (p->xinterval[0] + p->xc)/2.0;
    left_child->level = p->level + 1;
    left_child->parent = p;
    left_child->num_members = 0;
    for (int k=0;k<PP;k++){
        left_child->s[k] = 0.5 * ( left_child->xinterval[0] + left_child->xinterval[1] + std::cos(k*pi/PP)*( left_child->xinterval[1] - left_child->xinterval[0]  ));
    }

    right_child->xinterval[0] = p->xc;
    right_child->xinterval[1] = p->xinterval[1];
    right_child->xc = (p->xinterval[1] + p->xc)/2.0;
    right_child->level = p->level + 1;
    right_child->parent = p;
    right_child->num_members = 0;
    for (int k=0;k<PP;k++){
        right_child->s[k] = 0.5 * ( right_child->xinterval[0] + right_child->xinterval[1] + std::cos(k*pi/PP)*( right_child->xinterval[1] - right_child->xinterval[0]  ));
    }

#if TESTFLAG
    cout << "Members are " << p->members[0] << " through " << p->members[1] << endl;
    cout << endl;
    cout << "Splitting particles with " << p->num_members << " Members" << endl;
#endif

    left_child->members[0] = p->members[0];
    right_child->members[1] = p->members[1];
    for(size_t k=p->members[0];p->members[1];k++){
        if(source_particles[k] >= p->xc){
            left_child->members[1] = k-1;
            right_child->members[0] = k;
            break;
        }
    }
    left_child->num_members = left_child->members[1] - left_child->members[0] + 1;
    right_child->num_members = right_child->members[1] - right_child->members[0] + 1;

    
#if TESTFLAG
    cout << " Split panels sucessfully, found " << left_child->num_members << " left particles and " << right_child->num_members << " right particles" << endl;
#endif


    p->left_child = left_child;
    p->right_child = right_child;

    *tree_size += 2;

    if( left_child->num_members > N0  ){
        split_panel(p->left_child, source_particles, tree_size, leaf_size);
    }
    else{
        *leaf_size += 1;
    }
    if (right_child->num_members > N0 ){
        split_panel(p->right_child, source_particles, tree_size, leaf_size);
    }
    else{
        *leaf_size += 1;
    }

}

void free_tree_list(panel *panel){
    if(panel->left_child){free_tree_list(panel->left_child);}
    if(panel->right_child){free_tree_list(panel->right_child);}
    free(panel);
}

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
        double leaf_radius = leaf->xc - leaf->xinterval[0];
        double source_radius = source_panel->xc - source_panel->xinterval[0];
        double distance = std::fabs(leaf->xc - source_panel->xc);
        ///////////////////////////// THIS LINE FOR PERIODIC CONDITIONS //////////////////////////
//        distance = std::fmin(distance, period-distance);
//        /////////////////////////////////////////////////////////////////////////////////////
#if TESTFLAG
        cout << "Leaf Id, Source ID, Leaf radius, source radius, distance, ratio: " << endl;
        cout << leaf->id << "\t" << source_panel->id << "\t" << leaf_radius << "\t" << source_radius << "\t" << distance << "\t" << (leaf_radius + source_radius) / distance << endl;
#endif
        if ( (leaf_radius + source_radius) / distance < MAC){
            far_ids[leaf_id*leaf_size + *far_index] = source_panel->id;
            *far_index += 1;
            leaf->far_size += 1;
        }
        else{
            init_interaction_lists(leaf, source_panel->left_child, near_ids, far_ids, near_index, far_index, leaf_id, leaf_size, period);
            init_interaction_lists(leaf, source_panel->right_child, near_ids, far_ids, near_index, far_index, leaf_id, leaf_size, period);
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
        leaf_indicies[*leaf_id] = *current_id - 1;
        *leaf_id += 1;
    }
}

// This will eventually read L from the interface class
__device__ double kernel(double x, double y){
/*    const double eps = 1e-8;
    double z = x - y;
    z = z - round(z);
    return 0.5 * z * sqrt(1.0 + 4.0 * eps * eps) * rsqrt( z*z + eps*eps  ) - z;
    */
    return x*y;
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
    size_t source_indicies[source_size];
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

    const double L = 1.0;


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
    root.xinterval[0] = 0.0;
    root.xinterval[1] = L;
    root.xc = L/2;
    root.level = 0;
    root.num_members = source_size;
    for (int k=0;k<PP;k++){
        root.s[k] = 0.5 * ( root.xinterval[0] + root.xinterval[1] + std::cos(k*pi/PP)*( root.xinterval[1] - root.xinterval[0]  ));
    }
    
    int tree_size = 1;
    int leaf_size = 0;


    split_panel(&root, source_particles, &tree_size, &leaf_size);

    panel tree_list[tree_size];
    int leaf_indicies[leaf_size];

    int id = 0;
    int leaf_id = 0;


    init_tree_list(&root, tree_list, &id, leaf_indicies, &leaf_id);

#if TESTFLAG
    cout << endl;
    cout << "Tree list is" << endl;
    for (int k=0;k<tree_size;k++){
        cout << "Range: [" << tree_list[k].xinterval[0] << "," << tree_list[k].xinterval[1] << "]\t ID: " << tree_list[k].id << endl;
    }
    cout << endl;
    cout << "Leaf indicies are" << endl;
    for (int k=0;k<leaf_size;k++){
        cout << "Range: [" << tree_list[leaf_indicies[k]].xinterval[0] << "," << tree_list[leaf_indicies[k]].xinterval[1] << "]\t ID: " << leaf_indicies[k] << endl;
    }
    cout << endl;
#endif

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

#if TESTFLAG
    cout << endl;
    for(int k=0;k<leaf_size;k++){
        cout << "Leaf " << k << " has ID " << leaf_indicies[k] << " and near Ids" << endl;
        for(int i=0;i<leaf_size;i++){
            if (near_interactions[k*leaf_size + i] == -1){break;}
            cout << near_interactions[k*leaf_size + i] << "\t";
        }
        cout << endl;
        cout << "and far Ids" << endl;
        for(int i=0;i<leaf_size;i++){
            if (far_interactions[k*leaf_size + i] == -1){break;}
            cout << far_interactions[k*leaf_size + i] << "\t";
        }
        cout << endl;
    }
    cout << endl;
#endif

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
    int gridlen = (source_size + blocksize - 1) / blocksize;
    init_modified_weights<<<gridlen,blocksize>>>(d_tree_list, d_particles, d_weights, source_size, tree_size);

#if TESTFLAG
    errcode = cudaMemcpy(tree_list, d_tree_list, tree_size*sizeof(panel), cudaMemcpyDeviceToHost);
    if(errcode != 0){cout << "Failed to copy tree list to host" << endl;}

    for(int k=0;k<tree_size;k++){
        cout << "Modified weights for panel " << k << endl;
        for (int j=0;j<PP;j++){
            cout << tree_list[k].modified_weights[j] << "\t";
        }
        cout << endl;
        cout << "Chebyshev Points for panel " << k << endl;
        for (int j=0;j<PP;j++){
            cout << tree_list[k].s[j] << "\t";
        }
        cout << endl;
    }
    cout << endl;
#endif

    free_tree_list(root.left_child);
    free_tree_list(root.right_child);

    //Should set this dynamically eventually
    cudaDeviceSetLimit(cudaLimitDevRuntimePendingLaunchCount, 16384);

    gridlen = (leaf_size + blocksize - 1) / blocksize;
    computesum<<<gridlen,blocksize>>>(d_efield, d_tree_list, d_leaf_indicies, d_targets, d_particles, d_weights, d_near_list, d_far_list, leaf_size);

    double ordered_e_field[source_size];
    errcode = cudaMemcpy(ordered_e_field, d_efield, target_size*sizeof(double), cudaMemcpyDeviceToHost);
    if(checkcudaerr(errcode) != 0){cout << "Failed to copy e_field to host" << endl;}

    for(size_t k=0;k<source_size;k++){
        e_field[source_indicies[k]] = ordered_e_field[k];
    }

    cudaFree(d_tree_list);
    cudaFree(d_particles);
    cudaFree(d_weights);
    cudaFree(d_leaf_indicies);
    cudaFree(d_efield);

}
