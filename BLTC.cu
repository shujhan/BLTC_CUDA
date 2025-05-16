#include "BLTC.hpp"
#include<cmath>

// TODO 5/14 Notes
/*
 * - The bug we thought was in far interactions probably was not.  The issue comes into how
 *   we were indexing the interactions when creating near_data.  We used leaf_indicies[interaction_id]
 *   as the source panel id, but interaction_id is already the source panel id.  This happened to look like
 *   it worked in the direct sum case because everything interacts with everything anyways, so indexing them
 *   in a schizophernic way accidentally worked
 *
 *   - Still not completley correct, but much closer
 *
 * - A major speedup can be obtained using the symmetry of near interactions -- if panel A has a near interaction
 *   with panel B, then panel B has a near interaction with panel A.  So we should be able to roughly half the computation needed.
 *   
 * - A small catch is associated with the above; as written in code, near interactions are not always symmetric (look at 128x128 case,
 *   panel 10 has a near interaction with panel 6, but panel 6 does not have a near interaction with panel 10).  This is the same
 *   case as before where we get a seperation (almost) exactly equal to the MAC (I think).  We could try to build this symmetry
 *   directly into init_interaction_lists, maybe.
 *
 */


#include <iostream>
#include <iomanip>
#include<assert.h>
#include<vector>
using std::cout; 
using std::endl;
using namespace std;

#define UNROLLNUM 1

/*
const double L = 4*pi;
const double Linv = 1.0/L; // multiplying by Linv in kernel is much faster than dividing by L (~30% for direct sum on 128x128)
const double eps = 1e-1;
const double epsoverLsq = eps*eps*Linv*Linv;
*/

#define TESTFLAG 0 

#define FULL_MASK 0xffffffff

#define cdpErrchk(ans) { cdpAssert((ans), __FILE__, __LINE__); }
__device__ void cdpAssert(cudaError_t code, const char *file, int line, bool abort=true)
{
    if (code != cudaSuccess)
    {
        printf("GPU kernel assert: %s %s %d\n", cudaGetErrorString(code), file, line);
	if (abort) assert (0);
    }
}

//TODO This docstring is out of date - does not include the new sorting methdo
/* split panel
 *
 * Creates the left and right children of the passed panel p.  If these children are 
 * too large and must themselves be split again, the function is called recursivley.
 * The entire tree can be constructed by passing in the root panel.  In this case
 * tree_size will contain the total number of panels in the tree (minus the root) and
 * leaf_size will contain the number of leaves.  Empty panels are excluded.
 *
 * This function assumes that xinterval, xc, level, members, and num_members are set in p.
 * It will set these values, set the parent of both children to be p and left_child and
 * right_child of p, initilize near_ids and far_ids to -1, and set the Chebyshev points.  
 * The id and modified_weights attributes are not set.
 *
 * When this function is called on the root, almost always tree_size and leaf_size should
 * be set to 0.
 */
void split_panel(panel *p, double *source_particles, int *tree_size, int *leaf_size, size_t *indicies){
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


/* init_modified_weights
 *
 * Calculates the modified weights for the BLTC algorithm.
 *
 * This kernel is intended to be called with block sizes of one warp, and assumes that the interpolation degree
 * is less than one minus the size of one warp (31).  The expectation is that every Chebyshev point in every panel has
 * one thread associated with it.  Each warp will calculate the modified weights for one panel.
 *
 */
__global__ void init_modified_weights(panel *d_tree_list, double *d_particles, double *d_weights, int source_size, int tree_size){
//    unsigned mask = __ballot_sync(FULL_MASK, threadIdx.x < PP);
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

/* init_interaction_lists
 *
 * This function will recursivley initilize the interaction list for the given leaf.  
 * Near interactions will be inserted to the array near_ids and far interactions
 * will be inserted into far_ids.  Each of these arrays should be of length leaf_size*leaf_size.
 * Each sequential leaf_size subset of these id arrays will contain the ids of panels which the source
 * panel has a near/far interaction with, followed by -1's to pad to the maximum length of leaf_size.
 *
 * The intention is for this function to be called initially with the source_panel set to the root panel,
 * and near_index and far_index both set to 0.
 *
 * TODO: This could probably be wrapped into a class interface, and an overloaded method could be defined without
 *       the source_panel, near_index, and far_index which simply calls this method with source_panel=root, near_index
 *       and far_index set to 0.
 *
 */
void init_interaction_lists(panel* leaf, panel* source_panel, int *near_ids, int *far_ids, int *near_index, int *far_index, int leaf_id, int leaf_size, double period, int *total_nears){
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

/* init_tree_list
 *
 * Recursivley initilizes the tree list and leaf list.  
 * This is just an array which contains all panels in the tree, keeping track of which
 * indicies in this array correspond to leaves.  Essentially flattens the tree for easier
 * use on the GPU.
 *
 * This function is intended to be called with the panel set to the root, and current_id and leaf_id
 * set to 0.  leaf_indicies should have length leaf_size, and tree_list should have length tree_size.
 *
 */
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
__device__ inline double kernel(double x, double y, double3 kernel_params){
/*    double z = (x - y)*Linv;
    z -= round(z);
    return 0.5 * z * sqrt(1.0 + 4.0 * epsoverLsq) * rsqrt( z*z + epsoverLsq ) - z;
    */
    double z = (x - y) * kernel_params.x;
    z -= round(z);
    return z * kernel_params.z * rsqrt( z*z + kernel_params.y ) - z;
}

__global__ void computepanelsum_far(double *e_field, panel *leaf_panel, panel *tree_list, double *target_particles, double *source_particles, int *d_far_list, int leaf_id, int leaf_size, double3 kernel_params){

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

// TODO This might be a bit faster if we use a 2D grid to get memberidx and near_idx instead of flat indexing
// Can we split this kernel up even further?  Paralleize the inner loop as well?  Not sure that would be better.
// Might entail tracking the number of particles we have direct interactions with when constructing the interaction list
// (shouldn't be a big deal)
/*
__global__ void computepanelsum_near(double *e_field, panel *leaf_panel, panel *tree_list, double *target_particles, double *source_particles, double *weights, int *d_near_list, int leaf_id, int leaf_size){

    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    if (idx >= leaf_panel->num_members * leaf_panel->near_size){return;}

    int member_idx = leaf_panel->members[0] + idx / leaf_panel->near_size;
    int near_idx = idx % leaf_panel->near_size;

    double px = source_particles[member_idx];
    double local_e = 0.0;
    panel near_panel = tree_list[d_near_list[leaf_size * leaf_id + near_idx]];


    for (size_t j=near_panel.members[0];j<=near_panel.members[1];j++){
        // TODO try inlining kernel
        local_e += kernel(px, source_particles[j]) * weights[j];
    }

    // This can be sped up, either using shared memory or writing directly to a temporary array that we can
    // later reduce over
    atomicAdd(e_field + member_idx, local_e);
}
*/

// TODO might be faster to have each thread handle a couple interactions instead of just one for memory reasons
/*
__global__ void computepanelsum_near(double *e_field, panel *leaf_panel, size_t left_mem, size_t right_mem, double *target_particles, double *source_particles, double *weights, int leaf_size){

    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    int num_interactions = right_mem - left_mem + 1;
    if (idx >= leaf_panel->num_members * num_interactions){return;}

    int member_idx = leaf_panel->members[0] + idx / num_interactions;
    int near_idx = idx % num_interactions;

    double px = source_particles[member_idx];
    double local_e = kernel(px, source_particles[left_mem + near_idx]) * weights[left_mem + near_idx];

    atomicAdd(e_field + member_idx, local_e);

}
*/

/* 5/13/25
__global__ void computepanelsum_near(double *e_field, panel* tree_list, int* near_interactions, double *target_particles, double *source_particles, double *weights, int leaf_size, int* leaf_indicies){

    int idx = blockIdx.x * blockDim.x + threadIdx.x;

    // Might be better to swap these for memory acess patterns
    int panel_flat_id = idx % (leaf_size*leaf_size);
    int leaf_particle_id = idx / (leaf_size*leaf_size);


    int leaf_id = panel_flat_id % leaf_size;
    int source_id = panel_flat_id / leaf_size;

    int source_panel_id = near_interactions[leaf_id * leaf_size + source_id];
    if (source_panel_id == -1){return;}

    panel leaf = tree_list[leaf_indicies[leaf_id]];
    panel source = tree_list[source_panel_id];

    if(leaf_particle_id >= leaf.num_members){return;}

    double target_particle = target_particles[leaf.members[0] + leaf_particle_id];
    double source_particle;
    double local_e = 0.0;

    for(size_t i=0;i<source.num_members;i++){
        source_particle = source_particles[source.members[0] + i];
        local_e += kernel(target_particle, source_particle) * weights[source.members[0] + i];
    }

    atomicAdd(e_field + leaf.members[0] + leaf_particle_id, local_e);
    
}
*/

// Rather than indexing with one index for a leaf and one index for a source panel, have a single index representing
// a leaf-source pair.  We can count the total number of pairs needed, and effectivley give each an id from which the
// leaf and source seperatley can be recovered.
// Implementation-wise, it might be easiest to create a tuple structure and then just make an array of tuples 
// of (leaf_id, near_id) pairs which can be indexed into
// We could even use 3-tuples to specify (leaf_id, near_id, number of interactions) or 4-tuples for
// (leaf_id, near_id, num particles in leaf, num particles in source) to get an easier count of how many threads are needed
// One thread dimension indexes the tuple array (and the leaf-source pair by proxy), the other dimension is a flattened index
// for the leaf particle index and source particle index

__global__ void computepanelsum_near(double *e_field, double *target_particles, double *source_particles, double *weights, uint4* near_data, int total_nears, double3 kernel_params){

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

/*
__global__ void computepanelsum_near(double *e_field, double *target_particles, double *source_particles, double *weights, size_t* near_data, int total_nears){

    //int particle_flat_id = blockIdx.x * blockDim.x + threadIdx.x;
    int interaction_id = blockIdx.z * blockDim.z + threadIdx.z;
    if (interaction_id >= total_nears) {return;}

    //tuple4 interaction = near_data[interaction_id];
        //int target_id = particle_flat_id % N0;
    //int source_id = particle_flat_id / N0;
    int target_id = blockIdx.x * blockDim.x + threadIdx.x;
    int source_id = blockIdx.y * blockDim.y + threadIdx.y;


    size_t target_mem_0 = near_data[4*interaction_id];
    size_t source_mem_0 = near_data[4*interaction_id+1];
    size_t target_size = near_data[4*interaction_id+2];
    size_t source_size = near_data[4*interaction_id+3];
    if(target_id >= target_size){return;}
    double target_x = target_particles[target_mem_0 + target_id];

    
    double local_e = 0.0;
    double source_x;
    double source_weight;
    for(size_t i=0;i<source_size;i++){
        source_x = source_particles[source_mem_0 + i];
        source_weight = weights[source_mem_0 + i];
        local_e += kernel(target_x, source_x) * source_weight;
    }

    atomicAdd(e_field + target_mem_0 + target_id, local_e);
        
}
*/

/* 5/14/25
__global__ void computepanelsum_near(double *e_field, panel* tree_list, int* near_interactions, double *target_particles, double *source_particles, double *weights, int leaf_size, int* leaf_indicies, tuple6* near_data){

    int idx = blockIdx.x * blockDim.x + threadIdx.x;

    int particle_flat_id = idx % N02;
    int interaction_id = idx / N02;

    tuple6 interaction = near_data[interaction_id];
    int target_id = particle_flat_id % N0;
    int source_id = particle_flat_id / N0;

    if(target_id >= interaction.target_size){return;}
    if(source_id >= interaction.source_size){return;}

    double target_x = target_particles[interaction.target_mem_0 + target_id];
    double source_x = source_particles[interaction.source_mem_0 + source_id];

    double local_e = kernel(target_x, source_x) * weights[interaction.source_mem_0 + source_id];

    atomicAdd(e_field + interaction.target_mem_0 + target_id, local_e);
        
}
*/



/* 5/13/25
__global__ void computepanelsum_near(double *e_field, panel* tree_list, int* near_interactions, double *target_particles, double *source_particles, double *weights, int leaf_size, int* leaf_indicies){

    int idx = blockIdx.x * blockDim.x + threadIdx.x;

    // Might be better to swap these for memory acess patterns
    int panel_flat_id = idx % (leaf_size*leaf_size);
    int interaction_flat_id = idx / (leaf_size*leaf_size);

    int leaf_id = panel_flat_id % leaf_size;
    int source_id = panel_flat_id / leaf_size;

    int source_panel_id = near_interactions[leaf_id * leaf_size + source_id];
    if (source_panel_id == -1){return;}

    int leaf_particle_id = interaction_flat_id % (N0*N0);
    int source_particle_id = interaction_flat_id / (N0*N0);

    panel leaf = tree_list[leaf_indicies[leaf_id]];
    panel source = tree_list[source_panel_id];

    if(leaf_particle_id >= leaf.num_members){return;}
    if(source_particle_id >= source.num_members){return;}

    double target_particle = target_particles[leaf.members[0] + leaf_particle_id];
    double source_particle = source_particles[source.members[0] + source_particle_id];
    double local_e = kernel(target_particle, source_particle) * weights[source.members[0] + source_particle_id];

    atomicAdd(e_field + leaf.members[0] + leaf_particle_id, local_e);
    
}
*/
/*
__global__ void computepanelsum_near(double *e_field, panel* tree_list, int* near_interactions, double *target_particles, double *source_particles, double *weights, int leaf_size, int* leaf_indicies){

    int leaf_id = blockIdx.x * blockDim.x + threadIdx.x;
    int target_id = blockIdx.y * blockDim.y + threadIdx.y;

    if (leaf_id >= leaf_size) {return;}

    double local_e = 0.0;

    panel leaf_panel = tree_list[leaf_indicies[leaf_id]];

    if (target_id >= leaf_panel.num_members) {return;}

    double target_particle = target_particles[leaf_panel.members[0] + target_id];

    for (size_t i=0;i<leaf_panel.near_size;i++){
        panel near_panel = tree_list[near_interactions[leaf_id*leaf_size + i]];
        for (size_t j=0;j<near_panel.num_members;j++){
            double source_particle = source_particles[near_panel.members[0] + j];
            
            local_e += kernel(target_particle, source_particle) * weights[near_panel.members[0] + j];
        }
    }

    e_field[leaf_panel.members[0] + target_id] = local_e;
    //atomicAdd(e_field + leaf_panel.members[0] + target_id, local_e);

}
*/
/*
__global__ void computepanelsum_near(double *e_field, panel *leaf_panel, size_t left_mem, size_t right_mem, double *target_particles, double *source_particles, double *weights, int leaf_size){

    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    int num_interactions = (right_mem - left_mem + 1) / 16;
    if (idx >= leaf_panel->num_members * num_interactions){return;}

    int member_idx = leaf_panel->members[0] + idx / num_interactions;
    int near_idx = 16*(idx % num_interactions);

    double px = source_particles[member_idx];
    double local_e = kernel(px, source_particles[left_mem + near_idx]) * weights[left_mem + near_idx];
    local_e += kernel(px, source_particles[left_mem + near_idx + 1]) * weights[left_mem + near_idx + 1];
    local_e += kernel(px, source_particles[left_mem + near_idx + 2]) * weights[left_mem + near_idx + 2];
    local_e += kernel(px, source_particles[left_mem + near_idx + 3]) * weights[left_mem + near_idx + 3];
    local_e += kernel(px, source_particles[left_mem + near_idx + 4]) * weights[left_mem + near_idx + 4];
    local_e += kernel(px, source_particles[left_mem + near_idx + 5]) * weights[left_mem + near_idx + 5];
    local_e += kernel(px, source_particles[left_mem + near_idx + 6]) * weights[left_mem + near_idx + 6];
    local_e += kernel(px, source_particles[left_mem + near_idx + 7]) * weights[left_mem + near_idx + 7];
    local_e += kernel(px, source_particles[left_mem + near_idx + 8]) * weights[left_mem + near_idx + 8];
    local_e += kernel(px, source_particles[left_mem + near_idx + 9]) * weights[left_mem + near_idx + 9];
    local_e += kernel(px, source_particles[left_mem + near_idx + 10]) * weights[left_mem + near_idx + 10];
    local_e += kernel(px, source_particles[left_mem + near_idx + 11]) * weights[left_mem + near_idx + 11];
    local_e += kernel(px, source_particles[left_mem + near_idx + 12]) * weights[left_mem + near_idx + 12];
    local_e += kernel(px, source_particles[left_mem + near_idx + 13]) * weights[left_mem + near_idx + 13];
    local_e += kernel(px, source_particles[left_mem + near_idx + 14]) * weights[left_mem + near_idx + 14];
    local_e += kernel(px, source_particles[left_mem + near_idx + 15]) * weights[left_mem + near_idx + 15];

    atomicAdd(e_field + member_idx, local_e);

}
*/


// TODO Might be worth breaking out two seperate kernels, one for near interactions and one for far interactions
__global__ void computepanelsum(double *e_field, panel *leaf_panel, panel *tree_list, double *target_particles, double *source_particles, double *weights, int *d_near_list, int *d_far_list, int leaf_id, int leaf_size){
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

__global__ void computesum(double *e_field, panel *tree_list, int *leaf_indicies, double *target_particles, double *source_particles, double *weights, int *d_near_list, int *d_far_list, int leaf_size){
    int idx = blockIdx.x*blockDim.x + threadIdx.x;

    if(idx >= leaf_size) {return;}

    panel *leaf_panel = tree_list + leaf_indicies[idx];

    int blocksize = 128;
    int gridlen = (leaf_panel->num_members + blocksize - 1) / blocksize;
    computepanelsum<<<gridlen, blocksize>>>(e_field, leaf_panel, tree_list, target_particles, source_particles, weights, d_near_list, d_far_list, idx, leaf_size);

}

__global__ void re_order(double* e_field, double* ordered_efield, size_t* indicies, int source_size){
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    if (idx >= source_size){return;}

    e_field[indicies[idx]] = ordered_efield[idx];
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
// Avoid particle sort -- probably requires changing panel data structure a bit
//    - Could use a hash table to associate panels with particle indicies
//    - Might be able to get away with only storing particle indicies with leafs, then tracking which
//      leafs are children of a particular panel (maybe even a second hash table)
//      Advatange would be that we do not need to keep as much memory per value of the first hash table
//      cuCollections could be useful here, but it is still under pretty active development so might not be super stable
// Should check if anything can be accelerated with CCCL (Thrust/CUB)

void BLTC(double *e_field, double *source_particles, double *target_particles, double *weights, 
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
