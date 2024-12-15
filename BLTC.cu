#include "BLTC.hpp"
//#include "kernel.cuh"
#include<cmath>

#include <iostream>
#include <iomanip>
#include<cstring> // for memset
using std::cout, std::endl;
using namespace std;

#define TESTFLAG 1 

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
    for(int k=0;k<Nmax;k++){
        left_child->near_ids[k] = -1;
        left_child->far_ids[k] = -1;
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
    for(int k=0;k<Nmax;k++){
        right_child->near_ids[k] = -1;
        right_child->far_ids[k] = -1;
    }

#if TESTFLAG
    cout << "Members are ";
    for (size_t k=0;k<p->num_members;k++){
        cout << p->members[k] << "\t";
    }
    cout << endl;
    cout << "Splitting particles with " << p->num_members << " Members" << endl;
#endif
    for (size_t k=0;k<p->num_members;k++){
        if (source_particles[p->members[k]] < p->xc){
            left_child->members[left_child->num_members] = p->members[k];
            left_child->num_members += 1;
        }
        else{
            right_child->members[right_child->num_members] = p->members[k];
            right_child->num_members += 1;
        }
    }

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

__global__ void init_modified_weights(panel* d_tree_list, double *d_particles, double *d_weights, int source_size, int tree_size){
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
    for (int k = 0; k < p->num_members; k++) {
        double y = d_particles[p->members[k]]; // particles in cluster
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
            p->modified_weights[i] += a1[i] * D * d_weights[p->members[i]];
        }
    }


}

void init_interaction_lists(panel* leaf, panel* source_panel, int *near_index, int *far_index, double period){
    // Check if source panel is a leaf
    if (!source_panel->left_child){
        leaf->near_ids[*near_index] = source_panel->id;
        *near_index += 1;
        leaf->near_size += 1;
    }
    else{
        double leaf_radius = leaf->xc - leaf->xinterval[0];
        double source_radius = source_panel->xc - source_panel->xinterval[0];
        double distance = std::fabs(leaf->xc - source_panel->xc);
        distance = std::fmin(distance, period-distance);
#if TESTFLAG
        cout << "Leaf Id, Source ID, Leaf radius, source radius, distance, ratio: " << endl;
        cout << leaf->id << "\t" << source_panel->id << "\t" << leaf_radius << "\t" << source_radius << "\t" << distance << "\t" << (leaf_radius + source_radius) / distance << endl;
#endif
        if ( (leaf_radius + source_radius) / distance < MAC){
            leaf->far_ids[*far_index] = source_panel->id;
            *far_index += 1;
            leaf->far_size += 1;
        }
        else{
            init_interaction_lists(leaf, source_panel->left_child, near_index, far_index, period);
            init_interaction_lists(leaf, source_panel->right_child, near_index, far_index, period);
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
    const double eps = 1e-3;
    double z = x - y;
    z = z - round(z);
    return 0.5 * z * sqrt(1 + 4 * eps * eps) / sqrt( z*z + eps*eps  ) - z;
    //return sqrt(x*y);

}

__global__ void computesum(double *e_field, panel *tree_list, int *leaf_indicies, double *source_particles, double *weights, size_t leaf_size){
    int idx = blockIdx.x*blockDim.x + threadIdx.x;

    if(idx >= leaf_size) {return;}

    panel leaf_panel = tree_list[leaf_indicies[idx]];

   for (size_t particlei=0;particlei<leaf_panel.num_members;particlei++){
        double px = source_particles[leaf_panel.members[particlei]];
        for(size_t k=0;k<leaf_panel.far_size;k++){
            panel far_panel = tree_list[leaf_panel.far_ids[k]];
            for (size_t j=0;j<PP;j++){
                e_field[leaf_panel.members[particlei]] += kernel(px, far_panel.s[j]) * far_panel.modified_weights[j];
            }
       } 

        for(size_t k=0;k<leaf_panel.near_size;k++){
            panel near_panel = tree_list[leaf_panel.near_ids[k]];
            for (size_t j=0;j<near_panel.num_members;j++){
                e_field[leaf_panel.members[particlei]] += kernel(px, source_particles[near_panel.members[j]]) * weights[near_panel.members[j]];
            }
        }
    }
    
   printf("e[%d] = %f\n", idx, e_field[idx]);

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
// Set xinterval and xc of root automatically
// Deal with more particles
// Dynamic parallelism

void BLTC(double *e_field, double *source_particles, double *target_particles, double *weights, 
        size_t e_field_size, size_t source_size, size_t target_size){

    const double L = 1.0;

    panel root;
    double xmin = source_particles[0];
    double xmax = source_particles[0];
    for (size_t k=0; k<source_size; k++){
        root.members[k] = k;
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
    for(int k=0;k<Nmax;k++){
        root.near_ids[k] = -1;
        root.far_ids[k] = -1;
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

    for(int k=0;k<leaf_size;k++){
        int near_index = 0;
        int far_index = 0;
        init_interaction_lists(tree_list+leaf_indicies[k], &root, &near_index, &far_index, L); 
    }

#if TESTFLAG
    cout << endl;
    for(int k=0;k<leaf_size;k++){
        cout << "Leaf " << k << " has ID " << leaf_indicies[k] << " and near Ids" << endl;
        int j=0;
        while ( tree_list[leaf_indicies[k]].near_ids[j] != -1){
            cout << tree_list[leaf_indicies[k]].near_ids[j] << "\t";
            j++ ;
        }
        cout << endl;
        cout << "and far Ids" << endl;
        j = 0;
        while ( tree_list[leaf_indicies[k]].far_ids[j] != -1){
            cout << tree_list[leaf_indicies[k]].far_ids[j] << "\t";
            j++;
        }
        cout << endl;
    }
    cout << endl;
#endif

    // Initilize device vars
    panel *d_tree_list;
    double *d_particles;
    double *d_weights;
    double *d_efield;
    int *d_leaf_indicies;

    cudaError_t errcode;

    errcode = cudaMalloc(&d_tree_list, tree_size*sizeof(panel));
    if(checkcudaerr(errcode) != 0){cout << "Failed allocating tree list" << endl;}

    errcode = cudaMalloc(&d_particles, source_size*sizeof(double));
    if(checkcudaerr(errcode) !=0){cout << "Failed allocating source particles" << endl;}

    errcode = cudaMalloc(&d_weights, source_size*sizeof(double));
    if(checkcudaerr(errcode) != 0){cout << "Failed allocating source weights" << endl;}

    errcode = cudaMalloc(&d_efield, target_size*sizeof(double));
    if(checkcudaerr(errcode) != 0){cout << "Failed allocating e_field" << endl;}

    errcode = cudaMalloc(&d_leaf_indicies, leaf_size*sizeof(int));
    if(checkcudaerr(errcode) != 0){cout << "Failed allocating leaf indicies" << endl;}

    errcode = cudaMemcpy(d_tree_list, tree_list, tree_size*sizeof(panel), cudaMemcpyHostToDevice);
    if(checkcudaerr(errcode) != 0){cout << "Failed copying tree list to device" << endl;}
 
    errcode = cudaMemcpy(d_particles, source_particles, source_size*sizeof(double), cudaMemcpyHostToDevice);
    if(checkcudaerr(errcode) != 0){cout << "Failed copying source particles to device" << endl;}
 
    errcode = cudaMemcpy(d_weights, weights, source_size*sizeof(double), cudaMemcpyHostToDevice);
    if(checkcudaerr(errcode) != 0){cout << "Failed copying source weights to device" << endl;}

    errcode = cudaMemcpy(d_leaf_indicies, leaf_indicies, leaf_size*sizeof(int), cudaMemcpyHostToDevice);
    if(checkcudaerr(errcode) != 0){cout << "Failed copying leaf indicies to device" << endl;}

    int blocksize = 1024;
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

    gridlen = (leaf_size + blocksize - 1) / blocksize;
    computesum<<<gridlen,blocksize>>>(d_efield, d_tree_list, d_leaf_indicies, d_particles, d_weights, leaf_size);

    errcode = cudaMemcpy(e_field, d_efield, target_size*sizeof(double), cudaMemcpyDeviceToHost);
    if(checkcudaerr(errcode) != 0){cout << "Failed to copy e_field to host" << endl;}

}
