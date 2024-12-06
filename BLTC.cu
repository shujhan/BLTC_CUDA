#include "BLTC.hpp"
#include<cmath>

#include <iostream>
#include <iomanip>
using std::cout, std::endl;
using namespace std;

void split_panel(panel *p, double* source_particles, int *tree_size, int *leaf_size){

    panel *left_child = new panel();
    panel *right_child = new panel();

    left_child->xinterval[0] = p->xinterval[0];
    left_child->xinterval[1] = p->xc;
    left_child->xc = (p->xinterval[0] + p->xc)/2.0;
    left_child->level = p->level + 1;
    left_child->parent = p;
    left_child->num_members = 0;

    right_child->xinterval[0] = p->xc;
    right_child->xinterval[1] = p->xinterval[1];
    right_child->xc = (p->xinterval[1] + p->xc)/2.0;
    right_child->level = p->level + 1;
    right_child->parent = p;
    right_child->num_members = 0;


    cout << "Splitting particles with " << p->num_members << " Members" << endl;
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
    cout << " Split panels sucessfully, found " << left_child->num_members << " left particles and " << right_child->num_members << " right particles" << endl;


//    cout << left_child->members.size() << "\t" << left_child->xinterval[0] << "\t" << left_child->xinterval[1] << "\t" << left_child->xc << "\t" << left_child->level << "\t" << endl;
//    cout << right_child->members.size() << "\t" << right_child->xinterval[0] << "\t" << right_child->xinterval[1] << "\t" << right_child->xc << "\t" << right_child->level << "\t"  << endl;

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

    if (idx > tree_size){return;}

    panel *p = d_tree_list + idx;

    // Calculate Chebyshev points
    for (int k=0;k < PP;k++){
        p->s[k] = 0.5 * ( p->xinterval[0] + p->xinterval[1] + std::cos(k*pi/PP)*( p->xinterval[1] - p->xinterval[0]  ));
    }

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

void init_modified_weights(panel *p, double *source_particles, double *weights, size_t source_size){
    // Calculate Chebyshev points
    for (int k=0;k < PP;k++){
        p->s[k] = 0.5 * ( p->xinterval[0] + p->xinterval[1] + std::cos(k*pi/PP)*( p->xinterval[1] - p->xinterval[0]  ));
    }

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
        double y = source_particles[p->members[k]]; // particles in cluster
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
            p->modified_weights[i] += a1[i] * D * weights[p->members[i]];
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

void push_to_device(panel *tree, double *source_particles, double *source_weights, int num_panels, int num_particles, panel *d_tree_list, double *d_particles, double *d_weights){
   cudaError_t errcode;

   errcode = cudaMalloc(&d_tree_list, num_panels*sizeof(panel));
   if (errcode != cudaSuccess){
        cout << "Failed to allocate tree on device with code" << errcode << endl;
   }
   errcode = cudaMalloc(&d_particles, num_particles*sizeof(double));
   if (errcode != cudaSuccess){
        cout << "Failed to allocate particles on device with code" << errcode << endl;
   }
   errcode = cudaMalloc(&d_weights, num_particles*sizeof(double));
   if (errcode != cudaSuccess){
        cout << "Failed to allocate weights on device with code" << errcode << endl;
   }

   errcode = cudaMemcpy(d_tree_list, tree, num_panels*sizeof(panel), cudaMemcpyHostToDevice);
   if (errcode != cudaSuccess){
        cout << "Failed to transfer tree to device with code" << errcode << endl;
   }

   errcode = cudaMemcpy(d_particles, source_particles, num_particles*sizeof(double), cudaMemcpyHostToDevice);
   if (errcode != cudaSuccess){
        cout << "Failed to transfer particles to device with code" << errcode << endl;
   }

   errcode = cudaMemcpy(d_weights, source_weights, num_particles*sizeof(double), cudaMemcpyHostToDevice);
   if (errcode != cudaSuccess){
        cout << "Failed to transfer particle weights to device with code" << errcode << endl;
   }
}

void BLTC(double *e_field, double *source_particles, double *target_particles, double *weights, 
        size_t e_field_size, size_t source_size, size_t target_size){

    const double L = 1.0;

    panel root;
    for (size_t k=0; k<source_size; k++){
        root.members[k] = k;
    }
    root.xinterval[0] = 0.0;
    root.xinterval[1] = L;
    root.xc = L/2;
    root.level = 0;
    root.num_members = source_size;

    int tree_size = 1;
    int leaf_size = 0;

    split_panel(&root, source_particles, &tree_size, &leaf_size);

    panel tree_list[tree_size];
    int leaf_indicies[leaf_size];

    int id = 0;
    int leaf_id = 0;

    init_tree_list(&root, tree_list, &id, leaf_indicies, &leaf_id);

    // Initilize device vars
    panel *d_tree_list;
    double *d_particles;
    double *d_weights;

    push_to_device(tree_list, source_particles, weights, tree_size, source_size, d_tree_list, d_particles, d_weights);

}
