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

    right_child->xinterval[0] = p->xc;
    right_child->xinterval[1] = p->xinterval[1];
    right_child->xc = (p->xinterval[1] + p->xc)/2.0;
    right_child->level = p->level + 1;
    right_child->parent = p;

    for (size_t k=0;k<p->members.size();k++){
        if (source_particles[p->members[k]] < p->xc){
            left_child->members.push_back(p->members[k]);
        }
        else{
            right_child->members.push_back(p->members[k]);
        }
    }


//    cout << left_child->members.size() << "\t" << left_child->xinterval[0] << "\t" << left_child->xinterval[1] << "\t" << left_child->xc << "\t" << left_child->level << "\t" << endl;
//    cout << right_child->members.size() << "\t" << right_child->xinterval[0] << "\t" << right_child->xinterval[1] << "\t" << right_child->xc << "\t" << right_child->level << "\t"  << endl;

    p->left_child = left_child;
    p->right_child = right_child;

    *tree_size += 2;

    if( left_child->members.size() > N0  ){
        split_panel(p->left_child, source_particles, tree_size, leaf_size);
    }
    else{
        *leaf_size += 1;
    }
    if (right_child->members.size() > N0 ){
        split_panel(p->right_child, source_particles, tree_size, leaf_size);
    }
    else{
        *leaf_size += 1;
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
    for (int k = 0; k < p->members.size(); k++) {
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

void init_tree_list(panel *p, vector<panel> tree_list, int *current_id, vector<int> leaf_indicies){
    
    p->id = *current_id;
    tree_list.push_back(*p);
    *current_id += 1;

    if (p->left_child){
        init_tree_list(p->left_child, tree_list, current_id, leaf_indicies);
    }
    if (p->right_child){
        init_tree_list(p->right_child, tree_list, current_id, leaf_indicies);
    }
   // Handle leafs
    if (!(p->left_child) && !(p->right_child)){
        leaf_indicies.push_back(*current_id - 1 );
    }
}
