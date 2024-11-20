#include "BLTC.hpp"

#include <iostream>
using std::cout, std::endl;

void split_panel(panel *p, double* source_particles, size_t source_size){

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

    for (size_t k=0;k<source_size;k++){
        if (source_particles[k] < p->xc){
            left_child->members.push_back(k);
        }
        else{
            right_child->members.push_back(k);
        }
    }


/*    size_t member_cutoff;

    cout << p->xc << endl;
    for(size_t k=p->members[0];k<=p->members[1];k++){
        cout << source_particles[k] << "\t" << p->xc << endl;
        if (source_particles[k] >= p->xc){
            member_cutoff = k-1;
            break;
        }
    }

    cout << "Member cutoff: " << member_cutoff << endl;
    
    left_child->members[0] = p->members[0];
    left_child->members[1] = member_cutoff;
    right_child->members[0] = member_cutoff+1;
    right_child->members[1] = p->members[1];
*/

    cout << left_child->members.size() << "\t" << left_child->xinterval[0] << "\t" << left_child->xinterval[1] << "\t" << left_child->xc << "\t" << left_child->level << "\t" << endl;
    cout << right_child->members.size() << "\t" << right_child->xinterval[0] << "\t" << right_child->xinterval[1] << "\t" << right_child->xc << "\t" << right_child->level << "\t"  << endl;

//    p->children[0] = &(child[0]);
//    p->children[1] = &(child[1]);
    p->left_child = left_child;
    p->right_child = right_child;

/*
    cout << member_cutoff - p->members[0] << endl;
    cout << p->members[1] - member_cutoff << endl;
    if ( (member_cutoff != p->members[0]) &&  (member_cutoff - p->members[0] + 1) > N0){
        cout << endl << "Beginning left interval" << endl << endl;;
        split_panel(p->left_child, source_particles);
    }
    if ( (p->members[1] != member_cutoff) &&  (p->members[1] - member_cutoff) > N0){
        cout << endl << "Beginning right interval" << endl << endl;
        split_panel(p->right_child, source_particles);
    }
*/

    if( left_child->members.size() > N0  ){
        split_panel(p->left_child, source_particles, source_size);
    }
    if (right_child->members.size() > N0 ){
        split_panel(p->right_child, source_particles, source_size);
    }

}

