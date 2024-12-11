#include <iostream>
#include <iomanip>
#include <cmath>
#include <cstring>
#include <sys/times.h>
#include <vector>
#include <cfloat>
#include <cassert>
#include <fstream>

using namespace std;

#include "BLTC.hpp"
#include "directsum.hpp"

void treeprint(panel *p){
    cout << p->xinterval[0] << "\t" << p->xinterval[1] << "\t" << p->xc << "\t" << p->level << "\t";
    for (size_t k=0; k<p->num_members;k++){
        cout << p->members[k] << "\t";
    }
    cout << endl;
    if (p->left_child){
        treeprint(p->left_child);
    }
    if (p->right_child){
        treeprint(p->right_child);
    }
}

int main(int argc, char** argv) {
    const double L = 1.0;
    const size_t N = 10; 
    const double dx = L/N;
    

    // initialzied root using new 
    panel root;
    double locs[N];
    for (size_t k = 0; k < N; k++){
        // Maybe need to handle nondistinct particles
        //locs[k] = sin( 2*pi/L * k * dx ) + k*dx;
        locs[k] = ((int)k)*dx;
        root.members[k] = k;
//        cout << locs[k] << "\t";
    }
//    cout << endl;

    //root.members[0] = 0;
    //root.members[1] = N-1;
    root.xinterval[0] = 0;
    root.xinterval[1] = L;
    root.xc = L/2;
    root.level = 0;
    root.num_members = N;

    int tree_size = 1;
    int leaf_size = 0;

    cout << "Calling split panel" << endl;

    split_panel(&root, locs, &tree_size, &leaf_size);

    cout << "Tree size = " << tree_size << endl << " Leaf size = " << leaf_size << endl;


//    cout << "interval[0] \t interval[1] \t xc \t level  " << endl;

//    treeprint(&root);
    
    double e_field[N];
    double weights[N];
    for (int k=0;k<N;k++){
        weights[k] = 1.0;
    }


    cout << "Calling BLTC" << endl;
    BLTC(e_field, locs, locs, weights, N, N, N);
    cout << "Finished BLTC" << endl;

    cout << endl << "Calling direct sum" << endl;
    directsum(e_field, locs, locs, weights, N, N);
    cout << "Finished direct sum, result is" << endl;
    for(size_t k=0;k<N;k++){
        cout << "e[" << k << "] = " << e_field[k] << endl;
    }

    cout << endl << "Calling direct sum serial" << endl;
    directsum_serial(e_field, locs, locs, weights, N, N);
    cout << "Finished direct sum serial, result is" << endl;
    for(size_t k=0;k<N;k++){
        cout << "e[" << k << "] = " << e_field[k] << endl;
    }

}
