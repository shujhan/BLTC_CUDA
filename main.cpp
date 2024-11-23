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

void treeprint(panel *p){
    cout << p->xinterval[0] << "\t" << p->xinterval[1] << "\t" << p->xc << "\t" << p->level << "\t";
    for (size_t k=0; k<p->members.size();k++){
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
        locs[k] = sin( 2*pi/L * k * dx ) + k*dx;
        root.members.push_back(k);
//        cout << locs[k] << "\t";
    }
//    cout << endl;

    //root.members[0] = 0;
    //root.members[1] = N-1;
    root.xinterval[0] = -2;
    root.xinterval[1] = 2;
    root.xc = 0.0;
    root.level = 0;

    split_panel(&root, locs);


//    cout << "interval[0] \t interval[1] \t xc \t level  " << endl;

//    treeprint(&root);

}
