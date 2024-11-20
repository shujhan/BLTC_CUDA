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
    cout << p->xinterval[0] << "\t" << p->xinterval[1] << "\t" << p->xc << "\t" << p->level << endl;
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
    
    panel root;
    double locs[N];
    for (size_t k=0;k<N;k++){
        locs[k] = k*dx;//sin( 2*pi/L * k * dx  );
        root.members.push_back(k);
    }

    //root.members[0] = 0;
    //root.members[1] = N-1;
    root.xinterval[0] = 0.0;
    root.xinterval[1] = 1.0;
    root.xc = L/2;
    root.level = 0;

    split_panel(&root, locs, N);


    cout << "interval[0] \t interval[1] \t xc \t level  " << endl;

    treeprint(&root);

}
