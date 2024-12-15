#include <iostream>
#include <iomanip>
#include <cmath>
#include <cstring>
#include <sys/times.h>
#include <cfloat>
#include <cassert>
#include <fstream>
#include<chrono>

using namespace std;
using namespace std::chrono;

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
        locs[k] = std::fmod( 0.5 * cos( 2*pi/L * k * dx ) + k*dx, L );
        while(locs[k] < 0){ locs[k] += L;  }
        cout << k << "\t" << locs[k] << endl;
       // locs[k] = ((int)k)*dx;
    }
    
    double e_field[N];
    double weights[N];
    for (int k=0;k<N;k++){
        weights[k] = 1.0;
    }


    cout << "Calling BLTC" << endl;
    auto start = high_resolution_clock::now();
    BLTC(e_field, locs, locs, weights, N, N, N);
    auto end = high_resolution_clock::now();
    auto duration = duration_cast<milliseconds>(end - start);
    cout << "BLTC time (ms): " << duration.count() << endl;
    cout << "Finished BLTC, result is" << endl;
    for(size_t k=0;k<N;k++){
        cout << "e[" << k << "] = " << setprecision(16) << e_field[k] << endl;
    }

    double direct_e[N];
    double direct_e_par[N];

    cout << endl << "Calling direct sum" << endl;
    start = high_resolution_clock::now();
    directsum(direct_e_par, locs, locs, weights, N, N);
    end = high_resolution_clock::now();
    duration = duration_cast<milliseconds>(end - start);
    cout << "Direct sum parallel time (ms): " << duration.count() << endl;
    cout << "Finished direct sum, result is" << endl;
    for(size_t k=0;k<N;k++){
        cout << "e[" << k << "] = " << setprecision(16) << direct_e_par[k] << endl;
    }

    cout << endl << "Calling direct sum serial" << endl;
    start = high_resolution_clock::now();
    directsum_serial(direct_e, locs, locs, weights, N, N);
    end = high_resolution_clock::now();
    duration = duration_cast<milliseconds>(end - start);
    cout << "Direct sum serial time (ms): " << duration.count() << endl;
    cout << "Finished direct sum serial, result is" << endl;
    for(size_t k=0;k<N;k++){
        cout << "e[" << k << "] = " << setprecision(16) << direct_e[k] << endl;
    }

    cout << endl << "Error: ";
    double err = 0.0;
    double direct_e_norm = 0.0;
    for(size_t k=0;k<N;k++){
        err += (e_field[k] - direct_e[k]) * (e_field[k] - direct_e[k]);
        direct_e_norm += direct_e[k] * direct_e[k];
    }

    double rel_err = sqrt(err/direct_e_norm);

    cout << setprecision(16) << rel_err << endl;
}
