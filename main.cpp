#include <iostream>
#include <iomanip>
#include <cmath>
#include <cstring>
#include <sys/times.h>
#include <cfloat>
#include <cassert>
#include <fstream>
#include<chrono>
#include<stdlib.h>

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
    const double L = 4*pi;
    const size_t N = argc==2 ? atoi(argv[1]) : 128*128;
    const double dx = L/N;

    cout << "N = " << N << endl;
    cout << "Interp degree = " << P << endl;
    cout << "MAC = " << MAC << endl;
    cout << "L = " << L << endl;
    

    // initialzied root using new 
    double* locs = (double*)malloc(sizeof(double)*N);
    double* source_locs = (double*)malloc(sizeof(double)*N);
    for (size_t k = 0; k < N; k++){
        // Maybe need to handle nondistinct particles
        locs[k] = std::fmod( 0.5 * cos( 2*pi * (k+1) * dx ) + (k+1)*dx, L );
        source_locs[k] = std::fmod( 0.5 * cos( 2*pi * (k+1) * dx ) + (k+1)*dx, L );
        while(locs[k] < 0){ locs[k] += L;  }
        while(source_locs[k] < 0){ source_locs[k] += L;  }
        //cout << k << "\t" << locs[k] << endl;
       // locs[k] = ((int)k)*dx;
    }
    
    double* e_field = (double*)malloc(sizeof(double)*N);
    double* weights = (double*)malloc(sizeof(double)*N);
    for (int k=0;k<N;k++){
        weights[k] = 1.0;
    }


    cout << "Calling BLTC" << endl;
    auto start = high_resolution_clock::now();
    BLTC(e_field, source_locs, locs, weights, N, N, N);
    auto end = high_resolution_clock::now();
    auto duration = duration_cast<milliseconds>(end - start);
    cout << "BLTC time (ms): " << duration.count() << endl;
//#if TESTFLAG
    cout << "Finished BLTC, result is" << endl;
    for(size_t k=0;k<10;k++){
        cout << "e[" << k << "] = " << setprecision(16) << e_field[k] << endl;
    }
//#endif

    double* direct_e_par = (double*)malloc(sizeof(double)*N);

    cout << endl << "Calling direct sum parallel" << endl;
    start = high_resolution_clock::now();
    directsum(direct_e_par, locs, locs, weights, N, N, false);
    end = high_resolution_clock::now();
    duration = duration_cast<milliseconds>(end - start);
    cout << "Direct sum parallel time (ms): " << setprecision(16) << duration.count() << endl;
//#if TESTFLAG
    cout << "Finished direct sum parallel, result is" << endl;
    for(size_t k=0;k<10;k++){
        cout << "e[" << k << "] = " << setprecision(16) << direct_e_par[k] << endl;
    }
//#endif
    
    cout << endl << "Error: ";
    double err = 0.0;
    double perr = 0.0;
    double direct_e_norm = 0.0;
    double maxerr = 0.0;
    size_t maxk = 0;
    for(size_t k=0;k<N;k++){
        perr = (e_field[k] - direct_e_par[k]) * (e_field[k] - direct_e_par[k]);
        if(perr > maxerr){maxerr = perr; maxk = k;}
        err += (e_field[k] - direct_e_par[k]) * (e_field[k] - direct_e_par[k]);
        direct_e_norm += direct_e_par[k] * direct_e_par[k];
    }

    double rel_err = sqrt(err/direct_e_norm);

    cout << setprecision(16) << rel_err << endl;
    cout << "Max error was " << setprecision(16) << sqrt(maxerr) << " at particle " << maxk << endl;
    cout << "BLTC: e[" << maxk << "] = " << setprecision(16) << e_field[maxk] << endl;
    cout << "Direct sum: e[" << maxk << "] = " << setprecision(16) << direct_e_par[maxk] << endl;

    free(locs);
    free(source_locs);
    free(e_field);
    free(weights);
    free(direct_e_par);


/*
    cout << endl << "Calling direct sum parallel (dynamic)" << endl;
    start = high_resolution_clock::now();
    directsum(direct_e_par_dyn, locs, locs, weights, N, N, true);
    end = high_resolution_clock::now();
    duration = duration_cast<milliseconds>(end - start);
    cout << "Direct sum dynamic parallel time (ms): " << duration.count() << endl;
    cout << "Finished direct sum dynamic parallel, result is" << endl;
    for(size_t k=0;k<10;k++){
        cout << "e[" << k << "] = " << setprecision(16) << direct_e_par_dyn[k] << endl;
    }
*/
/*
    if (N <= 1000000){
    cout << endl << "Calling direct sum serial" << endl;
    start = high_resolution_clock::now();
    directsum_serial(direct_e, locs, locs, weights, N, N);
    end = high_resolution_clock::now();
    duration = duration_cast<milliseconds>(end - start);
    cout << "Direct sum serial time (ms): " << setprecision(16) << duration.count() << endl;
//#if TESTFLAG
    cout << "Finished direct sum serial, result is" << endl;
    for(size_t k=0;k<20;k++){
        cout << "e[" << k << "] = " << setprecision(16) << direct_e[k] << endl;
    }
//#endif
    }
*/
    
        
}
