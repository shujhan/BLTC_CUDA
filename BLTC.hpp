#ifndef BLTC_HPP
#define BLTC_HPP

#include<cstddef>
#include<cfloat>
#include<cstring>

// Set numerical paramters
static const size_t numpars_s = 10000;  // num particles
static const size_t N0 = 2; // Max points in a leaf
static const size_t Nmax = 10; // Max particles in a panel
//static const max_level = 10;
static const double pi = 3.14159265358979323846;
static const int P = 6; // Order of Far-field approximation
static const int PP = P + 1;
static const int Pflat = PP;
static const double MAC = 0.6;

static const double eps = 0.1;  // Regularization param
//

typedef struct panel panel;

struct panel
{
    size_t members[Nmax]; // indicies of particles in source_parirticles
    size_t num_members;
    double modified_weights[PP]; // modified weights
    double s[PP]; // mapped Chebyshev points in the cluster
    double xinterval[2]; // start and end x position
    double xc; // Panel center x coordinate
    panel *left_child; // These pointers will not work on the device b/c they point to host versions
    panel *right_child;
    panel *parent;
    int level;
    int id;
    // Assign dynamically in future
    int near_ids[Nmax];
    int far_ids[Nmax];
    
};

// weights size = target_size
void BLTC(double *e_field, double *source_particles, double *target_particles, double *weights, 
        size_t e_field_size, size_t source_size, size_t target_size);

// Called recursivley from root panel to build tree
void split_panel(panel *p, double* source_particles, int *tree_size, int *leaf_size);

void init_modified_weights(panel *p, double *source_particles, double *weights, size_t source_size);

//void init_tree_list(panel *p, std::vector<panel> tree_list, int *current_id, std::vector<int> leaf_indicies);

#endif
