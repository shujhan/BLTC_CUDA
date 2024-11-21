#ifndef BLTC_HPP
#define BLTC_HPP

#include<cstddef>
#include<vector>

// Set numerical paramters
static const size_t numpars_s = 10000;  // num particles
static const size_t N0 = 2; // Max points in a leaf
//static const max_level = 10;
static const double pi = 3.14159265358979323846;
static const int P = 6; // Order of Far-field approximation
static const int PP = P + 1;
static const int Pflat = PP;
static const double sq_theta = 0.36; // theta = 0.6

static const double eps = 0.1;  // Regularization param
//

typedef struct panel panel;

struct panel
{
//    size_t members[2]; // Start and end indicies of contained particles
    std::vector<size_t> members;
    double weights[PP];
    double s[PP];
    double xinterval[2];
    double xc; // Panel center x coordinate
//    panel *children[2];
    panel *left_child;
    panel *right_child;
    panel *parent;
    int level;
    
};

// weights size = target_size
void BLTC(double *e_field, double *source_particles, double *target_particles, double *weights, 
        size_t e_field_size, size_t source_size, size_t target_size);

// Called recursivley from root panel to build tree
void split_panel(panel *p, double* source_particles);

void init_modified_weights(panel *p, double *source_particles, double *weights, size_t source_size);

#endif
