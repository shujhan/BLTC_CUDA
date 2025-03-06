#ifndef BLTC_HPP
#define BLTC_HPP

#include<cstddef>
#include<cfloat>
#include<cstring>

typedef struct panel panel;

struct panel
{
    size_t members[2]; // indicies of particles in source_parirticles
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
    //int *near_ids;
    //int *far_ids;
    size_t near_size;
    size_t far_size;
    
};

class ElectricField {
    public: 
        virtual void operator()     (double* es, double* targets, int nt, 
                                    double* sources, double* q_ws, int ns) = 0;
        virtual void print_field_obj() = 0;
        virtual ~ElectricField();
};

class E_MQ_DirectSum : public ElectricField {
    public:
        E_MQ_DirectSum();
        E_MQ_DirectSum(double L, double epsilon);
        void operator() (double* es, double* targets, int nt, 
                        double* sources, double* q_ws, int ns) override;
        void print_field_obj();
        ~E_MQ_DirectSum();
    private:
        double epsilon;
        double L;
};

class E_MQ_Treecode : public ElectricField {
    public:
        E_MQ_Treecode();
        // E_MQ_Treecode(double L, double epsilon, double beta);
        E_MQ_Treecode(double L, double epsilon,
                  double mac, int degree, int max_source, int max_target,
                  int verbosity);
        void operator() (double* es, double* targets, int nt, 
                        double* sources, double* q_ws, int ns) override;
        void print_field_obj();
        ~E_MQ_Treecode() override;

    private:

    // Set numerical paramters
    // size_t N;  // num particles
    size_t N0; // Max points in a leaf
    size_t Nmin = 1; // Min points in a leaf
    double pi = 3.14159265358979323846;
    int P; // Order of Far-field approximation
    int PP; // = P + 1
    int Pflat; // = PP
    
    double L;
    double epsilon;
    double MAC;
    int max_source;
    int max_target;

    double sq_theta;
    double norm_epsL;
    double epsLsq;

// weights size = target_size
void BLTC(double *e_field, double *source_particles, double *target_particles, double *weights, 
        size_t e_field_size, size_t source_size, size_t target_size);

// void quicksort(source_particles, (int)source_size, source_indicies);

// Called recursivley from root panel to build tree
void split_panel(panel *p, double* source_particles, int *tree_size, int *leaf_size);

void init_tree_list(panel *p, panel *tree_list, int *current_id, int *leaf_indicies, int *leaf_id);

void init_interaction_lists(panel* leaf, panel* source_panel, int *near_ids, int *far_ids, int *near_index, int *far_index, int leaf_id, int leaf_size, double period);

__global__ void init_modified_weights(panel *d_tree_list, double *d_particles, double *d_weights, int source_size, int tree_size);

__global__ void computesum(double *e_field, panel *tree_list, int *leaf_indicies, double *target_particles, double *source_particles, double *weights, int *d_near_list, int *d_far_list, int leaf_size);

__global__ void computepanelsum(double *e_field, panel *leaf_panel, panel *tree_list, double *target_particles, double *source_particles, double *weights, int *d_near_list, int *d_far_list, int leaf_id, int leaf_size);

__device__ double kernel(double x, double y);

void free_tree_list(panel *panel);

}



































#endif
