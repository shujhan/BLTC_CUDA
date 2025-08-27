
#ifndef FIELD_STRUCTURE_HPP
#define FIELD_STRUCTURE_HPP

#include<cstddef>
#include<cfloat>
#include<cstring>

using std::cout;
using std::endl;

using namespace std;

// Helper struct definition for panel
struct panel
{
    size_t members[2];
    double xinterval[2];
    double xc; // Panel center x coordinate
    double MAC; // r^2 / theta^2
    std::vector<size_t> children;

    // Initialization
    panel() : xc(0.0), MAC(0.0)
    {
        std::memset(members, 0, sizeof(members));
        std::memset(xinterval, 0, sizeof(xinterval));
    }
};

class ElectricField {
    public: 
        virtual void operator()     (double* es, double* targets, int nt, 
                                    double* sources, double* q_ws, int ns) = 0;
        virtual ~ElectricField();
};

class E_MQ_DirectSum : public ElectricField {
    public:
        E_MQ_DirectSum();
        E_MQ_DirectSum(double L, double epsilon);
        double epsilon;
        double L;
        void operator() (double *e_field, double *source_particles, double *target_particles, double *weights,
        size_t source_size, size_t target_size, bool nested);
        ~E_MQ_DirectSum();


    __global__ void direct_e_sum(double *d_efield, double *d_particles, double *d_target, double *d_weights,
        const size_t source_size, size_t target_size, double3 kernel_params);

    __device__ inline double kernelp(double x, double y, double3 kernel_params);
};

class E_MQ_DirectSum_openbcs : public ElectricField {
    public:
        E_MQ_DirectSum_openbcs();
        E_MQ_DirectSum_openbcs(double epsilon);
        double epsilon;
        void operator() (double* es, double* targets, int nt, 
                        double* sources, double* q_ws, int ns);
        ~E_MQ_DirectSum_openbcs();
};

class E_MQ_Treecode : public ElectricField {
    public:
        E_MQ_Treecode();
        E_MQ_Treecode(double L, double epsilon, double beta);
        E_MQ_Treecode(double L, double epsilon,
                  double mac, int degree, int max_source, int max_target,
                  int verbosity);
        void operator() (double* es, double* targets, int nt, 
                        double* sources, double* q_ws, int ns) override;
        void print_field_obj();
        ~E_MQ_Treecode();

    // private:
    void BLTC(double *e_field, double *source_particles, double *target_particles, double *weights, 
        size_t e_field_size, size_t source_size, size_t target_size);

    void split_panel(panel *p, double *source_particles, int *tree_size, int *leaf_size, size_t *indicies);
    void init_tree_list(panel *p, panel *tree_list, int *current_id, int *leaf_indicies, int *leaf_id);
    void free_tree_list(panel *panel);
    __global__ void init_modified_weights(panel *d_tree_list, double *d_particles, double *d_weights, int source_size, int tree_size);
    void init_interaction_lists(panel* leaf, panel* source_panel, int *near_ids, int *far_ids, int *near_index, int *far_index, int leaf_id, int leaf_size, double period, int *total_nears);
    __device__ inline double kernel(double x, double y, double3 kernel_params); 
    __global__ void computesum(double *e_field, panel *tree_list, int *leaf_indicies, double *target_particles, double *source_particles, double *weights, int *d_near_list, int *d_far_list, int leaf_size);
    __global__ void computepanelsum(double *e_field, panel *leaf_panel, panel *tree_list, double *target_particles, double *source_particles, double *weights, int *d_near_list, int *d_far_list, int leaf_id, int leaf_size);
    __global__ void computepanelsum_far(double *e_field, panel *leaf_panel, panel *tree_list, double *target_particles, double *source_particles, int *d_far_list, int leaf_id, int leaf_size, double3 kernel_params);
    __global__ void computepanelsum_near(double *e_field, double *target_particles, double *source_particles, double *weights, uint4* near_data, int total_nears, double3 kernel_params);
    __global__ void re_order(double* e_field, double* ordered_efield, size_t* indicies, int source_size);

    int checkcudaerr(cudaError_t err);


};

class E_MQ_Treecode_openbcs : public ElectricField {
    public:
        // KERNEL kernel;
        // std::vector<double> kernelParams;
        // SINGULARITY singularity;
        // APPROXIMATION approximation;
        // COMPUTE_TYPE compute_type;
        // double beta, theta;
        // int interpDegree, maxPerSourceLeaf, maxPerTargetLeaf;
        // int verbosity;

        // E_MQ_Treecode_openbcs();
        // E_MQ_Treecode_openbcs(double epsilon, double beta);
        // E_MQ_Treecode_openbcs(double epsilon,
        //     double theta, int interpDegree, int maxPerSourceLeaf, int maxPerTargetLeaf,
        //     int verbosity);
        // void operator() (double* es, double* targets, int nt, 
        //                 double* sources, double* q_ws, int ns);
        // void print_field_obj();
        // ~E_MQ_Treecode_openbcs();
        E_MQ_Treecode_openbcs();
        ~E_MQ_Treecode_openbcs();


};
#endif /* FIELD_STRUCTURE_HPP */
