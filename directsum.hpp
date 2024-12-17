//void kernel(double *res, double x, double y);

//void direct_e_sum(double *d_efield, double *d_particles, double *d_target, double *d_weights,
//        size_t source_size, size_t target_size);

void directsum(double *e_field, double *source_particles, double *target_particles, double *weights,
        size_t source_size, size_t target_size, bool dynamic);

void directsum_serial(double *e_field, double *source_particles, double *target_particles, double *weights,
        size_t source_size, size_t target_size);

//double kernel_serial(double x, double y);
