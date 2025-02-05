#ifndef C_FORTRAN_CALLER_H
#define C_FORTRAN_CALLER_H

#ifdef __cplusplus
extern "C" {
#endif

extern void init_sparse_solver(double *sparse_matrix, int *n);
extern void solve_sparse_system(double *rhs_vect, double *solution, int *n);

#ifdef __cplusplus
}
#endif

#endif