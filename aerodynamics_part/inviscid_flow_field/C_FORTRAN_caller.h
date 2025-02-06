#ifndef C_FORTRAN_CALLER_H
#define C_FORTRAN_CALLER_H

#include <stdlib.h>

#ifdef __cplusplus
extern "C" {
#endif

extern void solve_sparse_system(double* sparse_matrix, double* rhs, double* solution, int *n);
extern void solve_dense_system(double *A, double *B, double* solution, int *n);

#ifdef __cplusplus
}
#endif

#endif