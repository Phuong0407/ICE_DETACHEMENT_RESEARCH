#ifndef CSR_SPARSE_SOLVER_H
#define CSR_SPARSE_SOLVER_H

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define TOL 1.0E-6
#define MAX_ITER 100000

typedef struct {
    int size;            // Number of rows/columns
    int nnz;             // Number of nonzero elements
    double *val;         // Nonzero values
    int *col_idx;        // Column indices
    int *row_ptr;        // Row pointers
} CSRMatrix;





void init_sparse_solver(CSRMatrix *csr, const double *matrix, int n);
void mult_sparse_matrix(const CSRMatrix *csr, const double *x, double *y);
void solve_sparse_system_from_rhs(const CSRMatrix *csr, const double *rhs, double *solution);
void solve_sparse_system(const double *matrix, const double *rhs, double *solution, int n);
void free_csr(CSRMatrix *csr);





void init_sparse_solver(CSRMatrix *csr, const double *matrix, int n) {
    int i, j, index, count = 0;
    double matrix_norm = 0.0, zero_tol;

    for (i = 0; i < n; i++) {
        double row_sum = 0.0;
        for (j = 0; j < n; j++) {
            row_sum += fabs(matrix[i * n + j]);
        }
        if (row_sum > matrix_norm) {
            matrix_norm = row_sum;
        }
    }
    zero_tol = matrix_norm * __DBL_EPSILON__;

    for (i = 0; i < n * n; i++) {
        if (fabs(matrix[i]) > zero_tol) count++;
    }
    csr->size = n;
    csr->nnz = count;

    csr->val = (double *)malloc(count * sizeof(double));
    csr->col_idx = (int *)malloc(count * sizeof(int));
    csr->row_ptr = (int *)malloc((n + 1) * sizeof(int));

    index = 0;
    csr->row_ptr[0] = 0;
    for (i = 0; i < n; i++) {
        for (j = 0; j < n; j++) {
            if (fabs(matrix[i * n + j]) > zero_tol) {
                csr->val[index] = matrix[i * n + j];
                csr->col_idx[index] = j;
                index++;
            }
        }
        csr->row_ptr[i + 1] = index;
    }
}

void mult_sparse_matrix(const CSRMatrix *csr, const double *x, double *y) {
    int i, j;
    for (i = 0; i < csr->size; i++) {
        y[i] = 0.0;
        for (j = csr->row_ptr[i]; j < csr->row_ptr[i + 1]; j++) {
            y[i] += csr->val[j] * x[csr->col_idx[j]];
        }
    }
}

void solve_sparse_system_from_rhs(const CSRMatrix *csr, const double *rhs, double *solution) {
    int i, iter;
    double alpha, beta, rnorm, bnorm, tol = TOL;
    int n = csr->size;

    double *r = (double *)malloc(n * sizeof(double));
    double *p = (double *)malloc(n * sizeof(double));
    double *Ap = (double *)malloc(n * sizeof(double));

    // Initialize
    for (i = 0; i < n; i++) solution[i] = 0.0;
    
    bnorm = 0.0;
    for (i = 0; i < n; i++) bnorm += rhs[i] * rhs[i];
    bnorm = sqrt(bnorm);

    mult_sparse_matrix(csr, solution, Ap);
    for (i = 0; i < n; i++) {
        r[i] = rhs[i] - Ap[i];
        p[i] = r[i];
    }
    rnorm = 0.0;
    for (i = 0; i < n; i++) rnorm += r[i] * r[i];

    iter = 0;
    while ((sqrt(rnorm) / bnorm > tol) && (iter < MAX_ITER)) {
        iter++;
        mult_sparse_matrix(csr, p, Ap);

        double pAp = 0.0;
        for (i = 0; i < n; i++) pAp += p[i] * Ap[i];

        if (pAp == 0.0) break;

        alpha = rnorm / pAp;
        for (i = 0; i < n; i++) solution[i] += alpha * p[i];
        for (i = 0; i < n; i++) r[i] -= alpha * Ap[i];

        double new_rnorm = 0.0;
        for (i = 0; i < n; i++) new_rnorm += r[i] * r[i];

        if (sqrt(new_rnorm) / bnorm < tol) break;

        beta = new_rnorm / rnorm;
        for (i = 0; i < n; i++) p[i] = r[i] + beta * p[i];

        rnorm = new_rnorm;
    }

    if (iter == MAX_ITER) {
        printf("Failed to converge in %d iterations in solve_sparse_system.\n", MAX_ITER);
    }

    free(r);
    free(p);
    free(Ap);
}

void solve_sparse_system(const double *matrix, const double *rhs, double *solution, int n) {
    CSRMatrix csr;
    init_sparse_solver(&csr, matrix, n);
    solve_sparse_system_from_rhs(&csr, rhs, solution);
    free_csr(&csr);
}

void free_csr(CSRMatrix *csr) {
    free(csr->val);
    free(csr->col_idx);
    free(csr->row_ptr);
}

#endif