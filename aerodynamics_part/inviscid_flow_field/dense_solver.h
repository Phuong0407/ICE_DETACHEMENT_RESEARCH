#ifndef DENSE_SOLVER_H
#define DENSE_SOLVER_H

#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#define TINY 1.5E-16

void decompose_LU(double *A, int *idx, int *parity, int *is_singular, int n) {
    int i, j, k, i_max;
    double sum, row_max_val, pivot, pivot_comparator, row_swap, pivot_divisor;
    double *implicit_scaling_factor = (double*)malloc(n * sizeof(double));

    if (!implicit_scaling_factor) {
        fprintf(stderr, "Memory allocation failed in decompose_LU.\n");
        exit(1);
    }

    *parity = 1;
    *is_singular = 0;

    for (i = 0; i < n; i++) {
        row_max_val = 0.0;
        for (j = 0; j < n; j++) {
            if (fabs(A[i * n + j]) > row_max_val)
                row_max_val = fabs(A[i * n + j]);
        }
        if (row_max_val < TINY) {
            *is_singular = 1;
            free(implicit_scaling_factor);
            return;
        }
        implicit_scaling_factor[i] = 1.0 / row_max_val;
    }

    for (j = 0; j < n; j++) {
        for (i = 0; i < j; i++) {
            sum = A[i * n + j];
            for (k = 0; k < i; k++)
                sum -= A[i * n + k] * A[k * n + j];
            A[i * n + j] = sum;
        }

        pivot = 0.0;
        for (i = j; i < n; i++) {
            sum = A[i * n + j];
            for (k = 0; k < j; k++)
                sum -= A[i * n + k] * A[k * n + j];
            A[i * n + j] = sum;
            pivot_comparator = implicit_scaling_factor[i] * fabs(sum);
            if (pivot_comparator >= pivot) {
                i_max = i;
                pivot = pivot_comparator;
            }
        }

        if (j != i_max) {
            for (k = 0; k < n; k++) {
                row_swap = A[i_max * n + k];
                A[i_max * n + k] = A[j * n + k];
                A[j * n + k] = row_swap;
            }
            *parity = -(*parity);
            implicit_scaling_factor[i_max] = implicit_scaling_factor[j];
        }

        idx[j] = i_max;
        if (fabs(A[j * n + j]) < TINY)
            A[j * n + j] = TINY;

        if (j != n - 1) {
            pivot_divisor = 1.0 / A[j * n + j];
            for (i = j + 1; i < n; i++)
                A[i * n + j] *= pivot_divisor;
        }
    }

    free(implicit_scaling_factor);
}

void substitute_backward_LU(double *A, int *idx, double *B, int n) {
    int i, j, first_non_null_idx = -1;
    double sum;

    for (i = 0; i < n; i++) {
        int row_perm_idx = idx[i];
        sum = B[row_perm_idx];
        B[row_perm_idx] = B[i];
        if (first_non_null_idx >= 0) {
            for (j = first_non_null_idx; j < i; j++)
                sum -= A[i * n + j] * B[j];
        } else if (sum != 0.0) {
            first_non_null_idx = i;
        }
        B[i] = sum;
    }

    for (i = n - 1; i >= 0; i--) {
        sum = B[i];
        for (j = i + 1; j < n; j++)
            sum -= A[i * n + j] * B[j];
        B[i] = sum / A[i * n + i];
    }
}

void solve_dense_system(double *A, double *B, double **solution, int n) {
    if (*solution == NULL) {
        *solution = (double*)malloc(n * sizeof(double));
        if (!(*solution)) {
            fprintf(stderr, "Memory allocation failed in solve_dense_system.\n");
            exit(1);
        }
    }

    int *idx = (int*)malloc(n * sizeof(int));
    int is_singular, parity;

    if (!idx) {
        fprintf(stderr, "Memory allocation failed in solve_dense_system.\n");
        exit(1);
    }

    for (int i = 0; i < n; i++)
        (*solution)[i] = B[i];

    decompose_LU(A, idx, &parity, &is_singular, n);
    if (is_singular == 0)
        substitute_backward_LU(A, idx, *solution, n);

    free(idx);
}

#endif