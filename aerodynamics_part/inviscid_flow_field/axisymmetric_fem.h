#ifndef AXISYMMETRIC_FINITE_ELEMENT_H
#define AXISYMMETRIC_FINITE_ELEMENT_H

#include "data_structures_fem.h"
#include "C_FORTRAN_caller.h"

#include <stdio.h>
#include <stddef.h>
#include <stdlib.h>
#include <math.h>
#include <iostream>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif



typedef struct {
    double point;
    double weight;
} GaussQuadraturePoint;

const GaussQuadraturePoint Gauss8Point[8] = {
    { -0.96028985649753623168, 0.10122853629037625915 },
    { -0.79666647741362673959, 0.22238103445337447054 },
    { -0.52553240991632898582, 0.31370664587788728734 },
    { -0.18343464249564980494, 0.36268378337836209114 },
    {  0.18343464249564980494, 0.36268378337836209114 },
    {  0.52553240991632898582, 0.31370664587788728734 },
    {  0.79666647741362673959, 0.22238103445337447054 },
    {  0.96028985649753623168, 0.10122853629037625915 }
};



static inline double triangle_area(double x1, double y1, double x2, double y2, double x3, double y3) {
    return 0.5 * fabs(x1 * (y2 - y3) + x2 * (y3 - y1) + x3 * (y1 - y2));
}





double** compute_local_stiffness_matrix(double x1, double y1, double x2, double y2, double x3, double y3, double *element_area);
double* compute_global_stiffness_matrix(mesh2D *mesh);
void apply_neumann_boundary_conditions(const mesh2D *mesh, double *force_matrix);
void apply_dirichlet_boundary_conditions(const mesh2D *mesh, double **global_stiffness_matrix, double *force_matrix);
void compute_shape_function_contributions_in_neumann_conditions(double x1, double y1, double x2, double y2, double x3, double y3, double flux, double *N1, double *N2, double element_area);





double** compute_local_stiffness_matrix(double x1, double y1, double x2, double y2, double x3, double y3, double *element_area) {
    double c[3], d[3];

    double **local_stiffness_matrix;
    local_stiffness_matrix = (double**)calloc(3, sizeof(double*));
    for (size_t i = 0; i < 3; ++i)
        local_stiffness_matrix[i] = (double*)calloc(3, sizeof(double));

    c[0] = y2 - y3, c[1] = y3 - y1, c[2] = y1 - y2;
    d[0] = x3 - x2, d[1] = x1 - x3, d[2] = x2 - x1;

    double r_centroid = (y1 + y2 + y3) / 3.0;
    (*element_area) = triangle_area(x1, y1, x2, y2, x3, y3);
    double coefficient = M_PI * r_centroid / (2 * (*element_area));

    for (size_t i = 0; i < 3; ++i) {
        for (size_t j = 0; j < 3; ++j) {
            local_stiffness_matrix[i][j] = coefficient * (c[i] * c[j] + d[i] * d[j]);
        }
    }
    return local_stiffness_matrix;
}



double* compute_global_stiffness_matrix(mesh2D *mesh) {
    size_t num_nodes = mesh->num_nodes;

    double *global_stiffness_matrix;
    global_stiffness_matrix = (double*)calloc(num_nodes * num_nodes, sizeof(double*));

    size_t  num_elements = mesh->num_elements;
    for (size_t idx_elem = 0; idx_elem < num_elements; ++idx_elem) {
        size_t idx1 = mesh->triangle_elements[idx_elem].node_inds[0];
        size_t idx2 = mesh->triangle_elements[idx_elem].node_inds[1];
        size_t idx3 = mesh->triangle_elements[idx_elem].node_inds[2];

        double x1 = mesh->nodes[idx1].x, y1 = mesh->nodes[idx1].y;
        double x2 = mesh->nodes[idx2].x, y2 = mesh->nodes[idx2].y;
        double x3 = mesh->nodes[idx3].x, y3 = mesh->nodes[idx3].y;

        double **local_stiffness_matrix = compute_local_stiffness_matrix(x1, y1, x2, y2, x3, y3, &mesh->element_areas[idx_elem]);
        
        for(size_t local_row_idx = 0; local_row_idx < 3; ++local_row_idx) {
            size_t glo_row_idx = mesh->triangle_elements[idx_elem].node_inds[local_row_idx];

            for(size_t local_col_idx = 0; local_col_idx < 3; ++local_col_idx) {
                size_t glo_col_idx = mesh->triangle_elements[idx_elem].node_inds[local_col_idx];
                global_stiffness_matrix[glo_row_idx * num_nodes + glo_col_idx] += local_stiffness_matrix[local_row_idx][local_col_idx];
            }
        }
        for(size_t i = 0; i < 3; ++i) free(local_stiffness_matrix[i]);
        free(local_stiffness_matrix);
    }
    return global_stiffness_matrix;
}



void apply_neumann_boundary_conditions(const mesh2D *mesh, double *force_matrix) {
    size_t n = mesh->num_neumann_edges;
    if (n == 0) return;

    for (size_t i = 0; i < n; ++i) {
        double flux = mesh->neuman_bound[i].flux;
        if (fabs(flux) < 1e-12) continue;
        size_t idx_node_1 = mesh->neuman_bound[i].node_inds[0];
        size_t idx_node_2 = mesh->neuman_bound[i].node_inds[1];

        size_t idx_elem = mesh->neuman_bound[i].element_indx;
        size_t idx_elem_1 = mesh->triangle_elements[idx_elem].node_inds[0];
        size_t idx_elem_2 = mesh->triangle_elements[idx_elem].node_inds[1];
        size_t idx_elem_3 = mesh->triangle_elements[idx_elem].node_inds[2];

        size_t idx_node_3 = 0;
        if (idx_elem_1 == idx_node_1 && idx_elem_2 == idx_node_2) {
            idx_node_3 = idx_elem_3;
        } else if (idx_elem_2 == idx_node_1 && idx_elem_3 == idx_node_2) {
            idx_node_3 = idx_elem_1;
        } else if (idx_elem_3 == idx_node_1 && idx_elem_1 == idx_node_2) {
            idx_node_3 = idx_elem_2;
        }

        double x1 = 0.0, y1 = 0.0, x2 = 0.0, y2 = 0.0, x3 = 0.0, y3 = 0.0;
        x1 = mesh->nodes[idx_node_1].x, y1 = mesh->nodes[idx_node_1].y;
        x2 = mesh->nodes[idx_node_2].x, y2 = mesh->nodes[idx_node_2].y;
        x3 = mesh->nodes[idx_node_3].x, y3 = mesh->nodes[idx_node_3].y;
        
        double N1, N2;
        compute_shape_function_contributions_in_neumann_conditions(x1, y1, x2, y2, x3, y3, flux, &N1, &N2, mesh->element_areas[idx_elem]);

        force_matrix[idx_node_1] += N1;
        force_matrix[idx_node_2] += N2;
    }
}



void apply_dirichlet_boundary_conditions(const mesh2D *mesh, double **global_stiffness_matrix, double *force_matrix) {
    size_t n = mesh->num_dirichlet_nodes;
    if (n == 0) return;

    for (size_t i = 0; i < n; ++i) {
        size_t idx_col = mesh->dirichlet_inds[i];
        double node_vals = mesh->dirichlet_bound[i];
        double m = mesh->num_nodes;

        for (size_t j = 0; j < m; ++j) {
            force_matrix[j] -= global_stiffness_matrix[j][idx_col] * node_vals;
            global_stiffness_matrix[j][idx_col] = 0.0;
            global_stiffness_matrix[idx_col][j] = 0.0;
        }
        global_stiffness_matrix[idx_col][idx_col] = 1.0;   
        force_matrix[idx_col] = node_vals;         
    }
}



void compute_shape_function_contributions_in_neumann_conditions(double x1, double y1, double x2, double y2, double x3, double y3, double flux, double *N1, double *N2, double element_area) {
    *N1 = 0.0, *N2 = 0.0;

    double b1 = x2 * y3 - x3 * y2, b2 = x3 * y1 - x1 * y3;
    double c1 = y2 - y3, c2 = y3 - y1;
    double d1 = x3 - x2, d2 = x1 - x3;

    double edge_length = sqrt((x1 - x2) * (x1 - x2) + (y1 - y2) * (y1 - y2));

    for (size_t i = 0; i < 8; ++i) {
        double xi = Gauss8Point[i].point;
        double weight = Gauss8Point[i].weight;

        double x_ = (1 - xi) / 2.0 * x1 + (1 + xi) / 2.0 * x2;
        double y_ = (1 - xi) / 2.0 * y1 + (1 + xi) / 2.0 * y2;
        *N1 += (b1 + c1 * x_ + d1 * y_) * y_ * weight;
        *N2 += (b2 + c2 * x_ + d2 * y_) * y_ * weight;
    }
    double factor = M_PI * flux / element_area;
    *N1 *= factor * edge_length / 2.0;
    *N2 *= factor * edge_length / 2.0;
}



double* laplace_solver(mesh2D *mesh) {
    int n = mesh->num_nodes;
    double *global_stiffness_matrix = compute_global_stiffness_matrix(mesh);
    double *force_matrix = (double*)calloc(n, sizeof(double));
    apply_neumann_boundary_conditions(mesh, force_matrix);
    double *solution = (double*)calloc(n, sizeof(double));
    solve_sparse_system(global_stiffness_matrix, force_matrix, solution, &n);
    free(global_stiffness_matrix);
    free(force_matrix);
    return solution;
}
#endif