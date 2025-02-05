#ifndef AXISYMMETRIC_FINITE_ELEMENT_H
#define AXISYMMETRIC_FINITE_ELEMENT_H

#include "geometry.h"

#include <stdio.h>
#include <stddef.h>
#include <stdlib.h>
#include <math.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

static inline double triangle_area(double x1, double y1, double x2, double y2, double x3, double y3) {
    return 0.5 * fabs(x1 * (y2 - y3) + x2 * (y3 - y1) + x3 * (y1 - y2));
}





double** compute_local_stiffness_matrix(double x1, double y1, double x2, double y2, double x3, double y3);
double** compute_global_stiffness_matrix(mesh2D *mesh);
// double** compute_neumann_boundary_conditions(mesh2D *mesh);





double** compute_local_stiffness_matrix(double x1, double y1, double x2, double y2, double x3, double y3) {
    double c[3], d[3];

    double **local_stiffness_matrix;
    local_stiffness_matrix = (double**)calloc(3, sizeof(double*));
    for (size_t i = 0; i < 3; ++i)
        local_stiffness_matrix[i] = (double*)calloc(3, sizeof(double));

    c[0] = y2 - y3, c[1] = y3 - y1, c[2] = y1 - y2;
    d[0] = x3 - x2, d[1] = x1 - x3, d[2] = x2 - x1;

    double r_centroid = (y1 + y2 + y3) / 3.0;
    double A = triangle_area(x1, y1, x2, y2, x3, y3);
    double coefficient = M_PI * r_centroid / (2 * A);

    for (size_t i = 0; i < 3; ++i) {
        for (size_t j = 0; j < 3; ++j) {
            local_stiffness_matrix[i][j] = coefficient * (c[i] * c[j] + d[i] * d[j]);
        }
    }
    return local_stiffness_matrix;
}



double** compute_global_stiffness_matrix(mesh2D *mesh) {
    size_t num_nodes = mesh->num_nodes;

    double **global_stiffness_matrix;
    global_stiffness_matrix = (double**)calloc(num_nodes, sizeof(double*));
    for (size_t i = 0; i < num_nodes; ++i)
        global_stiffness_matrix[i] = (double*)calloc(num_nodes, sizeof(double));

    size_t  num_elements = mesh->num_elements;
    for (size_t idx_elem = 0; idx_elem < num_elements; ++idx_elem) {
        size_t idx1 = mesh->triangle_elements[idx_elem].node_inds[0];
        size_t idx2 = mesh->triangle_elements[idx_elem].node_inds[1];
        size_t idx3 = mesh->triangle_elements[idx_elem].node_inds[2];

        double x1 = mesh->nodes[idx1].x, y1 = mesh->nodes[idx1].y;
        double x2 = mesh->nodes[idx2].x, y2 = mesh->nodes[idx2].y;
        double x3 = mesh->nodes[idx3].x, y3 = mesh->nodes[idx3].y;

        double **local_stiffness_matrix = compute_local_stiffness_matrix(x1, y1, x2, y2, x3, y3);
        
        for(size_t local_row_idx = 0; local_row_idx < 3; ++local_row_idx) {
            size_t glo_row_idx = mesh->triangle_elements[idx_elem].node_inds[local_row_idx];

            for(size_t local_col_idx = 0; local_col_idx < 3; ++local_col_idx) {
                size_t glo_col_idx = mesh->triangle_elements[idx_elem].node_inds[local_col_idx];
                // printf("FUCK YOU!\n");
                global_stiffness_matrix[glo_row_idx][glo_col_idx] += local_stiffness_matrix[local_row_idx][local_col_idx];
            }
        }
    }
    for (size_t i = 0; i < num_nodes; ++i) {
        for (size_t j = 0; j < num_nodes; ++j) {
            printf(" %f", global_stiffness_matrix[i][j]);
        } printf("\n");
    }
    return global_stiffness_matrix;
}

#endif