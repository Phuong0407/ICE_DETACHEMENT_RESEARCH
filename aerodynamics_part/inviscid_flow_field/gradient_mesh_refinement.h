#ifndef GRADIENT_MESH_REFINEMENT_H
#define GRADIENT_MESH_REFINEMENT_H

#include "data_structures_fem.h"

#include <math.h>

// typedef struct {
//     double xi;
//     double eta;
//     double weight;
// } GaussQuadraturePoint;

// const GaussQuadraturePoint Gauss13PointTriangle[13] = {
//     { 0.333333333333333, 0.333333333333333, -0.14957004446767 },
//     { 0.479308067841923, 0.260345966079038,  0.175615257433204 },
//     { 0.260345966079038, 0.479308067841923,  0.175615257433204 },
//     { 0.260345966079038, 0.260345966079038,  0.175615257433204 },
//     { 0.869739794195568, 0.065130102902216,  0.053347235608839 },
//     { 0.065130102902216, 0.869739794195568,  0.053347235608839 },
//     { 0.065130102902216, 0.065130102902216,  0.053347235608839 },
//     { 0.638444188569809, 0.312865496004875,  0.077113760890257 },
//     { 0.638444188569809, 0.048690315425316,  0.077113760890257 },
//     { 0.312865496004875, 0.638444188569809,  0.077113760890257 },
//     { 0.312865496004875, 0.048690315425316,  0.077113760890257 },
//     { 0.048690315425316, 0.638444188569809,  0.077113760890257 },
//     { 0.048690315425316, 0.312865496004875,  0.077113760890257 }
// };





void generate_elementwise_raw_gradient(gradient_refinement2D *gradient_refinement) {
    const mesh2D *mesh = gradient_refinement->mesh;
    const size_t num_nodes = mesh->num_nodes;
    const size_t num_elements = mesh->num_elements;
    const node *nodes = gradient_refinement->mesh->nodes;
    const triangle_element *elements = mesh->triangle_elements;

    gradient_refinement->raw_gradient = (gradient2D*)malloc(num_elements * sizeof(double));

    for (size_t i = 0; i < num_elements; ++i) {
        size_t idx_1 = elements[i].node_inds[0];
        size_t idx_2 = elements[i].node_inds[1];
        size_t idx_3 = elements[i].node_inds[2];

        double x1 = nodes[idx_1].x, y1 = nodes[idx_1].y;
        double x2 = nodes[idx_2].x, y2 = nodes[idx_2].y;
        double x3 = nodes[idx_3].x, y3 = nodes[idx_3].y;

        double phi1 = gradient_refinement->fem_solution[idx_1];
        double phi2 = gradient_refinement->fem_solution[idx_2];
        double phi3 = gradient_refinement->fem_solution[idx_3];

        double A = mesh->element_areas[i];
        double c1 = y2 - y3, c2 = y3 - y1, c3 = y1 - y3;
        double d1 = x3 - x2, d2 = x1 - x3, d3 = x2 - x1;

        double dphi_dx = (phi1 * c1 + phi2 * c2 + phi3 * c3) / (2 * A);
        double dphi_dy = (phi1 * d1 + phi2 * d2 + phi3 * d3) / (2 * A);
        gradient_refinement->raw_gradient[i].dev_x = dphi_dx;
        gradient_refinement->raw_gradient[i].dev_y = dphi_dy;
    }
}



void compute_relative_error_by_element(gradient_refinement2D *gradient_refinement) {
    size_t num_elements = gradient_refinement->mesh->num_elements;
    gradient_refinement->errors_sq = (double*)malloc(num_elements * sizeof(double));

    for (size_t i = 0; i < num_elements; ++i) {
        double error = 0.0;
        double raw_Dev_x = gradient_refinement->raw_gradient[i].dev_x;
        double raw_Dev_y = gradient_refinement->raw_gradient[i].dev_y;

        const triangle_element element = gradient_refinement->mesh->triangle_elements[i];
        size_t node1 = element.node_inds[0];
        size_t node2 = element.node_inds[1];
        size_t node3 = element.node_inds[2];
    
        double x1 = gradient_refinement->mesh->nodes[node1].x;
        double y1 = gradient_refinement->mesh->nodes[node1].y;
        double x2 = gradient_refinement->mesh->nodes[node2].x;
        double y2 = gradient_refinement->mesh->nodes[node2].y;
        double x3 = gradient_refinement->mesh->nodes[node3].x;
        double y3 = gradient_refinement->mesh->nodes[node3].y;
    
        double A = gradient_refinement->mesh->element_areas[i];
    
        double recovered_dev_x_1 = gradient_refinement->recovered_gradient[node1].dev_x;
        double recovered_dev_y_1 = gradient_refinement->recovered_gradient[node1].dev_y;
        double recovered_dev_x_2 = gradient_refinement->recovered_gradient[node2].dev_x;
        double recovered_dev_y_2 = gradient_refinement->recovered_gradient[node2].dev_y;
        double recovered_dev_x_3 = gradient_refinement->recovered_gradient[node3].dev_x;
        double recovered_dev_y_3 = gradient_refinement->recovered_gradient[node3].dev_y;
    
        double x_mid_12 = 0.5 * (x1 + x2);
        double y_mid_12 = 0.5 * (y1 + y2);
        double x_mid_23 = 0.5 * (x2 + x3);
        double y_mid_23 = 0.5 * (y2 + y3);
        double x_mid_31 = 0.5 * (x3 + x1);
        double y_mid_31 = 0.5 * (y3 + y1);

        double N1_12 = 0.5, N2_12 = 0.5, N3_12 = 0.0;
        double N1_23 = 0.0, N2_23 = 0.5, N3_23 = 0.5;
        double N1_31 = 0.5, N2_31 = 0.0, N3_31 = 0.5;

        double recovered_dev_x_12 = recovered_dev_x_1 * N1_12 + recovered_dev_x_2 * N2_12 + recovered_dev_x_3 * N3_12;
        double recovered_dev_y_12 = recovered_dev_y_1 * N1_12 + recovered_dev_y_2 * N2_12 + recovered_dev_y_3 * N3_12;
        double recovered_dev_x_23 = recovered_dev_x_1 * N1_23 + recovered_dev_x_2 * N2_23 + recovered_dev_x_3 * N3_23;
        double recovered_dev_y_23 = recovered_dev_y_1 * N1_23 + recovered_dev_y_2 * N2_23 + recovered_dev_y_3 * N3_23;
        double recovered_dev_x_31 = recovered_dev_x_1 * N1_31 + recovered_dev_x_2 * N2_31 + recovered_dev_x_3 * N3_31;
        double recovered_dev_y_31 = recovered_dev_y_1 * N1_31 + recovered_dev_y_2 * N2_31 + recovered_dev_y_3 * N3_31;

        error += (recovered_dev_x_12 - raw_Dev_x) * (recovered_dev_x_12 - raw_Dev_x) + (recovered_dev_y_12 - raw_Dev_y) * (recovered_dev_y_12 - raw_Dev_y);
        error += (recovered_dev_x_23 - raw_Dev_x) * (recovered_dev_x_23 - raw_Dev_x) + (recovered_dev_y_23 - raw_Dev_y) * (recovered_dev_y_23 - raw_Dev_y);
        error += (recovered_dev_x_31 - raw_Dev_x) * (recovered_dev_x_31 - raw_Dev_x) + (recovered_dev_y_31 - raw_Dev_y) * (recovered_dev_y_31 - raw_Dev_y);

        gradient_refinement->errors_sq[i] = error * A / 3.0;
    }
    for (size_t i = 0; i < num_elements; ++i) {
        gradient_refinement->error_sq_total += gradient_refinement->errors_sq[i];
    }
}



void refine_mesh(gradient_refinement2D *gradient_refinement) {
    size_t num_elements = gradient_refinement->mesh->num_elements;
    double raw_error = 0.0;
    for (size_t i = 0; i < num_elements; ++i) {
        double raw_Dev_x = gradient_refinement->raw_gradient[i].dev_x;
        double raw_Dev_y = gradient_refinement->raw_gradient[i].dev_y;
        raw_error += raw_Dev_x * raw_Dev_x + raw_Dev_y * raw_Dev_y;
    }
    double etol = sqrt(raw_error + gradient_refinement->error_sq_total/num_elements);

    for (size_t i = 0; i < num_elements; ++i) {
        double relative_error = sqrt(gradient_refinement->errors_sq[i]) / etol;
        if (relative_error <= 1.0) continue;

        
    }
}



#endif