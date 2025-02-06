#ifndef GRADIENT_RECOVERY_H
#define GRADIENT_RECOVERY_H

#include "data_structures_fem.h"
#include "dense_solver.h"

#include <stdlib.h>
#include <string.h>
#include <unordered_set>
#include <vector>
#include <iostream>
#include <cmath>

constexpr std::size_t order_of_patch = 6;





inline std::vector<std::unordered_set<std::size_t>> generate_first_level_patch(const triangle_element* element_connection, std::size_t num_elements, std::size_t num_nodes) {
    std::vector<std::unordered_set<std::size_t>> first_level_patch(num_nodes);

    for (std::size_t i = 0; i < num_elements; ++i) {
        std::size_t idx_elem_1 = element_connection[i].node_inds[0];
        std::size_t idx_elem_2 = element_connection[i].node_inds[1];
        std::size_t idx_elem_3 = element_connection[i].node_inds[2];

        first_level_patch[idx_elem_1].insert({idx_elem_2, idx_elem_3});
        first_level_patch[idx_elem_2].insert({idx_elem_3, idx_elem_1});
        first_level_patch[idx_elem_3].insert({idx_elem_1, idx_elem_2});
    }
    return first_level_patch;
}



inline std::unordered_set<std::size_t> find_all_internal_nodes(const std::vector<std::unordered_set<std::size_t>> &first_level_patch, const std::unordered_set<std::size_t> &patch_of_node) {
    std::unordered_set<std::size_t> internal_nodes;
    for (const auto& connected_node : patch_of_node)
        if (first_level_patch[connected_node].size() >= order_of_patch)
            internal_nodes.insert(connected_node);
    return internal_nodes;
}





std::vector<std::unordered_set<std::size_t>> generate_full_level_patch(const triangle_element* element_connection, std::size_t num_elements, std::size_t num_nodes, std::vector<std::unordered_set<std::size_t>> &first_level_patch) {
    first_level_patch.clear();
    first_level_patch = generate_first_level_patch(element_connection, num_elements, num_nodes);
    std::vector<std::unordered_set<std::size_t>> full_level_patch = first_level_patch;

    for (std::size_t i = 0; i < num_nodes; ++i) {
        const auto &node_patch = first_level_patch[i];
        if (node_patch.size() >= order_of_patch) continue;

        std::unordered_set<std::size_t> internal_nodes = find_all_internal_nodes(first_level_patch, node_patch);
        for (const auto &internal_node : internal_nodes) {
            full_level_patch[i].insert(first_level_patch[internal_node].begin(), first_level_patch[internal_node].end());
        } full_level_patch[i].erase(i);
    }

    for (std::size_t i = 0; i < num_nodes; ++i) {
        while (full_level_patch[i].size() < order_of_patch) {
            for (const auto &neighbor : full_level_patch[i]) {
                full_level_patch[i].insert(full_level_patch[neighbor].begin(), full_level_patch[neighbor].end());
            }
        } full_level_patch[i].erase(i);
    }
    return full_level_patch;
}



double* compute_AT_A(const double *A, std::size_t row, std::size_t col) {
    double *ATA = (double*)malloc(col * col * sizeof(double));
    memset(ATA, 0, col * col * sizeof(double));

    for (size_t i = 0; i < col; ++i) {
        for (size_t j = 0; j < col; ++j) {
            for (size_t k = 0; k < row; ++k) {
                ATA[i * col + j] += A[k * col + i] * A[k * col + j];
            }
        }
    }
    return ATA;
}



double* compute_AT_B(const double *A, const double *B, std::size_t row, std::size_t col) {
    double *AT_B = (double*)malloc(col * sizeof(double));
    memset(AT_B, 0, col * sizeof(double));

    for (size_t i = 0; i < col; ++i) {
        for (size_t k = 0; k < row; ++k) {
            AT_B[i] += A[k * col + i] * B[k];
        }
    }
    return AT_B;
}





std::vector<double> compute_longest_length_by_first_level_patch(const node* nodes, const std::vector<std::unordered_set<std::size_t>>& first_level_patchs) {
    std::size_t number_nodes = first_level_patchs.size();
    std::vector<double> longest_length(number_nodes, 0.0);

    for (std::size_t i = 0; i < number_nodes; ++i) {
        double max_length = 0.0;
        for (std::size_t j : first_level_patchs[i]) {
            double dx = nodes[i].x - nodes[j].x;
            double dy = nodes[i].y - nodes[j].y;
            double length = std::sqrt(dx * dx + dy * dy);
            if (length > max_length) max_length = length;
        }
        longest_length[i] = max_length;
    }
    return longest_length;
}



gradient2D* generate_recovered_gradient(const mesh2D* mesh, const double *fem_solution) {
    std::size_t num_nodes = mesh->num_nodes;
    std::vector<std::unordered_set<std::size_t>> first_level_patchs;
    std::vector<std::unordered_set<std::size_t>> full_level_patchs = generate_full_level_patch(mesh->triangle_elements, mesh->num_elements, mesh->num_nodes, first_level_patchs);
    std::vector<double> longest_length = compute_longest_length_by_first_level_patch(mesh->nodes, first_level_patchs);
    gradient2D* recovered_gradient = (gradient2D*)malloc(num_nodes * sizeof(gradient2D));

    for (std::size_t i = 0; i < num_nodes; ++i) {
        const auto &full_level_patch = full_level_patchs[i];
        std::size_t size_of_fitting_matrix = full_level_patch.size() + 1;

        double* A = (double*)calloc(size_of_fitting_matrix * order_of_patch, sizeof(double));
        A[0] = 1.0;
        double* local_solution = (double*)calloc(size_of_fitting_matrix, sizeof(double));
        local_solution[0] = fem_solution[i];
        
        double mlength = longest_length[i];
        node node_at_origin = mesh->nodes[i];
        std::size_t k = 1;

        for (std::size_t j : full_level_patch) {
            node patch_node = mesh->nodes[j];
            double xi, eta;
            xi = (patch_node.x - node_at_origin.x) / mlength;
            eta = (patch_node.y - node_at_origin.y) / mlength;
            local_solution[k] = fem_solution[j];
            A[k * order_of_patch] = 1.0;
            A[k * order_of_patch + 1] = xi;
            A[k * order_of_patch + 2] = eta;
            A[k * order_of_patch + 3] = xi * xi;
            A[k * order_of_patch + 4] = xi * eta;
            A[k * order_of_patch + 5] = eta * eta;
            k++;
        }
        double *sqA = compute_AT_A(A, size_of_fitting_matrix, order_of_patch);
        double *B = compute_AT_B(A, local_solution, size_of_fitting_matrix, order_of_patch);
        double *solution = nullptr;
        solve_dense_system(sqA, B, &solution, order_of_patch);
        recovered_gradient[i].ind = i;
        recovered_gradient[i].dev_x = solution[1] / mlength;
        recovered_gradient[i].dev_y = solution[2] / mlength;
        if(solution != NULL) free(solution);
        free(local_solution);
        free(sqA);
        free(B);
        free(A);
    }
    return recovered_gradient;
}

#endif