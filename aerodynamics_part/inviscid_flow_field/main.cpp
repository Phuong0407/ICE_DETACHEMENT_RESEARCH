#include "geometry.h"
#include "axisymmetric_fem.h"
#include "C_FORTRAN_caller.h"
#include "gradient_recovery.hpp"

#include <time.h>

#include <iostream>

int main() {

    // clock_t start, end;
    // double cpu_time_used;

    // start = clock();

    mesh2D *mesh = generate_algebraic_grid(4.0, 2.0, 3.0, 3.0, 3.0, 2, 2, 2, 2);
    // generate_elliptic_mesh(mesh);
    generate_grid_connection(mesh);
    generate_neumann_boundary_conditions(0.0, 4.0/3.0, 0.0, -1.0, mesh);
    // visualize_mesh_via_gnu_plot(mesh);

    // double *global_stiffness_matrix = compute_global_stiffness_matrix(mesh);
    // double *force_matrix = (double*)calloc(mesh->num_nodes, sizeof(double));
    // apply_neumann_boundary_conditions(mesh, force_matrix);

    // int n = mesh->num_nodes;
    // double *solution = (double*)calloc(n, sizeof(double));
    // init_sparse_solver(global_stiffness_matrix, &n);
    // solve_sparse_system(force_matrix, solution, &n);

    // generate_node_recovery_patch_connection(mesh->triangle_elements, mesh->num_elements, mesh->num_nodes);

    // std::vector<std::unordered_set<std::size_t>> first_level_patch;
    // std::vector<std::unordered_set<std::size_t>> a = generate_full_level_patch(mesh->triangle_elements, mesh->num_elements, mesh->num_nodes, first_level_patch);

    // for (const auto & b : a) {
    //     for (const auto & c : b) {
    //         std::cout << c << " ";
    //     } std::cout << std::endl;
    // }







    // end = clock();
    // cpu_time_used = ((double)(end - start)) / CLOCKS_PER_SEC;
    // printf("Execution time: %f seconds\n", cpu_time_used);

    return 0;
}