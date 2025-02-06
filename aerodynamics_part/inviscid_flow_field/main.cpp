#include "geometry.h"
#include "axisymmetric_fem.h"
#include "gradient_recovery.hpp"

#include <time.h>
#include <iostream>



void print_gradient2D_to_file(const char *filename, mesh2D* mesh, gradient2D *data, size_t num_nodes) {
    FILE *file = fopen(filename, "w");
    if (!file) {
        perror("Error opening file");
        return;
    }

    for (size_t i = 0; i < num_nodes; ++i) {
        fprintf(file, "gradient(%lf, %lf) = (%.10f, %.10f)\n", mesh->nodes[i].x, mesh->nodes[i].y, data[i].dev_x, data[i].dev_y);
    }

    fclose(file);
}



int main() {

    clock_t start, end;
    double cpu_time_used;
    start = clock();

    mesh2D *mesh = generate_algebraic_grid(4.0, 2.0, 3.0, 3.0, 3.0, 100, 5, 100, 5);
    generate_elliptic_mesh(mesh);
    generate_grid_connection(mesh);
    visualize_mesh_via_gnu_plot(mesh);
    generate_neumann_boundary_conditions(0.0, 4.0/3.0, 0.0, -1.0, mesh);
    int n = mesh->num_nodes;
    double *solution = laplace_solver(mesh);
    gradient2D *recovered_gradient = generate_recovered_gradient(mesh, solution);
    print_gradient2D_to_file("gradient_recovery_result.txt", mesh, recovered_gradient, n);
    free_gradient2D(recovered_gradient);
    free(solution);
    free_mesh2D(mesh);
    end = clock();
    cpu_time_used = ((double)(end - start)) / CLOCKS_PER_SEC;
    printf("Execution time: %f seconds\n", cpu_time_used);

    return 0;
}