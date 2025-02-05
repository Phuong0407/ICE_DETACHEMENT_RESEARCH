#include "geometry.h"
#include "axisymmetric_fem.h"
#include "C_FORTRAN_caller.h"

#include <time.h>

int main() {

    clock_t start, end;
    double cpu_time_used;

    start = clock();

    mesh2D *mesh = generate_algebraic_grid(4.0, 2.0, 3.0, 3.0, 3.0, 100, 100, 100, 100);
    // generate_elliptic_mesh(mesh);
    generate_grid_connection(mesh);
    generate_neumann_boundary_conditions(0.0, 4.0/3.0, 0.0, -1.0, mesh);
    // visualize_mesh_via_gnu_plot(mesh);

    double *global_stiffness_matrix = compute_global_stiffness_matrix(mesh);
    double *force_matrix = (double*)calloc(mesh->num_nodes, sizeof(double));
    apply_neumann_boundary_conditions(mesh, force_matrix);

    int n = mesh->num_nodes;
    double *solution = (double*)calloc(n, sizeof(double));
    init_sparse_solver(global_stiffness_matrix, &n);
    solve_sparse_system(force_matrix, solution, &n);

    // for (int i = 0; i < n; i++) {
    //     printf("force_mat[%d] = %lf\n", i, force_matrix[i]);
    // }

    // printf("COMPUTED FEM SOLUTION:\n");
    // for (int i = 0; i < n; i++) {
    //     printf("phi[%d] = %lf\n", i, solution[i]);
    // }

    return 0;
}