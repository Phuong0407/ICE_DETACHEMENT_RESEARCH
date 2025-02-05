#include "geometry.h"
#include "axisymmetric_fem.h"

int main() {

    // mesh2D *A = generate_algebraic_grid(4.0, 2.0, 3.0, 3.0, 3.0, 10, 10, 10, 10);
    // generate_elliptic_mesh(A);
    
    // node *bottom, *top, *left, *right;
    // generate_boundary_node(4.0, 2.0, 3.0, 3.0, 3.0, 10, 10, 10, 30, &bottom, &top, &left, &right);

    // mesh2D *A = generate_stretching_grid(0.5, 0.001, bottom, top, 30, 30);
    // generate_elliptic_mesh(A);

    // mesh2D *A = generation_by_transfinite_interpolation(30, 30, bottom, top, left, right);
    // generate_grid_connection(A);

    mesh2D *mesh = generate_algebraic_grid(4.0, 2.0, 3.0, 3.0, 3.0, 2, 2, 2, 2);
    generate_grid_connection(mesh);
    visualize_mesh_via_gnu_plot(mesh);

    double **global_stiffness_matrix = compute_global_stiffness_matrix(mesh);

    return 0;
}