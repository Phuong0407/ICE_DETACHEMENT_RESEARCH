#include "fem.h"

int main() {

    // mesh2D *A = generate_algebraic_grid(4.0, 2.0, 3.0, 3.0, 3.0, 10, 10, 10, 10);
    // generate_elliptic_mesh(A);
    // mesh2D *A = 
    node *bottom, *top, *left, *right;
    generate_boundary_node(4.0, 2.0, 3.0, 3.0, 3.0, 10, 10, 10, 30, &bottom, &top, &left, &right);

    // mesh2D *A = generate_stretching_grid(0.5, 0.001, bottom, top, 30, 30);
    // generate_elliptic_mesh(A);
    mesh2D *A = generation_by_transfinite_interpolation(30, 30, bottom, top, left, right);
    generate_grid_connection(A);

    visualize_mesh_via_gnu_plot(A);

    return 0;
}