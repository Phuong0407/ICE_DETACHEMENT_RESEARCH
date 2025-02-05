#include "fem.h"

int main() {

    mesh2D *A = generate_algebraic_grid(4.0, 2.0, 3.0, 3.0, 3.0, 10, 10, 10, 10);
    generate_elliptic_mesh(A);
    visualize_mesh_via_gnu_plot(A);

    return 0;
}