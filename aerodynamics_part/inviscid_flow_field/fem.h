#ifndef AXISYMMETRIC_FINITE_ELEMENT_H
#define AXISYMMETRIC_FINITE_ELEMENT_H

#include <stdio.h>
#include <stdlib.h>
#include <stddef.h>
#include <math.h>

void safe_free(void **ptr) {
    if (*ptr) {
        free(*ptr);
        *ptr = NULL;
    }
}

typedef struct {
    size_t id;
    double x, y;
} node;

typedef struct {
    size_t id;
    size_t node_inds[3];
} triangle_element;

typedef struct {
    size_t id;
    size_t node_inds[2];
    double flux;
} edge;

typedef struct {
    size_t num_nodes;
    size_t num_elements;
    size_t num_elements_hor, num_elements_ver;
    node *nodes;
    triangle_element *triangle_elements;
    edge *neuman_bound;
    size_t *dirichlet_inds;
} mesh2D;





mesh2D *init_mesh2D(size_t num_elements_hor, size_t num_elements_ver);
mesh2D* generate_algebraic_grid(double L1, double L2, double H1, double H2, double H3, size_t N, size_t M1, size_t M2, size_t M3);
void generate_grid_connection(mesh2D *mesh);
void generate_elliptic_mesh(mesh2D *mesh);
void generate_boundary_node(double L1, double L2, double H1, double H2, double H3, size_t M1, size_t M2, size_t M3, size_t N, node **bottom, node **top, node **left, node **right);
mesh2D* generate_stretching_grid(double alpha, double eta1, node *bottom, node *top, size_t M, size_t N);
mesh2D* generation_by_transfinite_interpolation(size_t N, size_t M, node *bottom, node *top, node *left, node *right);
void visualize_mesh_via_gnu_plot(mesh2D *mesh);





mesh2D *init_mesh2D(size_t num_elements_hor, size_t num_elements_ver) {
    mesh2D *mesh = (mesh2D*)malloc(sizeof(mesh2D));
    mesh->num_elements_hor = num_elements_hor;
    mesh->num_elements_ver = num_elements_ver;
    mesh->num_nodes = (num_elements_hor + 1) * (num_elements_ver + 1);
    mesh->num_elements = 2 * num_elements_hor * num_elements_ver;
    mesh->nodes = (node*)malloc(mesh->num_nodes * sizeof(node));
    mesh->triangle_elements = (triangle_element*)malloc(mesh->num_elements * sizeof(triangle_element));
    return mesh;
}



mesh2D* generate_algebraic_grid(double L1, double L2, double H1, double H2, double H3, size_t N, size_t M1, size_t M2, size_t M3) {
    size_t M = M1 + M2 + M3;
    size_t num_nodes = (M + 1) * (N + 1);
    mesh2D *mesh = init_mesh2D(M, N);

    double dy1 = L1/N, dy2 = L2/N;
    double dx1 = H1/M1, dx2 = H2/M2, dx3 = H3/M3;

    for (size_t i = 0; i <= N; ++i) {
        for (size_t j = 0; j <= M1; ++j) {
            mesh->nodes[i * (M + 1) + j].id = i * (M + 1) + j;
            mesh->nodes[i * (M + 1) + j].x = dx1 * j;
            mesh->nodes[i * (M + 1) + j].y = dy1 * i;
        }

        double a = ((dy2 - dy1) * i + L1 - L2) / H2;
        double b = dy1 * i - a * H1;
        for (size_t j = M1 + 1; j <= M1 + M2; ++j) {
            double x_1 = H1 + dx2 * (j - M1);
            double y_1 = a * x_1 + b;
            mesh->nodes[i * (M + 1) + j].id = i * (M + 1) + j;
            mesh->nodes[i * (M + 1) + j].x = x_1;
            mesh->nodes[i * (M + 1) + j].y = y_1;
        }

        for (size_t j = M1 + M2 + 1; j <= M; ++j) {
            mesh->nodes[i * (M + 1) + j].id = i * (M + 1) + j;
            mesh->nodes[i * (M + 1) + j].x = H1 + H2 + dx3 * (j - M1 - M2);
            mesh->nodes[i * (M + 1) + j].y = L1 - L2 + dy2 * i;
        }
    }
    return mesh;
}



void generate_grid_connection(mesh2D *mesh) {
    size_t N = mesh->num_elements_ver;
    size_t M = mesh->num_elements_hor;
    size_t elem_id = 0;
    for (size_t i = 0; i < N; ++i) {
        for (size_t j = 0; j < M; ++j) {
            size_t n1 = i * (M + 1) + j;
            size_t n2 = n1 + 1;
            size_t n3 = n2 + M + 1;
            size_t n4 = n1 + M + 1;

            mesh->triangle_elements[elem_id].id = elem_id;
            mesh->triangle_elements[elem_id].node_inds[0] = n1;
            mesh->triangle_elements[elem_id].node_inds[1] = n2;
            mesh->triangle_elements[elem_id].node_inds[2] = n3;
            elem_id++;

            mesh->triangle_elements[elem_id].id = elem_id;
            mesh->triangle_elements[elem_id].node_inds[0] = n1;
            mesh->triangle_elements[elem_id].node_inds[1] = n3;
            mesh->triangle_elements[elem_id].node_inds[2] = n4;
            elem_id++;
        }
    }
}



void generate_elliptic_mesh(mesh2D *mesh) {
    size_t num_nodes = mesh->num_nodes;
    double *new_x = (double*)calloc(num_nodes, sizeof(double));
    double *new_y = (double*)calloc(num_nodes, sizeof(double));
    double *alpha = (double*)calloc(num_nodes, sizeof(double));
    double *beta = (double*)calloc(num_nodes, sizeof(double));
    double *gamma = (double*)calloc(num_nodes, sizeof(double));
    
    double min_err = 1e-6;
    size_t max_iter = 10000;
    double *error_1 = (double*)calloc(max_iter, sizeof(double));
    double *error_2 = (double*)calloc(max_iter, sizeof(double));

    size_t N = mesh->num_elements_ver, M = mesh->num_elements_hor;
    for (size_t iter = 0; iter < max_iter; ++iter) {
        for (size_t i = 1; i < N; ++i) {
            for (size_t j = 1; j < M; ++j) {
                size_t ij = i * (M + 1) + j;
                size_t ijp = ij +1 ;
                size_t ijm = ij - 1;
                size_t ipj = ij + M + 1;
                size_t imj = ij - (M + 1);
                size_t ipjp = ipj + 1;
                size_t ipjm = ipj - 1;
                size_t imjp = imj + 1;
                size_t imjm = imj - 1;

                double dxi_x = (mesh->nodes[ijp].x - mesh->nodes[ijm].x);
                double dxi_y = (mesh->nodes[ijp].y - mesh->nodes[ijm].y);
                double deta_x = (mesh->nodes[ipj].x - mesh->nodes[imj].x);
                double deta_y = (mesh->nodes[ipj].y - mesh->nodes[imj].y);

                alpha[ij] = 0.25 * (pow(deta_x, 2.0) + pow(deta_y, 2.0));
                gamma[ij] = 0.25 * (pow(dxi_x, 2.0) + pow(dxi_y, 2.0));
                beta[ij] = 0.0625 * (deta_x * dxi_x + deta_y * dxi_y);
                double denominator = alpha[ij] + gamma[ij] + 1e-9;
                new_x[ij] = 0.5 / denominator * 
                            (-2.0 * beta[ij] * (mesh->nodes[ipjp].x - mesh->nodes[ipjm].x - mesh->nodes[imjp].x + mesh->nodes[imjp].x) +
                                    alpha[ij] * (mesh->nodes[ijp].x + mesh->nodes[ijm].x) +
                                    gamma[ij] * (mesh->nodes[ipj].x + mesh->nodes[imj].x));
                new_y[ij] = 0.5 / denominator * 
                            (-2.0 * beta[ij] * (mesh->nodes[ipjp].y - mesh->nodes[ipjm].y - mesh->nodes[imjp].y + mesh->nodes[imjm].y) +
                                    alpha[ij] * (mesh->nodes[ijp].y + mesh->nodes[ijm].y) +
                                    gamma[ij] * (mesh->nodes[ipj].y + mesh->nodes[imj].y));
            }
        }

        double max_error_x = 0.0, max_error_y = 0.0;
        for (size_t i = 1; i < N; ++i) {
            for (size_t j = 1; j < M; ++j) {
                size_t ij = i * (M + 1) + j;
                max_error_x = fmax(max_error_x, fabs(new_x[i] - mesh->nodes[ij].x));
                max_error_y = fmax(max_error_y, fabs(new_y[ij] - mesh->nodes[ij].y));
            }
        }
        error_1[iter] = max_error_x;
        error_2[iter] = max_error_y;
        for (size_t i = 1; i < N; ++i) {
            for (size_t j = 1; j < M; ++j) {
                size_t ij = i * (M + 1) + j;
                mesh->nodes[ij].x = new_x[ij];
                mesh->nodes[ij].y = new_y[ij];
            }
        }
        if (error_1[iter] < min_err && error_2[iter] < min_err)
            break;
    }
}



void generate_boundary_node(double L1, double L2, double H1, double H2, double H3,
                            size_t M1, size_t M2, size_t M3, size_t N, 
                            node **bottom, node **top, node **left, node **right) {
    size_t M = M1 + M2 + M3;

    double dy1 = L1 / N, dy2 = L2 / N;
    double dx1 = H1 / M1, dx2 = H2 / M2, dx3 = H3 / M3;

    *bottom = (node*)calloc(M + 1, sizeof(node));
    *top = (node*)calloc(M + 1, sizeof(node));
    *left = (node*)calloc(N + 1, sizeof(node));
    *right = (node*)calloc(N + 1, sizeof(node));

    for (size_t i = 0; i <= M1; ++i) {
        double x = H1 / M1 * i;
        (*bottom)[i].x = x;
        (*bottom)[i].y = 0.0;
        (*top)[i].x = x;
        (*top)[i].y = L1;
    }
    double a = (L1 - L2) / H2;
    for (size_t i = M1 + 1; i <= M1 + M2; ++i) {
        double x = H1 + dx2 * (i - M1);
        double y = a * dx2 * (i - M1);
        (*bottom)[i].x = x;
        (*bottom)[i].y = y;
        (*top)[i].x = x;
        (*top)[i].y = L1;
    }
    for (size_t i = M1 + M2 + 1; i <= M; ++i) {
        double x = H1 + H2 + dx3 * (i - M1 - M2);
        (*bottom)[i].x = x;
        (*bottom)[i].y = L1 - L2;
        (*top)[i].x = x;
        (*top)[i].y = L1;
    }
    for (size_t i = 0; i <= N; ++i) {
        (*left)[i].x = 0.0;
        (*right)[i].x = H1 + H2 + H3;
        (*left)[i].y = L1 / N * i;
        (*right)[i].y = L1 - L2 + L2 / N * i;
    }
}



mesh2D* generate_stretching_grid(double alpha, double eta1, node *bottom, node *top, size_t M, size_t N) {
    mesh2D *mesh = init_mesh2D(M, N);
    double dxi = 1.0 / M, deta = 1.0 / N;
    for (size_t i = 0; i <= N; ++i) {
        double eta = i * deta;
        for (size_t j = 0; j <= M; ++j) {
            size_t ij = i * (M + 1) + j;
            double xi = j * dxi;
            mesh->nodes[ij].id = ij;
            mesh->nodes[ij].x = bottom[j].x;
            if (eta <= eta1) {
                mesh->nodes[ij].y = (top[j].y - bottom[j].y) * eta1 * (exp(alpha * eta / eta1) - 1.0) / (exp(alpha) - 1) + bottom[j].y;
            }
            else {
                mesh->nodes[ij].y = (top[j].y - bottom[j].y) * (1.0 - (1.0 - eta1) * (exp(alpha * (1 - eta) / (1 - eta1)) - 1.0) / (exp(alpha) - 1)) + bottom[j].y;
            }
        }
    }
    return mesh;
}



mesh2D* generation_by_transfinite_interpolation(size_t N, size_t M, node *bottom, node *top, node *left, node *right) {
    mesh2D *mesh = init_mesh2D(M, N);
    double dxi = 1.0 / M;
    double deta = 1.0 / N;
    for (size_t i = 0; i <= N; ++i) {
        double eta = deta * i;
        for (size_t j = 0; j <= M; ++j) {
            double xi = dxi * j;
            size_t ij = i * (M + 1) + j;
            mesh->nodes[ij].id = ij;
            mesh->nodes[ij].x = (1 - xi) * left[i].x + xi * right[i].x
                     + (1 - eta) * bottom[j].x + eta * top[j].x
                     - (1- xi) * (1 - eta) * bottom[0].x - (1 - xi) * eta * top[0].x
                     - (1 - eta) * xi * bottom[M].x - xi * eta * top[M].x;
            mesh->nodes[ij].y = (1 - xi) * left[i].y + xi * right[i].y +
                      (1 - eta) * bottom[j].y + eta * top[j].y -
                      (1- xi) * (1 - eta) * bottom[0].y - (1 - xi) * eta * top[0].y -
                      (1 - eta) * xi * bottom[M].y - xi * eta * top[M].y;
        }        
    }
    return mesh;
}



void visualize_mesh_via_gnu_plot(mesh2D *mesh) {
    FILE *nodeFile = fopen("grid_points.dat", "w");
    FILE *connFile = fopen("grid_connections.dat", "w");

    if (!nodeFile || !connFile) {
        fprintf(stderr, "Error: Unable to open output files!\n");
        if (nodeFile) fclose(nodeFile);
        if (connFile) fclose(connFile);
        return;
    }

    for (size_t i = 0; i < mesh->num_nodes; i++) {
        fprintf(nodeFile, "%lf %lf\n", mesh->nodes[i].x, mesh->nodes[i].y);
    }
    fclose(nodeFile);

    for (size_t i = 0; i < mesh->num_elements; i++) {
        size_t idx1 = mesh->triangle_elements[i].node_inds[0];
        size_t idx2 = mesh->triangle_elements[i].node_inds[1];
        size_t idx3 = mesh->triangle_elements[i].node_inds[2];

        fprintf(connFile, "%lf %lf\n%lf %lf\n\n", mesh->nodes[idx1].x, mesh->nodes[idx1].y,
                mesh->nodes[idx2].x, mesh->nodes[idx2].y);
        fprintf(connFile, "%lf %lf\n%lf %lf\n\n", mesh->nodes[idx2].x, mesh->nodes[idx2].y,
                mesh->nodes[idx3].x, mesh->nodes[idx3].y);
        fprintf(connFile, "%lf %lf\n%lf %lf\n\n", mesh->nodes[idx3].x, mesh->nodes[idx3].y,
                mesh->nodes[idx1].x, mesh->nodes[idx1].y);
    }
    fclose(connFile);

    FILE *scriptFile = fopen("plot_grid.gnu", "w");
    if (!scriptFile) {
        fprintf(stderr, "Error: Unable to create Gnuplot script!\n");
        return;
    }
    fprintf(scriptFile, "set title 'Computational Mesh'\n");
    fprintf(scriptFile, "set xlabel 'X-axis'\n");
    fprintf(scriptFile, "set ylabel 'Y-axis'\n");
    fprintf(scriptFile, "set grid\n");
    fprintf(scriptFile, "set key off\n");
    fprintf(scriptFile, "plot 'grid_points.dat' with points pt 7 ps 1 lc rgb 'blue', \\\n");
    fprintf(scriptFile, "     'grid_connections.dat' with lines lc rgb 'black'\n");
    fclose(scriptFile);

    system("gnuplot -p plot_grid.gnu");
}



#endif