#ifndef DATA_STRUCTURES_FEM_H
#define DATA_STRUCTURES_FEM_H

#include <stddef.h>

typedef struct {
    double x, y;
} node;

typedef struct {
    size_t node_inds[3];
} triangle_element;

typedef struct {
    size_t node_inds[2];
    size_t element_indx;
    double flux;
} edge;

typedef struct {
    size_t num_nodes;
    size_t num_elements;
    size_t num_elements_hor, num_elements_ver;
    size_t num_neumann_edges, num_dirichlet_nodes;
    node *nodes;
    triangle_element *triangle_elements;
    edge *neuman_bound;
    size_t *dirichlet_inds;
    double *dirichlet_bound;
} mesh2D;

typedef struct {
    size_t ind;
    double dev_x;
    double dev_y;
} gradient2D;

#endif