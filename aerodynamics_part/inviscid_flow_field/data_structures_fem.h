#ifndef DATA_STRUCTURES_FEM_H
#define DATA_STRUCTURES_FEM_H

#include <stddef.h>
#include <stdlib.h>

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

void free_nodes(node *nodes) {
    if (nodes) {
        free(nodes);
        nodes = NULL;
    }
}

void free_triangle_elements(triangle_element *elements) {
    if (elements) {
        free(elements);
        elements = NULL;
    }
}

void free_edges(edge *edges) {
    if (edges) {
        free(edges);
        edges = NULL;
    }
}

void free_mesh2D(mesh2D *mesh) {
    if (!mesh) return;

    if (mesh->nodes) free(mesh->nodes);
    if (mesh->triangle_elements) free(mesh->triangle_elements);
    if (mesh->neuman_bound) free(mesh->neuman_bound);
    if (mesh->dirichlet_inds) free(mesh->dirichlet_inds);
    if (mesh->dirichlet_bound) free(mesh->dirichlet_bound);

    mesh->nodes = NULL;
    mesh->triangle_elements = NULL;
    mesh->neuman_bound = NULL;
    mesh->dirichlet_inds = NULL;
    mesh->dirichlet_bound = NULL;

    free(mesh);
}

void free_gradient2D(gradient2D *gradients) {
    if (gradients) {
        free(gradients);
    }
}

#endif