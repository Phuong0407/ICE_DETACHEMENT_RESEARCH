#ifndef DATA_STRUCTURES_FEM_H
#define DATA_STRUCTURES_FEM_H

#include <stddef.h>
#include <stdlib.h>

// BASIC MESH STRUCTURE
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
    double *element_areas;
    edge *neuman_bound;
    size_t *dirichlet_inds;
    double *dirichlet_bound;
} mesh2D;

// Gradient Structure to store gradient data
typedef struct {
    size_t ind;
    double dev_x;
    double dev_y;
} gradient2D;

typedef struct {
    mesh2D *mesh;
    double *fem_solution;
    gradient2D *raw_gradient;
    gradient2D *recovered_gradient;
    double* errors_sq;
    double error_sq_total;
} gradient_refinement2D;


typedef struct {
    size_t node_inds[2];
    size_t left_element;
    size_t right_element;
    size_t mid_node;
} edge2D;

typedef struct {
    size_t key;
    size_t mid_node_index;
} midpoint_map;

size_t edge_key(size_t node1, size_t node2) {
    return (node1 < node2) ? (node1 * 1000000 + node2) : (node2 * 1000000 + node1);
}

typedef struct quadtree_node {
    size_t element_index;
    size_t child_elements[4];
    bool is_refined;
    struct quadtree_node *children[4];
} quadtree_node;

typedef struct quadtree {
    size_t num_nodes;
    quadtree_node *root;
    quadtree_node *nodes;
} quadtree;

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
    if (mesh->element_areas) free(mesh->element_areas);
    if (mesh->neuman_bound) free(mesh->neuman_bound);
    if (mesh->dirichlet_inds) free(mesh->dirichlet_inds);
    if (mesh->dirichlet_bound) free(mesh->dirichlet_bound);

    mesh->nodes = NULL;
    mesh->triangle_elements = NULL;
    mesh->element_areas = NULL;
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