#ifndef GRID_HPP
#define GRID_HPP

#include "vertex.hpp"
#include "edge.hpp"
#include "cell.hpp"

// #include "./aerodynamics_part/inviscid_flow_field/base/cell.hpp"

#include <vector>

class grid {
private:
    std::vector<vertex> vertices;
    std::vector<edge> edges;
    std::vector<cell> cells;

private:
    std::vector<std::array<double, 2>> points;

public:
    grid() = default;

    void add_vertex(unsigned int id, bool is_dirichlet = false, bool is_neumann = false) {
        vertices.emplace_back(id, is_dirichlet, is_neumann);
    }
    void add_edge(unsigned int id, unsigned int vertex1, unsigned int vertex2, bool is_neumann = false) {
        edges.emplace_back(id, vertex1, vertex2, is_neumann);
    }
    void add_cell(unsigned int v1, unsigned int v2, unsigned int v3) {
        cells.emplace_back(v1, v2, v3);
    }

    unsigned int get_num_vertices() const { return vertices.size(); }
    unsigned int get_num_edges() const { return edges.size(); }
    unsigned int get_num_cells() const { return cells.size(); }

    const std::vector<vertex>& get_vertices() const { return vertices; }
    const std::vector<edge>& get_edges() const { return edges; }
    const std::vector<cell>& get_cells() const { return cells; }
};

#endif
