#ifndef grid_hpp
#define grid_hpp

#include "vertex.hpp"
#include "edge.hpp"
#include "cell.hpp"

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

    unsigned int get_num_vertices() const { return vertices.size(); }
    unsigned int get_num_edges() const { return edges.size(); }
    unsigned int get_num_cells() const { return cells.size(); }

    const std::vector<vertex>& get_vertices() const { return vertices; }
    const std::vector<edge>& get_edges() const { return edges; }
    const std::vector<cell>& get_cells() const { return cells; }
};

#endif