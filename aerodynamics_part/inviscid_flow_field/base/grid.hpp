#ifndef grid_hpp
#define grid_hpp

#include "vertex.hpp"
#include "edge.hpp"
#include "cell.hpp"

#include <set>
#include <unordered_map>
#include <vector>
#include <iostream>
#include <exception>
#include <limits>

struct ArrayHash {
    std::size_t operator()(const std::array<unsigned int, 2>& arr) const {
        return std::hash<unsigned int>{}(arr[0]) ^ (std::hash<unsigned int>{}(arr[1]) << 1);
    }
};

class grid {
private:
    unsigned int num_nodes;
    unsigned int num_cells;
    unsigned int num_edges;

    std::vector<vertex> vertices;
    std::vector<edge> edges;
    std::vector<cell> cells;

private:
    std::vector<double> cell_areas;
    std::vector<double> x, y;

public:
    grid() = default;

    void init_grid_coordinate(std::vector<double> &x, std::vector<double> &y) {
        if (x.size() != y.size())
            throw std::runtime_error("Mismatch between number of x coordinates and number of y cooridinates. The program terminates now.");

        this->num_nodes = x.size();
        this->x = std::move(x);
        this->y = std::move(y);
    }

    void init_cells(const std::vector<std::array<unsigned int, 3>> &connection) {
        this->num_cells = connection.size();
        cells.clear();
        cells.resize(num_cells);
        for (unsigned int i = 0; i < num_cells; ++i) {
            cells[i].init_id(i);
            cells[i].init_vertices(connection[i]);
        }
    }

    void init_edges() {    
        std::unordered_map<std::array<unsigned int, 2>, unsigned int, ArrayHash> edges_map;
        std::unordered_map<std::array<unsigned int, 2>, std::vector<unsigned int>, ArrayHash> edge_to_cell;
    
        unsigned int edge_id = 0;
    
        for (unsigned int cell_id = 0; cell_id < cells.size(); ++cell_id) {
            const auto& elem = cells[cell_id];
    
            std::array<unsigned int, 2> cell_edges[3] = {
                {std::min(elem.get_node_id(0), elem.get_node_id(1)), std::max(elem.get_node_id(0), elem.get_node_id(1))},
                {std::min(elem.get_node_id(1), elem.get_node_id(2)), std::max(elem.get_node_id(1), elem.get_node_id(2))},
                {std::min(elem.get_node_id(2), elem.get_node_id(0)), std::max(elem.get_node_id(2), elem.get_node_id(0))}
            };

            for (const auto& edge : cell_edges) {
                edge_to_cell[edge].push_back(cell_id);
                auto [it, inserted] = edges_map.try_emplace(edge, edge_id);
                if (inserted) {
                    ++edge_id;
                }
            }
        }
    
        num_edges = edges_map.size();
        this->edges.resize(num_edges);
    
        for (const auto& [edge_dat, id] : edges_map) {
            this->edges[id].init_id(id);
            this->edges[id].init_vertex_ids(edge_dat);
            this->edges[id].init_cell_ids(edge_to_cell[edge_dat]);
        }
        // for (const auto& cell_ : cells) {
        //     std::cout << "[" << cell_.get_edge_id(0) << "(" << edges[cell_.get_edge_id(0)].get_vertex_id(0) << "<--->" << edges[cell_.get_edge_id(0)].get_vertex_id(1) << "), "
        //               << cell_.get_edge_id(1) << "(" << edges[cell_.get_edge_id(1)].get_vertex_id(0) << "<--->" << edges[cell_.get_edge_id(1)].get_vertex_id(1) << "), " 
        //               << cell_.get_edge_id(2) << "(" << edges[cell_.get_edge_id(2)].get_vertex_id(0) << "<--->" << edges[cell_.get_edge_id(2)].get_vertex_id(1) << ")]" 
        //               << "---"
        //               << "[" << cell_.get_node_id(0) << ", " << cell_.get_node_id(1) << ", " << cell_.get_node_id(2) << "]" << std::endl;
        // }
    }

    void init_cells_to_edges() {
        std::vector<std::vector<unsigned int>> cell_ids_(cells.size());
    
        for (const auto& edge_ : edges) {
            const auto& cell_ids = edge_.get_cell_ids();
            for (const auto& id : cell_ids) {
                if (id != std::numeric_limits<unsigned int>::max()) {
                    cell_ids_[id].push_back(edge_.get_id());
                }
            }
        }

        for (unsigned int i = 0; i < num_cells; ++i) {
            cells[i].init_edges(cell_ids_[i]);
        }

        // for (size_t i = 0; i < cells.size(); ++i) {
        //     std::cout << "Cell " << i << ": [";
        //     std::cout << cells[i].get_edge_id(0) << "(" << edges[cells[i].get_edge_id(0)].get_vertex_id(0) << "---" << edges[cells[i].get_edge_id(0)].get_vertex_id(1) << ")"
        //               << ", "
        //               << cells[i].get_edge_id(1) << "(" << edges[cells[i].get_edge_id(1)].get_vertex_id(0) << "---" << edges[cells[i].get_edge_id(1)].get_vertex_id(1) << ")"
        //               << ", "
        //               << cells[i].get_edge_id(2) << "(" << edges[cells[i].get_edge_id(2)].get_vertex_id(0) << "---" << edges[cells[i].get_edge_id(2)].get_vertex_id(1) << ")";
        //     std::cout << std::endl;


        //         std::cout << cell_ids_[i][j] << " ";
        //     std::cout << "]---";
        //     for (size_t j = 0; j < cell_ids_[i].size(); ++j) {
        //         std::cout << "[" << edges[cell_ids_[i][j]].get_vertex_id(0) << ", " << edges[cell_ids_[i][j]].get_vertex_id(1) << "]";
        //         // ", " << edges[cell_ids_[i][j]].get_vertex_id(2) << 
        //     }
        //     std::cout << std::endl;
        // }

        for (size_t i = 0; i < cells.size(); ++i) {
            std::cout << "Cell " << i << ":" << std::endl;
            std::cout << "   ";
            edges[cells[i].get_edge_id(0)].print();
            std::cout << std::endl;
            std::cout << "   ";
            edges[cells[i].get_edge_id(1)].print();
            std::cout << std::endl;
            std::cout << "   ";
            edges[cells[i].get_edge_id(2)].print();
            std::cout << std::endl;
        }
    }

    void init_vertices() {
        vertices.clear();
        vertices.resize(num_nodes);
        std::vector<std::set<unsigned int>> vertices_set(num_nodes);

        unsigned int num_cells = cells.size();
        for (auto const &elem : cells) {
            unsigned int idx1 = elem.get_node_id(0);
            unsigned int idx2 = elem.get_node_id(1);
            unsigned int idx3 = elem.get_node_id(2);
            vertices_set[idx1].insert(idx2);
            vertices_set[idx1].insert(idx3);
            vertices_set[idx2].insert(idx3);
            vertices_set[idx2].insert(idx1);
            vertices_set[idx3].insert(idx1);
            vertices_set[idx3].insert(idx2);
            vertices[idx1].init_cell_ids(elem.get_cell_id());
            vertices[idx2].init_cell_ids(elem.get_cell_id());
            vertices[idx3].init_cell_ids(elem.get_cell_id());
        }
        for (unsigned int i = 0; i < num_nodes; ++i) {
            vertices[i].init_vertex_id(i);
            vertices[i].init_neighbor_ids(vertices_set[i]);
        }
    }

    // void init_cells_to_edges() {
    //     for (const auto &edge : edges) {
    //         cells[edge.get_cell_ids(0)].init_edges();
    //     }
    // }

    unsigned int get_num_vertices() const { return vertices.size(); }
    unsigned int get_num_edges() const { return edges.size(); }
    unsigned int get_num_cells() const { return cells.size(); }

    const std::vector<vertex>& get_vertices() const { return vertices; }
    const std::vector<edge>& get_edges() const { return edges; }
    const std::vector<cell>& get_cells() const { return cells; }
};

#endif