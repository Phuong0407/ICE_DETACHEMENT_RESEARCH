#ifndef cell_hpp
#define cell_hpp

#include <array>
#include <vector>
#include <optional>
#include <cassert>
#include <exception>
#include <algorithm>

class cell {
private:
    unsigned int id;
    std::array<unsigned int, 3> vertex_ids;
    std::array<unsigned int, 3> edge_ids;

    /**
     * Adaptive mesh refinement attributes
     */
    bool is_active = false;
    std::optional<unsigned int> parent_id;
    std::vector<unsigned int> children_ids;

public:
    cell() = default;
    void init_id(unsigned int id) { this->id = id; }
    void init_vertices(unsigned int v1, unsigned int v2, unsigned int v3) {
        vertex_ids = std::array<unsigned int, 3>{v1, v2, v3};
    }
    void init_vertices(std::array<unsigned int, 3> vertex_ids) {
        this->vertex_ids = std::move(vertex_ids);
    }
    void init_edges(unsigned int e1, unsigned int e2, unsigned int e3) {
        edge_ids = std::array<unsigned int, 3>{e1, e2, e3};
    }
    void init_edges(std::array<unsigned int, 3> edge_ids) {
        this->edge_ids = std::move(edge_ids);
    }
    void init_edges(std::vector<unsigned int> edge_ids) {
        if (edge_ids.size() != 3)
            throw std::runtime_error("Mismatch between number of x coordinates and number of y cooridinates. The program terminates now.");
        std::copy_n(edge_ids.begin(), 3, this->edge_ids.begin());
    }
    const std::array<unsigned int, 3>& get_node_ids() const { return vertex_ids; }
    unsigned int get_node_id(unsigned int idx) const { return vertex_ids[idx]; }
    const std::array<unsigned int, 3> get_edge_ids() const { return edge_ids; }
    unsigned int get_edge_id(unsigned int idx) const { return edge_ids[idx]; }
    
    /**
     * Adaptive mesh refinement functions
     */
    void set_parent(unsigned int parent) {
        parent_id = parent;
        is_active = false;
    }

    void add_child(unsigned int child) {
        children_ids.push_back(child);
        is_active = false;
    }

    bool is_refined() const { return !children_ids.empty(); }
    bool is_leaf() const { return is_active; }
    std::optional<unsigned int> get_parent() const { return parent_id; }
    const std::vector<unsigned int>& get_children() const { return children_ids; }
    unsigned int get_cell_id() const { return this->id; }
};

#endif