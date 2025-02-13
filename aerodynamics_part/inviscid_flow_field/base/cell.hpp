#ifndef cell_hpp
#define cell_hpp

#include <array>
#include <vector>
#include <optional>
#include <cassert>

class cell {
private:
    std::array<unsigned int, 3> vertex_ids;

    /**
     * Adaptive mesh refinement attributes
     */
    std::optional<unsigned int> parent_id;
    std::vector<unsigned int> children_ids;
    bool is_active = false;

public:
    cell() = default;
    cell(unsigned int v1, unsigned int v2, unsigned int v3) : vertex_ids{v1, v2, v3} {}

    const std::array<unsigned int, 3>& get_nodes() const { return vertex_ids; }
    // void add_edge(unsigned int edge_id) { edge_ids.push_back(edge_id);}
    // const std::vector<unsigned int>& get_edges() const { return edge_ids; }

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
};

#endif