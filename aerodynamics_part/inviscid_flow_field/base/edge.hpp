#ifndef edge_hpp
#define edge_hpp

#include <array>

class edge {
private:
    std::array<unsigned int, 2> cell_ids;
    std::array<unsigned int, 2> vertex_ids;
    bool is_neumann;

public:
    edge() = default;
    edge(unsigned int edge_id, unsigned int cell_id_1, unsigned int cell_id_2, unsigned int vertex_id_1, unsigned int vertex_id_2, bool neumann = false) : vertex_ids{vertex_id_1, vertex_id_2}, is_neumann(neumann) {
        if (is_neumann) {
            cell_ids[0] = cell_id_1;
            cell_ids[1] = -1;
        } else
            cell_ids = {cell_id_1, cell_id_2};
    }

public:
    const std::array<unsigned int, 2>& get_cell_ids() const { return cell_ids; }
    const std::array<unsigned int, 2>& get_vertices() const { return vertex_ids; }
    bool is_neumann_bound() const { return is_neumann; }
    bool is_internal() const { return !is_neumann; }
};

#endif
