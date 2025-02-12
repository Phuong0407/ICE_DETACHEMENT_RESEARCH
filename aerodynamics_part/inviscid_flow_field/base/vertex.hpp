#ifndef VERTEX_HPP
#define VERTEX_HPP

#include <vector>
#include <array>

class vertex{
private:
    unsigned int id;
    std::vector<unsigned int> cell_ids;
    std::vector<unsigned int> neighbor_ids;
    bool is_dirichlet = false;
    bool is_neumann = false;

public:
    vertex() = default;
    vertex(unsigned int id, bool dirichlet = false, bool neumann = false) : id(id), is_dirichlet(is_dirichlet), is_neumann(is_neumann) {}
    void add_cell(unsigned int cell_id) { cell_ids.push_back(cell_id); }
    void add_neighbor(unsigned int neighbor_id) { neighbor_ids.push_back(neighbor_id); }

public:
    unsigned int get_id() const { return id; }
    bool is_bound() const { return is_dirichlet || is_neumann; }
    bool is_dirichlet_bound() const { return is_dirichlet; }
    bool is_neumann_bound() const { return is_neumann; }
    std::vector<unsigned int> get_cell_ids() const { return cell_ids; }
};

#endif