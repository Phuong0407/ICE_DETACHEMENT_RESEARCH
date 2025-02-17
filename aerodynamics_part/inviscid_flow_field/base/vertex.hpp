#ifndef vertex_hpp
#define vertex_hpp

#include <vector>
#include <array>
#include <set>

class vertex{
private:
    unsigned int id;
    std::vector<unsigned int> cell_ids;
    std::vector<unsigned int> neighbor_ids;
    bool is_dirichlet = false;
    bool is_neumann = false;

public:
    vertex() = default;
    void init_vertex_id(unsigned int id) { this->id = id; }
    void init_cell_ids(std::vector<unsigned int> cell_ids) { this->cell_ids = std::move(cell_ids); }
    void init_cell_ids(unsigned int cell_id) { cell_ids.push_back(cell_id); }
    void init_neighbor_ids(std::vector<unsigned int> neighbor_ids) { this->neighbor_ids = std::move(neighbor_ids); }
    void init_neighbor_ids(const std::set<unsigned int>& neighbors) { neighbor_ids.assign(neighbors.begin(), neighbors.end()); }

public:
    unsigned int get_id() const { return id; }
    bool is_bound() const { return is_dirichlet || is_neumann; }
    bool is_dirichlet_bound() const { return is_dirichlet; }
    bool is_neumann_bound() const { return is_neumann; }
    std::vector<unsigned int> get_cell_ids() const { return cell_ids; }
    std::vector<unsigned int> get_neighbor_ids() const { return neighbor_ids; }
};

#endif