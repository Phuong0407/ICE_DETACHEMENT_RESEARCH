#ifndef edge_hpp
#define edge_hpp

#include <array>
#include <exception>
#include <stdexcept>
#include <limits>

class edge {
private:
    unsigned int id;
    std::array<unsigned int, 2> cell_ids;
    std::array<unsigned int, 2> vertex_ids;
    bool is_neumann = false;
    bool is_boundary = false;

public:
    edge() = default;
    edge(unsigned int edge_id, unsigned int cell_id_1, unsigned int cell_id_2, unsigned int vertex_id_1, unsigned int vertex_id_2, bool neumann = false) : vertex_ids{vertex_id_1, vertex_id_2}, is_neumann(neumann) {
        if (is_neumann) {
            cell_ids[0] = cell_id_1;
            cell_ids[1] = -1;
        } else
            cell_ids = {cell_id_1, cell_id_2};
    }
    void init_id(unsigned int id) { this->id = id; }
    void init_cell_ids(std::array<unsigned int, 2> cell_ids) { this->cell_ids = std::move(cell_ids); }
    void init_cell_ids(unsigned int c1, unsigned int c2) { this->cell_ids = std::array<unsigned int, 2>{c1, c2}; }
    void init_cell_ids(std::vector<unsigned int> cell_ids) {
        const unsigned int num_cells = cell_ids.size();
        if (num_cells != 1 && num_cells != 2) {
            std::runtime_error("There is more than two cells connected to the edges. The mesh is abnormal.");
        } else if (num_cells == 1) {
            is_boundary = true;
            this->cell_ids[0] = cell_ids[0];
            this->cell_ids[1] = std::numeric_limits<unsigned int>::max();
        } else {
            this->cell_ids[0] = cell_ids[0];
            this->cell_ids[1] = cell_ids[1];
        }
    }
    void init_vertex_ids(std::array<unsigned int, 2> vertex_ids) { this->vertex_ids = std::move(vertex_ids); }
    void init_vertex_ids(unsigned int v1, unsigned int v2) { this->vertex_ids = std::array<unsigned int, 2>{v1, v2}; }

public:
    const std::array<unsigned int, 2>& get_cell_ids() const { return cell_ids; }
    const std::array<unsigned int, 2>& get_vertices() const { return vertex_ids; }
    unsigned int get_cell_ids(unsigned int idx) const { return cell_ids[idx]; }
    unsigned int get_vertex_id(unsigned int idx) const { return vertex_ids[idx]; }
    unsigned int get_id() const { return this->id; }
    bool is_neumann_bound() const { return is_neumann; }
    bool is_internal() const { return !is_boundary; }
    void print() const {
        std::cout << "edge id = " << id << " (" << vertex_ids[0] << "<--->" << vertex_ids[1] << ")";
        std::cout << " ";
        std::cout << "in cell " << cell_ids[0] << " ";
        if (cell_ids[1] != std::numeric_limits<unsigned int>::max())
            std::cout << cell_ids[1];
        else
            std::cout << "(boundary edge)";
        std::cout << "   ";
    }
};

#endif
