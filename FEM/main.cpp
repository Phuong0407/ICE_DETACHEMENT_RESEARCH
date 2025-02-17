#include "./include/entity/vertex.hpp"
#include "./include/entity/edge.hpp"
#include "./include/entity/face.hpp"
#include "./include/entity/cell.hpp"

#define index_t unsigned int

int main() {
    mesh_entity::vertex<2, index_t> v(200);
    v.print();
    std::cout << v.nv() << std::endl;

    mesh_entity::edge<2, index_t> e(4);
    e.print();
    std::cout << e.nv() << std::endl;
    
    mesh_entity::face<3, index_t> f(4, mesh_entity::facetype_t::triangle, mesh_entity::faceorder_t::linear);
    f.print();
    std::cout << f.nv() << std::endl;
    
    mesh_entity::cell<3, index_t> c(5, mesh_entity::celltype_t::tetrahedron, mesh_entity::cellorder_t::quintic);
    c.print();
    std::cout << c.nv() << std::endl;
}