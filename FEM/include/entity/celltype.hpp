#ifndef CELLTYPE_HPP
#define CELLTYPE_HPP

#include <string>

namespace mesh_entity {

    enum class celltype_t {
        undefined,
        triangle,
        quadrilateral,
        tetrahedron,
        hexahedron
    };

    class celltype_helper {
    public:
        static std::string to_string(celltype_t type) {
            if (type == celltype_t::triangle)
                return "triangle";
            else if (type == celltype_t::quadrilateral)
                return "quadrilateral";
            else if (type == celltype_t::tetrahedron)
                return "tetrahedron";
            else if (type == celltype_t::hexahedron)
                return "hexahedron";
            else
                return "undefined";
        }
    };

} // end namepsace mesh_entity

#endif