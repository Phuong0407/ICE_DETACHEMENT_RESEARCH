#ifndef FACETYPE_HPP
#define FACETYPE_HPP

#include <string>

namespace mesh_entity
{
    // face is alway 3D
    enum class facetype_t {
        triangle,
        quadrilateral
    };

    class facetype_helper {
        public:
            static std::string to_string(facetype_t type) {
                if (type == facetype_t::triangle)
                    return "triangle";
                else if (type == facetype_t::quadrilateral)
                    return "quadrilateral";
                else
                    return "undefined";
        }
    };
    
} // namespace mesh_entity


#endif