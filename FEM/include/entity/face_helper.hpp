#ifndef FACE_HELPER_HPP
#define FACE_HELPER_HPP

#include "facetype.hpp"
#include "faceorder.hpp"

namespace mesh_entity
{
    class face_helper {
        public:
            static constexpr unsigned int __nse(facetype_t type) {
                if (type == facetype_t::triangle)
                    return 3;
                else if (type == facetype_t::quadrilateral)
                    return 4;
                else
                    return 0;
            }
            
            static constexpr bool __valid_celltype(unsigned int spdim, facetype_t type) {
                return ((spdim == 3 && (type == facetype_t::triangle || type == facetype_t::quadrilateral)));
            }

            static constexpr unsigned int __nv(unsigned int spdim, facetype_t type, faceorder_t order) {
                unsigned int _nth = faceorder_helper::nth(order);
                if (type == facetype_t::triangle)
                    return (_nth + 1) * (_nth + 2) / 2;
                else if (type == facetype_t::quadrilateral)
                    return (_nth + 1) * (_nth + 1);
                else
                    return 0;
            }
    };

} // namespace mesh_entity


#endif