#ifndef CELL_HELPER_HPP
#define CELL_HELPER_HPP

#include "celltype.hpp"
#include "cellorder.hpp"

namespace mesh_entity
{
    class cell_helper 
    {
    public:

        static constexpr unsigned int __nse(unsigned int _spdim, celltype_t type) {
            if (_spdim == 2) {
                if (type == celltype_t::triangle)
                    return 3;
                else if (type == celltype_t::quadrilateral)
                    return 4;
                else
                    return 0;
            } else if (_spdim == 3) {
                if (type == celltype_t::tetrahedron)
                    return 4;
                else if (type == celltype_t::hexahedron)
                    return 6;
                else
                    return 0;
            } else {
                return 0;
            }
        }
        
        static constexpr unsigned int __nv(unsigned int _spdim, celltype_t type, cellorder_t order) {
            unsigned int _nth = cellorder_helper::nth(order);
            if (_spdim == 2) {
                if (type == celltype_t::triangle)
                    return (_nth + 1) * (_nth + 2) / 2;
                else if (type == celltype_t::quadrilateral)
                    return (_nth + 1) * (_nth + 1);
                else
                    return 0;
            } else if (_spdim == 3) {
                if (type == celltype_t::tetrahedron)
                    return (_nth + 1) * (_nth + 2) * (_nth + 3) / 6;
                else if (type == celltype_t::hexahedron)
                    return (_nth + 1) * (_nth + 1) * (_nth + 1);
                else
                    return 0;
            } else {
                return 0;
            }
        }

        static constexpr unsigned int __spdim(celltype_t type) {
            if ((type == celltype_t::triangle) || (type == celltype_t::quadrilateral))
                return 2;
            else if ((type == celltype_t::tetrahedron) || (type == celltype_t::hexahedron))
                return 3;
            else
                return 0;
        }
        
        static constexpr bool __valid_celltype(unsigned int _spdim, celltype_t type) {
            return ((_spdim == 2 && (type == celltype_t::triangle || type == celltype_t::quadrilateral)) ||
                    (_spdim == 3 && (type == celltype_t::tetrahedron || type == celltype_t::hexahedron)));
        }

    };

} // end namespace mesh_entity

#endif