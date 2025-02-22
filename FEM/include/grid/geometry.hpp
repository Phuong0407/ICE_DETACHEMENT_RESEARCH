#ifndef MESH_GEOMETRY_HPP
#define MESH_GEOMETRY_HPP

#include <vector>

namespace mesh_geometry {
    
    template<unsigned int spdim, typename index_t, typename number_traits>
    class mesh_geometry {
        private:
            static_assert(spdim == 1 || spdim == 2 || spdim == 3, "spatial dimension must be 1, 2, or 3.");

            index_t _n_node;

            // there is spacedim array, and each of them has size num_nodes
            std::vector<number_traits> x;
            std::vector<number_traits> y;
            // z is only available in 3d
            std::vector<number_traits> z;

    };

} // namespace mesh_geometry


#endif