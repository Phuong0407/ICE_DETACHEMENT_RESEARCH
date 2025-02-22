#ifndef MESH_GENERATION_DELAUNAY_DELAUNAY_TRIANGULATAION_HPP
#define MESH_GENERATION_DELAUNAY_DELAUNAY_TRIANGULATAION_HPP

#include <geometry/polygon.hpp>


namespace mesh_generation
{
    template<unsigned int spdim = 2, typename real_t = double, typename index_t = unsigned int>
    class delaunay_triangulation{
        private:
            geometry::polygon<spdim, real_t> domain;
            

        public:
            delaunay_triangulation() = default;
            ~delaunay_triangulation() = default;

    };

} // namespace mesh_generation


#endif