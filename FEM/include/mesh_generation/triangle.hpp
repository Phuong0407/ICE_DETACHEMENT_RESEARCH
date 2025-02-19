#ifndef MESH_GENERATION_TRIANGLE_HPP
#define MESH_GENERATION_TRIANGLE_HPP

#include "polygon.hpp"

namespace mesh_generation
{

    template<unsigned int _spdim, typename double_t = double>
    class triangle : public polygon<_spdim, double_t> 
    {
    public:
        triangle() = default;

        explicit triangle(const std::array<point<_spdim, double_t>, 3>& points)
            : polygon<_spdim, double_t>(std::vector<point<_spdim, double_t>>(points.begin(), points.end())) {}

        explicit triangle(const point<_spdim, double_t>& p1,
                          const point<_spdim, double_t>& p2,
                          const point<_spdim, double_t>& p3)
            : polygon<_spdim, double_t>({p1, p2, p3}) {}

        [[nodiscard]] double_t area() const {
            if constexpr (_spdim == 2) 
            {
                const auto& p1 = this->vertices[0];
                const auto& p2 = this->vertices[1];
                const auto& p3 = this->vertices[2];

                return std::abs((p1.x() * (p2.y() - p3.y()) +
                                 p2.x() * (p3.y() - p1.y()) +
                                 p3.x() * (p1.y() - p2.y())) * 0.5);
            }
            else if constexpr (_spdim == 3) 
            {
                vector<3, double_t> v1 = this->vertices[1] - this->vertices[0];
                vector<3, double_t> v2 = this->vertices[2] - this->vertices[0];

                return v1.cross(v2).length() * 0.5;
            }
            else 
            {
                throw std::logic_error("area calculation is only available for 2D and 3D triangles.");
            }
        }

        [[nodiscard]] bool __valid__() const 
        {
            return this->vertices.size() == 3;
        }

        [[nodiscard]] const point<_spdim, double_t>>& first() const {
            return this->get_vertex(0);
        }
        [[nodiscard]] const point<_spdim, double_t>>& second() const {
            return this->get_vertex(1);
        }
        [[nodiscard]] const point<_spdim, double_t>>& third() const {
            return this->get_vertex(2);
        }
    };

} // namespace mesh_generation

#endif