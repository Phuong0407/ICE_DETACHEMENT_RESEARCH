#ifndef MESH_GENERATION_POLYGON_HPP
#define MESH_GENERATION_POLYGON_HPP

#include "geometry.hpp"
#include "point.hpp"
#include "vector.hpp"
#include "edge.hpp"

#include <algorithm>
#include <stdexcept>
#include <optional>

namespace mesh_generation
{
    template<unsigned int _spdim, typename double_t = double>
    class polygon
    {
    private:
        std::vector<point<_spdim, double_t>> vertices;
        mutable std::optional<bool> is_positive_orientation;
        mutable std::optional<vector<_spdim, double_t>> v_norm;

    public:
        polygon() = default;

        polygon(const std::vector<point<_spdim, double_t>>& points)
        {
            if (points.size() < 3)
                throw std::invalid_argument("a polygon must have at least three vertices.");
            if (constexpr (_spdim == 3)) {
                if (!__coplanar__(points))
                    throw std::invalid_argument("the given points are not coplanar.");
            }
            vertices = std::move(points);
        }

        void add_vertex(const point<_spdim, double_t>& p) {
            vertices.push_back(p);
        }

        [[nodiscard]] static bool __coplanar__(const std::vector<point<_spdim, double_t>>& points) const
        {
            if constexpr (_spdim == 2) return true;
            if (vertices.size() == 3) return true;

            vector<_spdim, double_t> v1 = points[1] - points[0];
            vector<_spdim, double_t> v2 = points[2] - points[0];
            vector<_spdim, double_t> normal = v1.cross(v2);

            for (size_t i = 3; i < points.size(); ++i)
            {
                vector<_spdim, double_t> v = points[i] - points[0];
                if (std::abs(normal.dot(v)) > etol)
                    return false;
            }
            return true;
        }

        [[nodiscard]] vector<3, double_t> normal() const
        {
            if constexpr (_spdim != 3)
                throw std::logic_error("normal vector is only defined for 3D polygons.");

            vector<3, double_t> v1 = vertices[1] - vertices[0];
            vector<3, double_t> v2 = vertices[2] - vertices[0];
            return v1.cross(v2).normalize();
        }

        [[nodiscard]] double_t area() const
        {
            if constexpr (_spdim != 2)
                throw std::logic_error("area is only defined for 2D polygons.");

            double_t sum = 0.0;
            size_t n = vertices.size();

            for (size_t i = 0; i < n; ++i)
            {
                const point<2, double_t>& p1 = vertices[i];
                const point<2, double_t>& p2 = vertices[(i + 1) % n];
                sum += (p1.x() * p2.y() - p2.x() * p1.y());
            }

            return std::abs(sum) * 0.5;
        }

        [[nodiscard]] vector<3, double_t> norm() const
        {
            if constexpr (_spdim != 3)
                throw std::logic_error("normal vector is only defined for 3D polygons.");

            if (!v_norm.has_value())
            {
                vector<3, double_t> v1 = vertices[1] - vertices[0];
                vector<3, double_t> v2 = vertices[2] - vertices[0];
                v_norm = v1.cross(v2).normalize();
            }
            return v_norm.value();
        }


        [[nodiscard]] bool __oriented_clockwise__(std::optional<vector<3, double_t>> reference_normal = std::nullopt) const
        {
            if constexpr (_spdim == 2) 
            {
                double_t sum = 0.0;
                size_t n = vertices.size();
        
                for (size_t i = 0; i < n; ++i)
                {
                    const point<2, double_t>& p1 = vertices[i];
                    const point<2, double_t>& p2 = vertices[(i + 1) % n];
                    sum += (p2.x() - p1.x()) * (p2.y() + p1.y());
                }
        
                return sum > 0;
            }
            else if constexpr (_spdim == 3)
            {
                if (!this->v_norm.has_value()) {
                    this->v_norm = this->normal();
                    std::cerr << "[Warning] no reference normal found. Computed default normal.\n";
                }
        
                if (reference_normal.has_value()) {
                    return this->v_norm.value().dot(reference_normal.value()) > 0;
                }
                else {
                    return this->v_norm.value().z() < 0;
                }
            }
            else
            {
                throw std::logic_error("orientation check is only available for 2D and 3D polygons.");
            }
        }
        

        void flip_orientation() { std::reverse(vertices.begin(), vertices.end()); }

        

    };
}

#endif