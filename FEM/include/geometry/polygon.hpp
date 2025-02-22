#ifndef GEOMETRY_POLYGON_HPP
#define GEOMETRY_POLYGON_HPP

#include <geometry/geometry.hpp>

#include <algorithm>
#include <stdexcept>
#include <optional>

namespace geometry
{
    template<unsigned int spdim, typename real_t = double>
    class polygon
    {
    private:
        std::vector<point<spdim, real_t>> vertices;
        mutable std::optional<bool> is_positive_orientation;
        mutable std::optional<vector<spdim, real_t>> v_norm;

    public:
        polygon() = default;

        polygon(const std::vector<point<spdim, real_t>>& points)
        {
            if (points.size() < 3)
                throw std::invalid_argument("a polygon must have at least three vertices.");
            if constexpr (spdim == 3) {
                if (!__coplanar__(points))
                    throw std::invalid_argument("the given points are not coplanar.");
            }
            vertices = std::move(points);
        }

        void add_vertex(const point<spdim, real_t>& p) {
            vertices.push_back(p);
        }

        template<typename... Args, typename = std::enable_if_t<(sizeof...(Args) == spdim) && (std::conjunction_v<std::is_convertible<Args, real_t>...>)>>
        void add_vertex(Args... args) {
            point<spdim, real_t> p(args...);
            vertices.push_back(p);
        }

        [[nodiscard]] static bool __coplanar__(const std::vector<point<spdim, real_t>>& points)
        {
            if constexpr (spdim == 2) return true;
            if (points.size() == 3) return true;

            vector<spdim, real_t> v1 = points[1] - points[0];
            vector<spdim, real_t> v2 = points[2] - points[0];
            vector<spdim, real_t> normal = v1.cross(v2);

            for (size_t i = 3; i < points.size(); ++i)
            {
                vector<spdim, real_t> v = points[i] - points[0];
                if (std::abs(normal.dot(v)) > ETOL<real_t>)
                    return false;
            }
            return true;
        }

        [[nodiscard]] vector<3, real_t> normal() const
        {
            if constexpr (spdim != 3)
                throw std::logic_error("normal vector is only defined for 3D polygons.");

            vector<3, real_t> v1 = vertices[1] - vertices[0];
            vector<3, real_t> v2 = vertices[2] - vertices[0];
            return v1.cross(v2).normalize();
        }

        [[nodiscard]] real_t area() const
        {
            if constexpr (spdim == 2) 
            {
                real_t sum = 0.0;
                size_t n = vertices.size();
        
                for (size_t i = 0; i < n; ++i)
                {
                    const point<2, real_t>& p1 = vertices[i];
                    const point<2, real_t>& p2 = vertices[(i + 1) % n];
                    sum += (p1(0) * p2(1) - p2(0) * p1(1));
                }
                return std::abs(sum) * 0.5;
            }
            else if constexpr (spdim == 3) 
            {
                vector<3, real_t> normal_vec = normal();
                real_t sum = 0.0;
                size_t n = vertices.size();
        
                for (size_t i = 1; i < n - 1; ++i)
                {
                    vector<3, real_t> v1 = vertices[i] - vertices[0];
                    vector<3, real_t> v2 = vertices[i + 1] - vertices[0];
                    sum += (v1.cross(v2)).dot(normal_vec);
                }
        
                return std::abs(sum) * 0.5;
            }
            else 
            {
                throw std::logic_error("area calculation is only available for 2D and 3D polygons.");
            }
        }        

        [[nodiscard]] bool __oriented_clockwise__(std::optional<vector<3, real_t>> reference_normal = std::nullopt) const
        {
            if constexpr (spdim == 2) 
            {
                real_t sum = 0.0;
                size_t n = vertices.size();
        
                for (size_t i = 0; i < n; ++i)
                {
                    const point<2, real_t>& p1 = vertices[i];
                    const point<2, real_t>& p2 = vertices[(i + 1) % n];
                    sum += (p1(0) * p2(1) - p2(0) * p1(1));
                }
        
                return sum > 0;
            }
            else if constexpr (spdim == 3)
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

        void print() const {
            std::cout << "polygon: ";
            for (const auto& p : vertices) {
                p.print();
                std::cout << " ";
            }
            std::cout << std::endl;
        }

        bool __contains__(const point<spdim, real_t>& P) const {
            int count = 0;
            int n = vertices.size();
            
            for (size_t i = 0, j = n - 1; i < n; j = i++) {
                const point<spdim, real_t>& A = vertices[i];
                const point<spdim, real_t>& B = vertices[j];

                if ((A[1] > P[1]) != (B[1] > P[1]) &&
                    (P[0] < (B[0] - A[0]) * (P[1] - A[1]) / (B[1] - A[1]) + A[0])) {
                    count++;
                }
            }
            return (count % 2) == 1;
        }

        [[nodiscard]] point<spdim, real_t> nearest_point_on_boundary(const point<spdim, real_t>& P) const {
            point<spdim, real_t> closestPoint;
            real_t minDistance = std::numeric_limits<real_t>::max();

            for (size_t i = 0; i < vertices.size(); ++i) {
                const point<spdim, real_t>& A = vertices[i];
                const point<spdim, real_t>& B = vertices[(i + 1) % vertices.size()];

                point<spdim, real_t> proj = P.project_to_edge(A, B);
                
                real_t distance = P.dist(proj);

                if (distance < minDistance) {
                    minDistance = distance;
                    closestPoint = proj;
                }
            }
            return closestPoint;
        }

        [[nodiscard]] const std::vector<edge<spdim, real_t>> get_boundary_edge() const {
            std::vector<edge<spdim, real_t>> edges;
            for (size_t i = 0; i < vertices.size(); ++i) {
                edges.emplace_back(vertices[i], vertices[(i + 1) % vertices.size()]);
            }
            return edges;
        }

        [[nodiscard]] const point<spdim, real_t>& get_vertex(unsigned int i) const {
            if (i >= spdim) {
                throw std::out_of_range("index out of bounds in point coordinate component-access.");
            }
            return vertices[i];
        }
    };
}

#endif