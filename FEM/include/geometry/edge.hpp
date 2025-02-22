#ifndef GEOMETRY_EDGE_HPP
#define GEOMETRY_EDGE_HPP

#include <geometry/geometry.hpp>

#include <type_traits>
#include <iostream>
#include <algorithm>

namespace geometry
{
    template <unsigned int spdim, typename real_t = double>
    class edge
    {
        static_assert(spdim == 1 || spdim == 2 || spdim == 3, "spatial dimension must be 1, 2, or 3.");
        static_assert(std::is_floating_point_v<real_t>, "coordinate type must be a floating point type.");

    private:
        point<spdim, real_t> p1;
        point<spdim, real_t> p2;

    public:
        edge() : p1(), p2() {}
        ~edge() = default;

        explicit edge(const point<spdim, real_t>& p1, const point<spdim, real_t>& p2) : p1(p1), p2(p2) {}

        template <typename... Args1, typename... Args2, 
                  typename = std::enable_if_t<(sizeof...(Args1) == spdim) && (sizeof...(Args2) == spdim) &&
                                              (std::conjunction_v<std::is_convertible<Args1, real_t>...>) &&
                                              (std::conjunction_v<std::is_convertible<Args2, real_t>...>)>>
        edge(Args1... args1, Args2... args2) : p1(args1...), p2(args2...) {}

        edge(const edge& other) = default;
        edge(edge&& other) noexcept = default;
        edge& operator=(const edge& other) = default;
        edge& operator=(edge&& other) noexcept = default;

        bool operator==(const edge& other) const {
            return (p1 == other.p1 && p2 == other.p2) || (p1 == other.p2 && p2 == other.p1);
        }
        bool operator!=(const edge& other) const {
            return !(*this == other);
        }

        [[nodiscard]] const point<spdim, real_t>& first() const { return p1; }
        [[nodiscard]] const point<spdim, real_t>& second() const { return p2; }

        constexpr unsigned int dim() const { return spdim; }

        real_t length() const {
            return p1.dist(p2);
        }

        void print() const {
            p1.print();
            std::cout << "<--->";
            p2.print();
        }

        [[nodiscard]] const vector<spdim, real_t>& perpendicular() const {
            static_assert(spdim == 2, "This method is only valid for 2D edges.");

            real_t dx = p2(0) - p1(0);
            real_t dy = p2(1) - p1(1);

            return vector<spdim, real_t>(-dy, dx);
        }

        [[nodiscard]] bool is_intersected_with(const edge& other) const {
            static_assert(spdim == 2, "intersection is only available in 2D");
            auto orientation = [](const point& p, const point& q, const point& r) -> int {
                double val = (q.y - p.y) * (r.x - q.x) - (q.x - p.x) * (r.y - q.y);
                if (val == 0) return 0;
                return (val > 0) ? 1 : -1;
            };

            int O1 = orientation(A, B, other.A);
            int O2 = orientation(A, B, other.B);
            int O3 = orientation(other.A, other.B, A);
            int O4 = orientation(other.A, other.B, B);

            if (std::abs(O1 - O2) <= ETOL<real_t> && std::abs(O3 - O4) <= ETOL<real_t>)
                return true;

            auto on_segment = [](const point& p, const point& q, const point& r) -> bool {
                return std::min(p.x, r.x) <= q.x && q.x <= std::max(p.x, r.x) &&
                        std::min(p.y, r.y) <= q.y && q.y <= std::max(p.y, r.y);
            };
    
            if (O1 == 0 && on_segment(A, other.A, B)) return true;
            if (O2 == 0 && on_segment(A, other.B, B)) return true;
            if (O3 == 0 && on_segment(other.A, A, other.B)) return true;
            if (O4 == 0 && on_segment(other.A, B, other.B)) return true;
        
            return false;
        }
    };

} // namespace geometry

#endif // GEOMETRY_EDGE_HPP
