#ifndef MESH_GENERATION_EDGE_HPP
#define MESH_GENERATION_EDGE_HPP

#include "point.hpp"
#include "vector.hpp"

#include <type_traits>
#include <iostream>
#include <algorithm>

namespace mesh_generation
{
    template <unsigned int _spdim, typename double_t = double>
    class edge
    {
        static_assert(_spdim == 1 || _spdim == 2 || _spdim == 3, "spatial dimension must be 1, 2, or 3.");
        static_assert(std::is_floating_point_v<double_t>, "coordinate type must be a floating point type.");

    private:
        point<_spdim, double_t> p1;
        point<_spdim, double_t> p2;

    public:
        edge() : p1(), p2() {}
        ~edge() = default;

        edge(const point<_spdim, double_t>& p1, const point<_spdim, double_t>& p2) : p1(p1), p2(p2) {}

        template <typename... Args1, typename... Args2, 
                  typename = std::enable_if_t<(sizeof...(Args1) == _spdim) && (sizeof...(Args2) == _spdim) &&
                                              (std::conjunction_v<std::is_convertible<Args1, double_t>...>) &&
                                              (std::conjunction_v<std::is_convertible<Args2, double_t>...>)>>
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

        [[nodiscard]] const point<_spdim, double_t>& first() const { return p1; }
        [[nodiscard]] const point<_spdim, double_t>& second() const { return p2; }

        constexpr unsigned int dim() const { return _spdim; }

        double_t length() const {
            return p1.dist(p2);
        }

        void print() const {
            p1.print();
            std::cout << "<--->";
            p2.print();
        }
    };

} // namespace mesh_generation

#endif // MESH_GENERATION_EDGE_HPP
