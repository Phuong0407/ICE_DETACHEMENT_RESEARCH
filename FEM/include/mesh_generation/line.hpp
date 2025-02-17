#ifndef MESH_GENERATION_LINE_HPP
#define MESH_GENERATION_LINE_HPP

#include "geometry.hpp"
#include "point.hpp"
#include "vector.hpp"

#include <optional>
#include <cmath>

namespace mesh_generation
{
    template<unsigned int _spdim, typename double_t>
    class line {
    private:
        static_assert(_spdim == 1 || _spdim == 2 || _spdim == 3, "spatial dimension must be 1, 2, or 3.");
        static_assert(std::is_floating_point_v<double_t>, "coordinate type must be a floating point type.");
        
        point<_spdim, double_t> p;
        vector<_spdim, double_t> v;

    public:
        line() : p(), v() {}
        ~line() = default;

        line(point<_spdim, double_t> p, vector<_spdim, double_t> v) : p(p), v(v) {}
        line(point<_spdim, double_t> p1, point<_spdim, double_t> p2) : p(p1), v(p2 - p1) {}

        line(const line& other) = default;
        line(line&& other) noexcept = default;
        line& operator=(const line& other) = default;
        line& operator=(line&& other) noexcept = default;

        [[nodiscard]] const point<_spdim, double_t>& __point__() const { return p; }
        [[nodiscard]] const vector<_spdim, double_t>& __dir__() const { return v; }
        [[nodiscard]] point<_spdim, double_t> point_at(double_t t) const { return p + v * t; }

        [[nodiscard]] bool operator==(const line& other) const {
            if (!v.__parallel__(other.dir()))
                return false;
            vector<_spdim, double_t> diff = other.__point__() - p;
            return diff.__parallel__(v);
        }

        [[nodiscard]] bool operator!=(const line& other) const { return !(*this == other); }
        [[nodiscard]] bool __parallel__(const line& other) const { return v.__parallel__(other.dir()); }

        [[nodiscard]] bool __intersect__(const line& other) const {
            if (__parallel__(other))
                return false;

            if constexpr (_spdim == 2) {
                double_t det = v.cross(other.dir());
                return std::abs(det) > ETOL<double_t>;
            }
            if constexpr (_spdim == 3) {
                return !__skew__(other);
            }
            return false;
        }

        [[nodiscard]] std::optional<point<_spdim, double_t>> intersect_pont(const line& other) const {
            if constexpr (_spdim == 2) {
                double_t det = v.cross(other.dir());
                if (std::abs(det) < ETOL<double_t>) {
                    return std::nullopt;
                }

                vector<2, double_t> diff = other.__point__() - p;
                double_t t = diff.cross(other.dir()) / det;
                return p + v * t;
            }
            return std::nullopt;
        }

        [[nodiscard]] bool __skew__(const line& other) const {
            if constexpr (_spdim == 3) {
                if (__parallel__(other))
                    return false;

                vector<3, double_t> cross_dir = v.cross(other.dir());
                vector<3, double_t> diff = other.__point__() - p;
                return std::abs(diff.dot(cross_dir)) > ETOL<double_t>;
            }
            return false;
        }
    };

} // namespace mesh_generation

#endif