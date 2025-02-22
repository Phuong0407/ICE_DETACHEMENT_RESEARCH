#ifndef GEOMETRY_LINE_HPP
#define GEOMETRY_LINE_HPP

#include <geometry/geometry.hpp>

#include <optional>
#include <cmath>

namespace geometry
{
    template<unsigned int spdim, typename real_t>
    class line {
    private:
        static_assert(spdim == 1 || spdim == 2 || spdim == 3, "spatial dimension must be 1, 2, or 3.");
        static_assert(std::is_floating_point_v<real_t>, "coordinate type must be a floating point type.");
        
        point<spdim, real_t> p;
        vector<spdim, real_t> v;

    public:
        line() : p(), v() {}
        ~line() = default;

        line(point<spdim, real_t> p, vector<spdim, real_t> v) : p(p), v(v) {}
        line(point<spdim, real_t> p1, point<spdim, real_t> p2) : p(p1), v(p2 - p1) {}

        line(const line& other) = default;
        line(line&& other) noexcept = default;
        line& operator=(const line& other) = default;
        line& operator=(line&& other) noexcept = default;

        [[nodiscard]] const point<spdim, real_t>& __point__() const { return p; }
        [[nodiscard]] const vector<spdim, real_t>& __dir__() const { return v; }
        [[nodiscard]] point<spdim, real_t> point_at(real_t t) const { return p + v * t; }

        [[nodiscard]] bool operator==(const line& other) const {
            if (!v.__parallel__(other.dir()))
                return false;
            vector<spdim, real_t> diff = other.__point__() - p;
            return diff.__parallel__(v);
        }

        [[nodiscard]] bool operator!=(const line& other) const { return !(*this == other); }
        [[nodiscard]] bool __parallel__(const line& other) const { return v.__parallel__(other.dir()); }

        [[nodiscard]] bool __intersect__(const line& other) const {
            if (__parallel__(other))
                return false;

            if constexpr (spdim == 2) {
                real_t det = v.cross(other.dir());
                return std::abs(det) > ETOL<real_t>;
            }
            if constexpr (spdim == 3) {
                return !__skew__(other);
            }
            return false;
        }

        [[nodiscard]] std::optional<point<spdim, real_t>> intersect_pont(const line& other) const {
            if constexpr (spdim == 2) {
                real_t det = v.cross(other.dir());
                if (std::abs(det) < ETOL<real_t>) {
                    return std::nullopt;
                }

                vector<2, real_t> diff = other.__point__() - p;
                real_t t = diff.cross(other.dir()) / det;
                return p + v * t;
            }
            return std::nullopt;
        }

        [[nodiscard]] bool __skew__(const line& other) const {
            if constexpr (spdim == 3) {
                if (__parallel__(other))
                    return false;

                vector<3, real_t> cross_dir = v.cross(other.dir());
                vector<3, real_t> diff = other.__point__() - p;
                return std::abs(diff.dot(cross_dir)) > ETOL<real_t>;
            }
            return false;
        }
    };

} // namespace geometry

#endif