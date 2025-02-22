#ifndef GEOMETRY_VECTOR_HPP
#define GEOMETRY_VECTOR_HPP

#include <geometry/geometry.hpp>

#include <array>
#include <cmath>
#include <type_traits>
#include <iostream>

namespace geometry
{
    template <unsigned int spdim, typename real_t = double>
    class   vector
    {
        static_assert(spdim == 1 || spdim == 2 || spdim == 3, "spatial dimension must be 1, 2, or 3.");
        static_assert(std::is_floating_point_v<real_t>, "coordinate type must be a floating point type.");

    private:
        std::array<real_t, spdim> components;

    public:
        vector() : components{} {}

        explicit vector(const std::array<real_t, spdim>& vals) : components(vals) {}

        explicit vector(const point<spdim, real_t>& P) {
            for (unsigned int i = 0; i < spdim; ++i)
                components[i] = P(i);
        }

        template <typename... Args, typename = std::enable_if_t<(sizeof...(Args) == spdim) &&
                                      (std::conjunction_v<std::is_convertible<Args, real_t>...>)>>
        explicit vector(Args... args) : components{static_cast<real_t>(args)...} {}

        vector(const vector&) = default;
        vector(vector&&) noexcept = default;
        vector& operator=(const vector&) = default;
        vector& operator=(vector&&) noexcept = default;
        ~vector() = default;

        void set(unsigned int i, real_t val) {
            if (i >= spdim) {
                throw std::out_of_range("index out of bounds in point coordinate component-access.");
            }
            components[i] = val;
        }

        bool is_parallel(const vector& other) const {
            real_t ratio = 0;
            bool first_valid_ratio = true;
            for (unsigned int i = 0; i < spdim; ++i) {
                real_t a = components[i];
                real_t b = other(i);
                if (std::abs(a) < ETOL<real_t> && std::abs(b) >= ETOL<real_t>)
                    return false;
                if (std::abs(a) >= ETOL<real_t> && std::abs(b) < ETOL<real_t>)
                    return false;
                else if (std::abs(a) <= ETOL<real_t> && std::abs(b) <= ETOL<real_t>)
                    continue;
                real_t current_ratio = a / b;
                if (first_valid_ratio) {
                    ratio = current_ratio;
                    first_valid_ratio = false;
                } else {
                    if (std::abs(current_ratio - ratio) > ETOL<real_t>)
                        return false;
                }
            }
            return true;
        }

        [[nodiscard]] vector operator+(const vector& other) const {
            vector result;
            for (unsigned int i = 0; i < spdim; ++i)
                result.components[i] = components[i] + other.components[i];
            return result;
        }

        [[nodiscard]] vector operator-(const vector& other) const {
            vector result;
            for (unsigned int i = 0; i < spdim; ++i)
                result.components[i] = components[i] - other.components[i];
            return result;
        }

        [[nodiscard]] vector operator*(real_t scalar) const {
            vector result;
            for (unsigned int i = 0; i < spdim; ++i)
                result.components[i] = components[i] * scalar;
            return result;
        }

        [[nodiscard]] vector operator/(real_t scalar) const {
            if (scalar == 0.0)
                throw std::domain_error("division by zero");
            vector result;
            for (unsigned int i = 0; i < spdim; ++i)
                result.components[i] = components[i] / scalar;
            return result;
        }

        [[nodiscard]] vector& operator+=(const vector& other) {
            for (unsigned int i = 0; i < spdim; ++i)
                components[i] += other.components[i];
            return *this;
        }

        [[nodiscard]] vector& operator-=(const vector& other) {
            for (unsigned int i = 0; i < spdim; ++i)
                components[i] -= other.components[i];
            return *this;
        }

        [[nodiscard]] vector& operator*=(real_t scalar) {
            for (unsigned int i = 0; i < spdim; ++i)
                components[i] *= scalar;
            return *this;
        }

        [[nodiscard]] vector& operator/=(real_t scalar) {
            if (scalar == 0.0)
                throw std::domain_error("Division by zero");
            for (unsigned int i = 0; i < spdim; ++i)
                components[i] /= scalar;
            return *this;
        }

        real_t dot(const vector& other) const {
            real_t sum = 0.0;
            for (unsigned int i = 0; i < spdim; ++i)
                sum += components[i] * other.components[i];
            return sum;
        }

        template <unsigned int _dim = spdim, typename std::enable_if_t<_dim == 3>* = nullptr>
        [[nodiscard]] vector cross(const vector& other) const {
            return vector(
                components[1] * other.components[2] - components[2] * other.components[1],
                components[2] * other.components[0] - components[0] * other.components[2],
                components[0] * other.components[1] - components[1] * other.components[0]
            );
        }
        template <unsigned int _dim = spdim, typename std::enable_if_t<_dim == 2>* = nullptr>
        real_t cross(const vector& other) const {
            return components[0] * other.components[1] - components[1] * other.components[0];
        }

        real_t magnitude() const {
            return std::sqrt(dot(*this));
        }

        [[nodiscard]] vector normalize() const {
            real_t mag = magnitude();
            if (mag == 0.0)
                throw std::domain_error("cannot normalize a zero vector");
            return *this / mag;
        }

        real_t operator()(unsigned int i) const {
            if (i >= spdim)
                throw std::out_of_range("index out of bounds");
            return components[i];
        }

        [[nodiscard]] const std::array<real_t, spdim>& operator()() const { return components; }
        [[nodiscard]] const std::array<real_t, spdim>& get_data() const { return components; }

        void print() const {
            std::cout << "(";
            for (size_t i = 0; i < spdim; ++i)
                std::cout << components[i] << (i < spdim - 1 ? ", " : "");
            std::cout << ")\n";
        }

    };

    template <unsigned int spdim, typename real_t>
    vector<spdim, real_t> operator*(real_t scalar, const vector<spdim, real_t>& v) {
        return v * scalar;
    }

} // namespace geometry

#endif