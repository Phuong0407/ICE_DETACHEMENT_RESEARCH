#ifndef MESH_GENERATION_VECTOR_HPP
#define MESH_GENERATION_VECTOR_HPP

#include "geometry.hpp"
#include "point.hpp"

#include <array>
#include <cmath>
#include <type_traits>
#include <iostream>

namespace mesh_generation
{
    template <unsigned int _spdim, typename double_t = double>
    class vector
    {
        static_assert(_spdim == 1 || _spdim == 2 || _spdim == 3, "spatial dimension must be 1, 2, or 3.");
        static_assert(std::is_floating_point_v<double_t>, "coordinate type must be a floating point type.");

    private:
        std::array<double_t, _spdim> components;

    public:
        vector() : components{} {}

        explicit vector(const std::array<double_t, _spdim>& vals) : components(vals) {}

        explicit vector(const point<double_t, _spdim>& P) {
            for (unsigned int i = 0; i < _spdim; ++i)
                components[i] = P(i);
        }

        template <typename... Args, typename = std::enable_if_t<(sizeof...(Args) == _spdim) &&
                                      (std::conjunction_v<std::is_convertible<Args, double_t>...>)>>
        explicit vector(Args... args) : components{static_cast<double_t>(args)...} {}

        vector(const vector&) = default;
        vector(vector&&) noexcept = default;
        vector& operator=(const vector&) = default;
        vector& operator=(vector&&) noexcept = default;
        ~vector() = default;

        void set(unsigned int i, double_t val) {
            if (i >= _spdim) {
                throw std::out_of_range("index out of bounds in point coordinate component-access.");
            }
            components[i] = val;
        }

        bool __parallel__(const vector& other) const {
            double_t ratio = 0;
            bool first_valid_ratio = true;
            for (unsigned int i = 0; i < _spdim; ++i) {
                double_t a = components[i];
                double_t b = other(i);
                if (std::abs(a) < ETOL<double_t> && std::abs(b) >= ETOL<double_t>)
                    return false;
                if (std::abs(a) >= ETOL<double_t> && std::abs(b) < ETOL<double_t>)
                    return false;
                else if (std::abs(a) <= ETOL<double_t> && std::abs(b) <= ETOL<double_t>)
                    continue;
                double_t current_ratio = a / b;
                if (first_valid_ratio) {
                    ratio = current_ratio;
                    first_valid_ratio = false;
                } else {
                    if (std::abs(current_ratio - ratio) > ETOL<double_t>)
                        return false;
                }
            }
            return true;
        }

        [[nodiscard]] vector operator+(const vector& other) const {
            vector result;
            for (unsigned int i = 0; i < _spdim; ++i)
                result.components[i] = components[i] + other.components[i];
            return result;
        }

        [[nodiscard]] vector operator-(const vector& other) const {
            vector result;
            for (unsigned int i = 0; i < _spdim; ++i)
                result.components[i] = components[i] - other.components[i];
            return result;
        }

        [[nodiscard]] vector operator*(double_t scalar) const {
            vector result;
            for (unsigned int i = 0; i < _spdim; ++i)
                result.components[i] = components[i] * scalar;
            return result;
        }

        [[nodiscard]] vector operator/(double_t scalar) const {
            if (scalar == 0.0)
                throw std::domain_error("division by zero");
            vector result;
            for (unsigned int i = 0; i < _spdim; ++i)
                result.components[i] = components[i] / scalar;
            return result;
        }

        [[nodiscard]] vector& operator+=(const vector& other) {
            for (unsigned int i = 0; i < _spdim; ++i)
                components[i] += other.components[i];
            return *this;
        }

        [[nodiscard]] vector& operator-=(const vector& other) {
            for (unsigned int i = 0; i < _spdim; ++i)
                components[i] -= other.components[i];
            return *this;
        }

        [[nodiscard]] vector& operator*=(double_t scalar) {
            for (unsigned int i = 0; i < _spdim; ++i)
                components[i] *= scalar;
            return *this;
        }

        [[nodiscard]] vector& operator/=(double_t scalar) {
            if (scalar == 0.0)
                throw std::domain_error("Division by zero");
            for (unsigned int i = 0; i < _spdim; ++i)
                components[i] /= scalar;
            return *this;
        }

        double_t dot(const vector& other) const {
            double_t sum = 0.0;
            for (unsigned int i = 0; i < _spdim; ++i)
                sum += components[i] * other.components[i];
            return sum;
        }

        template <unsigned int _dim = _spdim, typename std::enable_if_t<_dim == 3>* = nullptr>
        [[nodiscard]] vector cross(const vector& other) const {
            return vector(
                components[1] * other.components[2] - components[2] * other.components[1],
                components[2] * other.components[0] - components[0] * other.components[2],
                components[0] * other.components[1] - components[1] * other.components[0]
            );
        }
        template <unsigned int _dim = _spdim, typename std::enable_if_t<_dim == 2>* = nullptr>
        double_t cross(const vector& other) const {
            return components[0] * other.components[1] - components[1] * other.components[0];
        }

        double_t magnitude() const {
            return std::sqrt(dot(*this));
        }

        [[nodiscard]] vector normalize() const {
            double_t mag = magnitude();
            if (mag == 0.0)
                throw std::domain_error("cannot normalize a zero vector");
            return *this / mag;
        }

        double_t operator()(unsigned int i) const {
            if (i >= _spdim)
                throw std::out_of_range("index out of bounds");
            return components[i];
        }

        [[nodiscard]] std::array<double_t, _spdim>& operator()() const { return components; }

        void print() const {
            std::cout << "(";
            for (size_t i = 0; i < _spdim; ++i)
                std::cout << components[i] << (i < _spdim - 1 ? ", " : "");
            std::cout << ")\n";
        }

    };

    template <unsigned int _spdim, typename double_t>
    vector<_spdim, double_t> operator*(double_t scalar, const vector<_spdim, double_t>& v) {
        return v * scalar;
    }

} // namespace mesh_generation

#endif