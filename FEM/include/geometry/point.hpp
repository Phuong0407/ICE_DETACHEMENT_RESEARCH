/**
 * @file point.hpp
 * @author Phuong Diep (diepthanhphuong0407@gmail.com)
 * @brief 
 * @version 0.1
 * @date 2025-02-16
 * 
 * @copyright Copyright (c) 2025
 * 
 */

#ifndef GEOMETRY_POINT_HPP
#define GEOMETRY_POINT_HPP

#include <geometry/geometry.hpp>

#include <cstdarg>
#include <type_traits>
#include <array>
#include <vector>
#include <initializer_list>
#include <cmath>
#include <stdexcept>
#include <iostream>

namespace geometry
{
    template<unsigned int spdim, typename real_t = double>
    class point
    {
        static_assert(spdim == 1 || spdim == 2 || spdim == 3, "spatial dimension must be 1, 2, or 3.");
        static_assert(std::is_floating_point_v<real_t>, "coordinate type must be a floating point type.");

        private:
            std::array<real_t, spdim> x;

        public:
            point() : x{} {}

            explicit point(const real_t (&coords)[spdim]) {
                for (unsigned int i = 0; i < spdim; ++i) {
                    x[i] = coords[i];
                }
            }

            explicit point(const std::vector<real_t>& coords) {
                if (coords.size() != spdim) {
                    throw std::invalid_argument("incorrect number of coordinates.");
                }
                for (unsigned int i = 0; i < spdim; ++i) {
                    x[i] = coords[i];
                }
            }

            explicit point(const vector<spdim, real_t>& v) {
                for (unsigned int i = 0; i < spdim; ++i)
                    x[i] = v(i);
            }

            template<typename... Args,
                    typename = std::enable_if_t<(sizeof...(Args) == spdim)
                        && (std::conjunction_v<std::is_convertible<Args, real_t>...>)>>
            explicit point(Args... args) : x{static_cast<real_t>(args)...} {}

            point(const point& other) = default;
            point(point&& other) noexcept = default;
            point& operator=(const point& other) = default;
            point& operator=(point&& other) noexcept = default;

            bool operator==(const point& other) const {
                for (unsigned int i = 0; i < spdim; ++i) {
                    if (x[i] != other.x[i]) {
                        return false;
                    }
                }
                return true;
            }

            bool operator!=(const point& other) const {
                return !(*this == other);
            }

            void set(unsigned int i, real_t val) {
                if (i >= spdim) {
                    throw std::out_of_range("index out of bounds in point coordinate component-access.");
                }
                x[i] = val;
            }

            [[nodiscard]] const std::array<real_t, spdim>& operator()()const { return x; }

            real_t dist(const point& other) const {
                real_t sum = 0.0;
                for (unsigned int i = 0; i < spdim; ++i) {
                    real_t diff = x[i] - other.x[i];
                    sum += diff * diff;
                }
                return std::sqrt(sum);
            }

            real_t operator()(unsigned int i) const {
                if (i >= spdim) {
                    throw std::out_of_range("index out of bounds in point coordinate component-access.");
                }
                return x[i];
            }
            
            constexpr unsigned int dim() const { return spdim; }

            void print() const {
                std::cout << "(";
                for (unsigned int i = 0; i < spdim - 1; ++i) {
                    std::cout << x[i] << " ";
                } std::cout << x[spdim - 1];
                std::cout << ")";
            }


            point operator+(const vector<spdim, real_t>& vec) const {
                point result;
                for (unsigned int i = 0; i < spdim; ++i)
                    result.set(i, x[i] + vec(i));
                return result;
            }
    
            point operator-(const vector<spdim, real_t>& vec) const {
                point result;
                for (unsigned int i = 0; i < spdim; ++i)
                    result.set(i, x[i] - vec(i));
                return result;
            }

            point& operator+=(const vector<spdim, real_t>& vec) {
                for (unsigned int i = 0; i < spdim; ++i)
                    this->x[i] += vec(i);
                return *this;
            }
    
            point& operator-=(const vector<spdim, real_t>& vec) {
                for (unsigned int i = 0; i < spdim; ++i)
                    this->x[i] -= vec(i);
                return *this;
            }
            
            vector<spdim, real_t> operator-(const point<spdim, real_t>& other) const {
                vector<spdim, real_t> result;
                for (unsigned int i = 0; i < spdim; ++i)
                    result.set(i, x[i] - other(i));
                return result;
            }

            point<spdim, real_t> project_to_edge(const edge<spdim, real_t>& e) const
            {
                const point<spdim, real_t>& A = e.first();
                const point<spdim, real_t>& B = e.second();
            
                vector<spdim, real_t> AB = B - A;
                vector<spdim, real_t> AP = (*this) - A;
            
                real_t dot_AB_AB = AB.dot(AB);
                real_t dot_AP_AB = AP.dot(AB);
            
                if (std::abs(dot_AB_AB) < ETOL<real_t>) return A;
            
                real_t t = dot_AP_AB / dot_AB_AB;
                t = std::max(real_t(0), std::min(static_cast<real_t>(1), t));
                return A + AB * t;
            }
            
            point<spdim, real_t> project_to_edge(const point<spdim, real_t>& A, const point<spdim, real_t>& B) const
            {
                vector<spdim, real_t> AB = B - A;
                vector<spdim, real_t> AP = (*this) - A;
            
                real_t dot_AB_AB = AB.dot(AB);
                real_t dot_AP_AB = AP.dot(AB);
            
                if (std::abs(dot_AB_AB) < ETOL<real_t>) return A;
            
                real_t t = dot_AP_AB / dot_AB_AB;
                t = std::max(real_t(0), std::min(static_cast<real_t>(1), t));
                return A + AB * t;
            }

            static point<spdim, real_t> mid_point(const point<spdim, real_t>& A, const point<spdim, real_t>& B) {
                point<spdim, real_t> _mid_point;
                for (unsigned int i = 0; i < spdim; ++i) {
                    _mid_point.set(i, (A(i) + B(i)) / 2.0);
                }
                return _mid_point;
            }

    };

} // namespace geometry

#endif