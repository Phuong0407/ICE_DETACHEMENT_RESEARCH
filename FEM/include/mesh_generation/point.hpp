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

#ifndef MESH_GENERATION_POINT_HPP
#define MESH_GENERATION_POINT_HPP

#include "geometry.hpp"
#include "vector.hpp"
#include "edge.hpp"

#include <cstdarg>
#include <type_traits>
#include <array>
#include <vector>
#include <initializer_list>
#include <cmath>
#include <stdexcept>
#include <iostream>

namespace mesh_generation
{
    template<unsigned int _spdim, typename double_t = double>
    class point
    {
        static_assert(_spdim == 1 || _spdim == 2 || _spdim == 3, "spatial dimension must be 1, 2, or 3.");
        static_assert(std::is_floating_point_v<double_t>, "coordinate type must be a floating point type.");

        private:
            std::array<double_t, _spdim> x;

        public:
            point() : x{} {}

            explicit point(const double_t (&coords)[_spdim]) {
                for (unsigned int i = 0; i < _spdim; ++i) {
                    x[i] = coords[i];
                }
            }

            explicit point(const std::vector<double_t>& coords) {
                if (coords.size() != _spdim) {
                    throw std::invalid_argument("incorrect number of coordinates.");
                }
                for (unsigned int i = 0; i < _spdim; ++i) {
                    x[i] = coords[i];
                }
            }

            explicit point(const vector<double_t, _spdim>& v) {
                for (unsigned int i = 0; i < _spdim; ++i)
                    x[i] = v(i);
            }

            template<typename... Args,
                    typename = std::enable_if_t<(sizeof...(Args) == _spdim)
                        && (std::conjunction_v<std::is_convertible<Args, double_t>...>)>>
            explicit point(Args... args) : x{static_cast<double_t>(args)...} {}

            point(const point& other) = default;
            point(point&& other) noexcept = default;
            point& operator=(const point& other) = default;
            point& operator=(point&& other) noexcept = default;

            bool operator==(const point& other) const {
                for (unsigned int i = 0; i < _spdim; ++i) {
                    if (x[i] != other.x[i]) {
                        return false;
                    }
                }
                return true;
            }

            bool operator!=(const point& other) const {
                return !(*this == other);
            }

            void set(unsigned int i, double_t val) {
                if (i >= _spdim) {
                    throw std::out_of_range("index out of bounds in point coordinate component-access.");
                }
                x[i] = val;
            }

            [[nodiscard]] const std::array<double_t, _spdim>& operator()()const { return x; }

            double_t dist(const point& other) const {
                double_t sum = 0.0;
                for (unsigned int i = 0; i < _spdim; ++i) {
                    double_t diff = x[i] - other.x[i];
                    sum += diff * diff;
                }
                return std::sqrt(sum);
            }

            double_t operator()(unsigned int i) const {
                if (i >= _spdim) {
                    throw std::out_of_range("index out of bounds in point coordinate component-access.");
                }
                return x[i];
            }
            
            constexpr unsigned int dim() const { return _spdim; }

            void print() const {
                std::cout << "(";
                for (unsigned int i = 0; i < _spdim - 1; ++i) {
                    std::cout << x[i] << " ";
                } std::cout << x[_spdim - 1];
                std::cout << ")";
            }


            point operator+(const vector<_spdim, double_t>& vec) const {
                point result;
                for (unsigned int i = 0; i < _spdim; ++i)
                    result.set(i, x[i] + vec(i));
                return result;
            }
    
            point operator-(const vector<_spdim, double_t>& vec) const {
                point result;
                for (unsigned int i = 0; i < _spdim; ++i)
                    result.set(i, x[i] - vec(i));
                return result;
            }

            point& operator+=(const vector<_spdim, double_t>& vec) {
                for (unsigned int i = 0; i < _spdim; ++i)
                    this->x[i] += vec(i);
                return *this;
            }
    
            point& operator-=(const vector<_spdim, double_t>& vec) {
                for (unsigned int i = 0; i < _spdim; ++i)
                    this->x[i] -= vec(i);
                return *this;
            }
            
            vector<_spdim, double_t> operator-(const point<_spdim, double_t>& other) const {
                vector<_spdim, double_t> result;
                for (unsigned int i = 0; i < _spdim; ++i)
                    result.set(i, x[i] - other(i));
                return result;
            }

            point<_spdim, double_t> project_to_edge(const edge<_spdim, double_t>& e) const
            {
                const point<_spdim, double_t>& A = e.first();
                const point<_spdim, double_t>& B = e.second();
            
                vector<_spdim, double_t> AB = B - A;
                vector<_spdim, double_t> AP = (*this) - A;
            
                double_t dot_AB_AB = AB.dot(AB);
                double_t dot_AP_AB = AP.dot(AB);
            
                if (std::abs(dot_AB_AB) < ETOL<double_t>) return A;
            
                double_t t = dot_AP_AB / dot_AB_AB;
                t = std::max(double_t(0), std::min(static_cast<double_t>(1), t));
                return A + AB * t;
            }
            
            point<_spdim, double_t> project_to_edge(const point<_spdim, double_t>& A, const point<_spdim, double_t>& B) const
            {
                vector<_spdim, double_t> AB = B - A;
                vector<_spdim, double_t> AP = (*this) - A;
            
                double_t dot_AB_AB = AB.dot(AB);
                double_t dot_AP_AB = AP.dot(AB);
            
                if (std::abs(dot_AB_AB) < ETOL<double_t>) return A;
            
                double_t t = dot_AP_AB / dot_AB_AB;
                t = std::max(double_t(0), std::min(static_cast<double_t>(1), t));
                return A + AB * t;
            }

            static point<_spdim, double_t> mid_point(const point<_spdim, double_t>& A, const point<_spdim, double_t>& B) {
                point<_spdim, double_t> _mid_point;
                for (unsigned int i = 0; i < _spdim; ++i) {
                    _mid_point.set(i, (A(i) + B(i)) / 2.0);
                }
                return _mid_point;
            }

    };

} // namespace mesh_generation

#endif