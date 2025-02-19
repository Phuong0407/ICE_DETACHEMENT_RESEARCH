#ifndef MESH_GENERATION_MATRIX_HPP
#define MESH_GENERATION_MATRIX_HPP

#include <array>
#include <stdexcept>
#include <iostream>
#include "vector.hpp"

namespace mesh_generation {

    template<unsigned int N, typename double_t = double>
    class matrix {
    private:
        std::array<std::array<double_t, N>, N> data{};

    public:
        matrix() {
            for (auto &row : data) row.fill(0);
        }

        static matrix<N, double_t> identity() {
            matrix<N, double_t> I;
            for (unsigned int i = 0; i < N; ++i) {
                I(i, i) = 1;
            }
            return I;
        }

        double_t& operator()(unsigned int row, unsigned int col) {
            return data[row][col];
        }

        const double_t& operator()(unsigned int row, unsigned int col) const {
            return data[row][col];
        }

        vector<N, double_t> operator*(const vector<N, double_t>& v) const {
            vector<N, double_t> result;
            for (unsigned int i = 0; i < N; ++i) {
                for (unsigned int j = 0; j < N; ++j) {
                    result[i] += data[i][j] * v[j];
                }
            }
            return result;
        }

        matrix<N, double_t> operator*(const matrix<N, double_t>& other) const {
            matrix<N, double_t> result;
            for (unsigned int i = 0; i < N; ++i) {
                for (unsigned int j = 0; j < N; ++j) {
                    for (unsigned int k = 0; k < N; ++k) {
                        result(i, j) += data[i][k] * other(k, j);
                    }
                }
            }
            return result;
        }

        matrix<2, double_t> inverse() const {
            static_assert(N == 2, "Inverse function is implemented only for 2x2 matrices.");
            matrix<2, double_t> result;
            double_t det = data[0][0] * data[1][1] - data[0][1] * data[1][0];

            if (det == 0) {
                throw std::runtime_error("Matrix is singular and cannot be inverted.");
            }

            double_t inv_det = 1.0 / det;
            result(0, 0) = data[1][1] * inv_det;
            result(0, 1) = -data[0][1] * inv_det;
            result(1, 0) = -data[1][0] * inv_det;
            result(1, 1) = data[0][0] * inv_det;

            return result;
        }

        void print() const {
            for (const auto &row : data) {
                for (const auto &elem : row) {
                    std::cout << elem << " ";
                }
                std::cout << std::endl;
            }
        }
    };

} // namespace mesh_generation

#endif