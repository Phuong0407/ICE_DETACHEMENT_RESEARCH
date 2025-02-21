#ifndef DATA_STRUCTURE_MATRIX_MATRIX_HPP
#define DATA_STRUCTURE_MATRIX_MATRIX_HPP

#include <config.h>

#include <array>
#include <vector>
#include <iostream>
#include <iomanip>
#include <stdexcept>
#include <initializer_list>

namespace data_structure
{
    
    template<unsigned int M, unsigned int N, typename real_t = double, typename index_t = unsigned int>
    class matrix_base {
        private:
            std::array<real_t, M * N> _data;

        public:
            matrix_base() { _data.fill(0); }
            ~matrix_base() = default;

            matrix_base(std::array<real_t, M * N> data) : _data(data) {}

            matrix_base(std::initializer_list<real_t> vals) {
                index_t i = 0;
                for (real_t val : vals) {
                    if (i >= M * N)
                        break;
                    _data[i++] = val;
                }
            }

            matrix_base(std::vector<real_t> vec) {
                if (vec.size() != M * N)
                    throw std::invalid_argument("the number of elements in input vector mismatches the size of matrix in matrix(vector).");
                for (index_t i = 0; i < M * N; ++i) {
                    _data[i] = vec[i];
                }
            }

            matrix_base(const matrix_base& other) = default;
            matrix_base& operator=(const matrix_base& other) = default;
            matrix_base(matrix_base&& other) noexcept = default;
            matrix_base& operator=(matrix_base&& other) noexcept = default;

            inline real_t& operator()(index_t i, index_t j) {
                if (i >= M || j >= N)
                    throw std::out_of_range("matrix index out of bounds");
                return _data[i * N + j];
            }

            inline const real_t operator()(index_t i, index_t j) const {
                if (i >= M || j >= N)
                    throw std::out_of_range("matrix index out of bounds");
                return _data[i * N + j];
            }

            inline real_t& operator()(index_t i) {
                if (i >= M * N)
                    throw std::out_of_range("matrix index out of bounds");
                return _data[i];
            }

            const real_t operator()(index_t i) const {
                if (i >= M * N)
                    throw std::out_of_range("matrix index out of bounds");
                return _data[i];
            }

            inline constexpr index_t row_size() const { return M; }
            inline constexpr index_t col_size() const { return N; }
            const std::array<real_t, M * N> & get_data() const { return _data; }
            std::array<real_t, M * N> & get_data() { return _data; }

            matrix_base<M, N, real_t, index_t> operator+(const matrix_base &other) const {
                const unsigned int size = M * N;
                matrix_base<M, N, real_t, index_t> result;
            
                std::array<real_t, size> &result_data = result.get_data();
                const std::array<real_t, size> &this_data = this->get_data();
                const std::array<real_t, size> &other_data = other.get_data();
            
                #ifdef USING_BLAS
                    std::copy(other_data.begin(), other_data.end(), result_data.begin());
            
                    if constexpr (std::is_same_v<real_t, float>) {
                        cblas_saxpy(static_cast<int>(size), 1.0f, this_data.data(), 1, result_data.data(), 1);
                        return result;                        
                    } else if constexpr (std::is_same_v<real_t, double>) {
                        cblas_daxpy(static_cast<int>(size), 1.0, this_data.data(), 1, result_data.data(), 1);
                        return result;
                    } else {
                        #ifdef USING_OPENMP
                            goto OPENMP_ADDITION;
                        #endif
                    }
                #endif

                #if defined(USING_OPENMP) && !defined(USING_BLAS)
                OPENMP_ADDITION:
                    #pragma omp parallel for
                    for (index_t i = 0; i < size; ++i) {
                        result_data[i] = this_data[i] + other_data[i];
                    }
                    return result;
                #endif
                for (index_t i = 0; i < size; ++i) {
                    result_data[i] = this_data[i] + other_data[i];
                }
                return result;
            }

            matrix_base<M, N, real_t, index_t> operator-(const matrix_base &other) const {
                const unsigned int size = M * N;
                matrix_base<M, N, real_t, index_t> result;
            
                std::array<real_t, size> &result_data = result.get_data();
                const std::array<real_t, size> &this_data = this->get_data();
                const std::array<real_t, size> &other_data = other.get_data();
            
                #ifdef USING_BLAS
                    std::copy(other_data.begin(), other_data.end(), result_data.begin());
            
                    if constexpr (std::is_same_v<real_t, float>) {
                        cblas_saxpy(static_cast<int>(size), -1.0f, this_data.data(), 1, result_data.data(), 1);
                        return result;
                    } else if constexpr (std::is_same_v<real_t, double>) {
                        cblas_daxpy(static_cast<int>(size), -1.0, this_data.data(), 1, result_data.data(), 1);
                        return result;
                    } else {
                        #ifdef USING_OPENMP
                            goto OPENMP_SUBTRACTION;
                        #endif
                    }
                #endif

                #if defined(USING_OPENMP) && !defined(USING_BLAS)
                    OPENMP_SUBTRACTION:
                    #pragma omp parallel for
                    for (index_t i = 0; i < size; ++i) {
                        result_data[i] = this_data[i] - other_data[i];
                    }
                    return result;
                #endif
                for (index_t i = 0; i < size; ++i) {
                    result_data[i] = this_data[i] - other_data[i];
                }
                return result;
            }
            
            template <unsigned int P>
            matrix_base<M, P, real_t, index_t> operator*(const matrix_base<N, P, real_t, index_t> &other) const {
                matrix_base<M, P, real_t, index_t> result;

                std::array<real_t, M * P>& result_data = result.get_data();
                const std::array<real_t, M * N>& this_data = this->get_data();
                const std::array<real_t, N * P>& other_data = other.get_data();

                #ifdef USING_BLAS
                    if constexpr (std::is_same_v<real_t, float>) {
                        cblas_sgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans,
                                    M, P, N,
                                    1.0f, this_data.data(), N, other_data.data(), P,
                                    0.0f, result_data.data(), P);
                        return result;
                    } 
                    else if constexpr (std::is_same_v<real_t, double>) {
                        cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans,
                                    M, P, N,
                                    1.0, this_data.data(), N, other_data.data(), P,
                                    0.0, result_data.data(), P);
                        return result;
                    } 
                    else {
                        #ifdef USING_OPENMP
                            goto OPENMP_MULTIPLICATION;
                        #else
                            goto 
                        #endif
                    }
                #endif

                #if defined(USING_OPENMP) && !defined(USING_BLAS)
                OPENMP_MULTIPLICATION:
                #pragma omp parallel for collapse(2)
                    for (index_t i = 0; i < M; ++i) {
                        for (index_t j = 0; j < P; ++j) {
                            real_t sum = 0;
                            for (index_t k = 0; k < N; ++k) {
                                sum += this_data[i * N + k] * other_data[k * P + j];
                            }
                            result_data[i * P + j] = sum;
                        }
                    }
                    return result;
                #endif

                for (index_t i = 0; i < M; ++i) {
                    for (index_t j = 0; j < P; ++j) {
                        real_t sum = 0.0;
                        for (index_t k = 0; k < N; ++k) {
                            sum += this_data[i * N + k] * other_data[k * P + j];
                        }
                        result_data[i * P + j] = sum;
                    }
                }
                return result;
            }



    };

} // namespace data_structure    

#endif