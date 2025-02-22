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

    extern "C" {
        void dgetrf_(int* M, int* N, double* A, int* lda, int* ipiv, int* info);
    }

    template<unsigned int M, unsigned int N, typename real_t = double, typename index_t = unsigned int>
    class matrix {
        private:
            std::array<real_t, M * N> _data;

        public:
            matrix() { _data.fill(0); }
            ~matrix() = default;

            matrix(std::array<real_t, M * N> data) : _data(data) {}

            matrix(std::initializer_list<real_t> vals) {
                index_t i = 0;
                for (real_t val : vals) {
                    if (i >= M * N)
                        break;
                    _data[i++] = val;
                }
            }

            matrix(std::vector<real_t> vec) {
                if (vec.size() != M * N)
                    throw std::invalid_argument("the number of elements in input vector mismatches the size of matrix in matrix(vector).");
                for (index_t i = 0; i < M * N; ++i) {
                    _data[i] = vec[i];
                }
            }

            matrix(const matrix& other) = default;
            matrix& operator=(const matrix& other) = default;
            matrix(matrix&& other) noexcept = default;
            matrix& operator=(matrix&& other) noexcept = default;

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

            void print() const {
                for (unsigned int i = 0; i < M; ++i) {
                    for (unsigned int j = 0; j < N - 1; ++j) {
                        std::cout << this->operator()(i, j) << "   ";
                    }
                    std::cout << this->operator()(i, N - 1) << std::endl;
                }
            }

            inline constexpr index_t row_size() const { return M; }
            inline constexpr index_t col_size() const { return N; }
            const std::array<real_t, M * N> & get_data() const { return _data; }
            std::array<real_t, M * N> & get_data() { return _data; }

            matrix<M, N, real_t, index_t> operator+(const matrix &other) const {
                const unsigned int size = M * N;
                matrix<M, N, real_t, index_t> result;
            
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
                    }
                #endif

                #if defined(USING_OPENMP)
                    #pragma omp parallel for
                #endif
                for (index_t i = 0; i < size; ++i) {
                    result_data[i] = this_data[i] + other_data[i];
                }
                return result;
            }

            matrix<M, N, real_t, index_t>& operator+=(const matrix<M, N, real_t, index_t>& other) {
                const std::array<real_t, M * N>& other_data = other.get_data();
        
                #ifdef USING_BLAS
                    if constexpr (std::is_same_v<real_t, float>) {
                        cblas_saxpy(static_cast<int>(M * N), 1.0f, other_data.data(), 1, _data.data(), 1);
                        return *this;
                    } else if constexpr (std::is_same_v<real_t, double>) {
                        cblas_daxpy(static_cast<int>(M * N), 1.0, other_data.data(), 1, _data.data(), 1);
                        return *this;                        
                    }
                #endif
                
                #ifdef USING_OPENMP
                    #pragma omp parallel for
                #endif
                for (index_t i = 0; i < M * N; ++i) {
                    _data[i] += other_data[i];
                }
                return *this;
            }

            matrix<M, N, real_t, index_t> operator-(const matrix &other) const {
                const unsigned int size = M * N;
                matrix<M, N, real_t, index_t> result;
            
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
                    }
                #endif

                #ifdef USING_OPENMP
                    #pragma omp parallel for
                #endif
                for (index_t i = 0; i < size; ++i) {
                    result_data[i] = this_data[i] - other_data[i];
                }
                return result;
            }
            
            matrix<M, N, real_t, index_t>& operator-=(const matrix<M, N, real_t, index_t>& other) {
                const std::array<real_t, M * N>& other_data = other.get_data();
        
                #ifdef USING_BLAS
                    if constexpr (std::is_same_v<real_t, float>) {
                        cblas_saxpy(static_cast<int>(M * N), -1.0f, other_data.data(), 1, _data.data(), 1);
                        return *this;
                    } else if constexpr (std::is_same_v<real_t, double>) {
                        cblas_daxpy(static_cast<int>(M * N), -1.0, other_data.data(), 1, _data.data(), 1);
                        return *this;                        
                    }
                #endif
                
                #ifdef USING_OPENMP
                    #pragma omp parallel for
                #endif
                for (index_t i = 0; i < M * N; ++i) {
                    _data[i] -= other_data[i];
                }
                return *this;
            }

            template <unsigned int P>
            matrix<M, P, real_t, index_t> operator*(const matrix<N, P, real_t, index_t> &other) const {
                matrix<M, P, real_t, index_t> result;

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
                #endif

                #if defined(USING_OPENMP)
                    #pragma omp parallel for collapse(2)
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

            matrix<M, N, real_t, index_t> operator*(real_t val) const {
                matrix<M, N, real_t, index_t> result;
                std::array<real_t, M * N>& result_data = result.get_data();
                const std::array<real_t, M * N>& this_data = this->get_data();

                #ifdef USING_BLAS
                    std::copy(this_data.begin(), this_data.end(), result_data.begin());

                    if constexpr (std::is_same_v<real_t, float>) {
                        cblas_sscal(static_cast<int>(M * N), val, result_data.data(), 1);
                        return result;
                    } 
                    else if constexpr (std::is_same_v<real_t, double>) {
                        cblas_dscal(static_cast<int>(M * N), val, result_data.data(), 1);
                        return result;
                    }
                #endif

                #ifdef USING_OPENMP
                    #pragma omp parallel for
                #endif

                for (index_t i = 0; i < M * N; ++i) {
                    result_data[i] = this_data[i] * val;
                }
                return result;
            }

            real_t trace() const {
                static_assert(M == N, "trace is only defined for square matrices.");

                const std::array<real_t, M * N>& data = this->get_data();
                real_t sum = 0;

                #ifdef USING_BLAS
                    if constexpr (std::is_same_v<real_t, float>) {
                        return cblas_sasum(static_cast<int>(M), data.data(), N + 1);
                    } else if constexpr (std::is_same_v<real_t, double>) {
                        return cblas_dasum(static_cast<int>(M), data.data(), N + 1);
                    }
                #endif

                #ifdef USING_OPENMP
                    #pragma omp parallel for reduction(+:sum)
                #endif
                for (index_t i = 0; i < M; ++i) {
                    sum += data[i * N + i];
                }

                return sum;
            }

            matrix<N, N, real_t, index_t> gram_matrix() const {
                matrix<N, N, real_t, index_t> result;
                const std::array<real_t, M * N>& data = this->get_data();
                std::array<real_t, N * N>& result_data = result.get_data();

                #ifdef USING_BLAS
                    if constexpr (std::is_same_v<real_t, float>) {
                        cblas_sgemm(CblasRowMajor, CblasTrans, CblasNoTrans,
                                    N, N, M,
                                    1.0f, data.data(), N, data.data(), N,
                                    0.0f, result_data.data(), N);
                        return result;                        
                    } else if constexpr (std::is_same_v<real_t, double>) {
                        cblas_dgemm(CblasRowMajor, CblasTrans, CblasNoTrans,
                                    N, N, M,
                                    1.0, data.data(), N, data.data(), N,
                                    0.0, result_data.data(), N);
                        return result;                        
                    }
                #endif

                #ifdef USING_OPENMP
                    #pragma omp parallel for collapse(2)
                #endif
                for (index_t i = 0; i < N; ++i) {
                    for (index_t j = 0; j < N; ++j) {
                        real_t sum = 0;
                        for (index_t k = 0; k < M; ++k) {
                            sum += data[k * N + i] * data[k * N + j];
                        }
                        result_data[i * N + j] = sum;
                    }
                }
                return result;
            }

            matrix<M, M, real_t, index_t> covariance_matrix() const {
                matrix<M, M, real_t, index_t> result;
                const std::array<real_t, M * N>& data = this->get_data();
                std::array<real_t, M * M>& result_data = result.get_data();

                #ifdef USING_BLAS
                    if constexpr (std::is_same_v<real_t, float>) {
                        cblas_sgemm(CblasRowMajor, CblasNoTrans, CblasTrans,
                                    M, M, N,
                                    1.0f, data.data(), N, data.data(), N,
                                    0.0f, result_data.data(), M);
                        return result;                        
                    } else if constexpr (std::is_same_v<real_t, double>) {
                        cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasTrans,
                                    M, M, N,
                                    1.0, data.data(), N, data.data(), N,
                                    0.0, result_data.data(), M);
                        return result;                        
                    }
                #endif

                #ifdef USING_OPENMP
                    #pragma omp parallel for collapse(2)
                #endif
                for (index_t i = 0; i < M; ++i) {
                    for (index_t j = 0; j < M; ++j) {
                        real_t sum = 0;
                        for (index_t k = 0; k < N; ++k) {
                            sum += data[i * N + k] * data[j * N + k];
                        }
                        result_data[i * M + j] = sum;
                    }
                }
                return result;
            }

            matrix<N, M, real_t, index_t> transpose() const {
                matrix<N, M, real_t, index_t> result;
                const std::array<real_t, M * N>& data = this->get_data();
                std::array<real_t, N * M>& result_data = result.get_data();
        
                #ifdef USING_OPENMP
                        #pragma omp parallel for collapse(2)
                #endif
                for (index_t i = 0; i < M; ++i) {
                    for (index_t j = 0; j < N; ++j) {
                        result_data[j * M + i] = data[i * N + j];
                    }
                }
                return result;
            }
            
            real_t det() {
                static_assert(M == N, "determinant is only defined for square matrices.");

                std::array<int, M> pivot;
                int info;
                int size = static_cast<int>(M);
                matrix<M, M, real_t, index_t> temp(*this);
                real_t det = 1.0;

                #ifdef USING_BLAS
                    dgetrf_(&size, &size, temp.get_data().data(), &size, pivot.data(), &info);

                    if (info != 0)
                        return 0.0;
                    for (index_t i = 0; i < M; ++i) {
                        if (pivot[i] != static_cast<int>(i + 1)) det = -det;
                        det *= temp(i, i);
                    }
                    return det;
                #endif

                #ifdef USING_OPENMP
                    for (index_t i = 0; i < M; ++i) {
                        index_t pivot = i;
                        for (index_t j = i + 1; j < M; ++j) {
                            if (std::abs(temp(j, i)) > std::abs(temp(pivot, i))) {
                                pivot = j;
                            }
                        }
                        if (std::abs(temp(pivot, i)) < 1e-9) return 0.0;

                        if (pivot != i) {
                            #pragma omp parallel for
                            for (index_t k = 0; k < M; ++k) std::swap(temp(i, k), temp(pivot, k));
                            det = -det;
                        }

                        det *= temp(i, i);
                        real_t diag = temp(i, i);
                        #pragma omp parallel for
                        for (index_t j = i + 1; j < M; ++j) {
                            real_t factor = temp(j, i) / diag;
                            for (index_t k = i; k < M; ++k) {
                                temp(j, k) -= factor * temp(i, k);
                            }
                        }
                    }
                    return det;
                #endif

                for (index_t i = 0; i < M; ++i) {
                    index_t pivot = i;
                    for (index_t j = i + 1; j < M; ++j) {
                        if (std::abs(temp(j, i)) > std::abs(temp(pivot, i))) {
                            pivot = j;
                        }
                    }
                    if (std::abs(temp(pivot, i)) < 1e-9) return 0.0;

                    if (pivot != i) {
                        for (index_t k = 0; k < M; ++k) std::swap(temp(i, k), temp(pivot, k));
                        det = -det;
                    }

                    det *= temp(i, i);
                    real_t diag = temp(i, i);
                    for (index_t j = i + 1; j < M; ++j) {
                        real_t factor = temp(j, i) / diag;
                        for (index_t k = i; k < M; ++k) {
                            temp(j, k) -= factor * temp(i, k);
                        }
                    }
                }
                return det;
            }

            bool is_symmetric() const {
                static_assert(M == N, "determinant is only defined for square matrices.");
                bool _is_symmetric = true;
                #ifdef USING_OPENMP
                    #pragma omp parallel for collapse(2)
                #endif

                for (index_t i = 0; i < M; ++i) {
                    for (index_t j = 0; j < M; ++j) {
                        if (std::abs(_data[i * M + j] - _data[j * M + i]) < ETOL<real_t>)
                            continue;
                        else
                            #ifdef USING_OPENMP
                                #pragma omp critical
                            #endif
                            _is_symmetric = false;
                    }
                }
                return _is_symmetric;
            }

            matrix<M, M, real_t, index_t> inverse() {
                if (std::abs(this->det()) <= ETOL<real_t>)
                    throw std::runtime_error("Matrix is not invertible.");
            
                Eigen::Matrix<real_t, M, M> A;
                #ifdef USING_OPENMP
                    #pragma omp parallel for collapse(2)
                #endif
                for (unsigned int i = 0; i < M; ++i) {
                    for (unsigned int j = 0; j < M; ++j) {
                        A(i, j) = this->operator()(i, j);
                    }
                }
            
                Eigen::Matrix<real_t, M, M> A_inverse = A.inverse();
            
                matrix<M, M, real_t, index_t> inv_matrix;
                #ifdef USING_OPENMP
                    #pragma omp parallel for collapse(2)
                #endif
                for (unsigned int i = 0; i < M; ++i) {
                    for (unsigned int j = 0; j < M; ++j) {
                        inv_matrix(i, j) = A_inverse(i, j);
                    }
                }
                return inv_matrix;
            }
            
    };

} // namespace data_structure    

#endif