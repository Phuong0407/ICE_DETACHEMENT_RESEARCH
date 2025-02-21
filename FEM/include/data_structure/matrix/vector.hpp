#ifndef DATA_STRUCTURE_MATRIX_VECTOR_HPP
#define DATA_STRUCTURE_MATRIX_VECTOR_HPP

#include <data_structure/matrix/matrix.hpp>

#include <cmath>

namespace data_structure
{
    
    template<unsigned int M, typename real_t = double, typename index_t = unsigned int>
    class vector : protected matrix<M, 1, real_t, index_t> {
        public:
            using matrix<M, 1, real_t, index_t>::operator+;
            using matrix<M, 1, real_t, index_t>::operator-;
            using matrix<M, 1, real_t, index_t>::operator*;

            vector() { this->get_data().fill(0); }
            ~vector() = default;
            
            vector(const vector& other) = default;
            vector& operator=(const vector& other) = default;
            vector(vector&& other) noexcept = default;
            vector& operator=(vector&& other) noexcept = default;

            vector(std::array<real_t, M> data) {
                std::copy(data.begin(), data.end(), this->get_data().begin());
            }
        
            vector(std::initializer_list<real_t> vals) {
                index_t i = 0;
                auto& data = this->get_data();
                for (real_t val : vals) {
                    if (i >= M)
                        break;
                    data[i++] = val;
                }
            }
        
            vector(std::vector<real_t> vec) {
                if (vec.size() != M)
                    throw std::invalid_argument("vector size does not match the expected size.");
                auto& data = this->get_data();
                for (index_t i = 0; i < M; ++i) {
                    data[i] = vec[i];
                }
            }
        
            vector(const matrix<M, 1, real_t, index_t>& vec) : matrix<M, 1, real_t, index_t>(vec) {}

            real_t& operator()(index_t i) { return this->get_data()[i]; }
            const real_t& operator()(index_t i) const { return this->get_data()[i]; }
            constexpr index_t size() const { return M; }
            
            vector operator+(const vector& other) const {
                return this->matrix<M, 1, real_t, index_t>::operator+(other);
            }
        
            vector operator-(const vector& other) const {
                return this->matrix<M, 1, real_t, index_t>::operator-(other);
            }
        
            vector operator*(real_t scalar) const {
                return this->matrix<M, 1, real_t, index_t>::operator*(scalar);
            }

            real_t norm() const {
                const std::array<real_t, M> data = this->get_data();
                #ifdef USING_BLAS
                    if constexpr (std::is_same_v<real_t, float>) {
                        return cblas_snrm2(static_cast<int>(M), data.data(), 1);
                    } else if constexpr (std::is_same_v<real_t, double>) {
                        return cblas_dnrm2(static_cast<int>(M), data.data(), 1);
                    }
                #endif

                real_t sum = 0;
                #ifdef USING_OPENMP
                    #pragma omp parallel for reduction(+:sum)
                    for (index_t i = 0; i < M; ++i) {
                        sum += data[i] * data[i];
                    }
                    return std::sqrt(sum);
                #endif
                for (index_t i = 0; i < M; ++i) {
                    sum += data[i] * data[i];
                }
                return std::sqrt(sum);
            }

            template<unsigned int N>
            matrix<M, N, real_t, index_t> outer(const vector<N, real_t, index_t>& u) {
                matrix<M, N, real_t, index_t> result;

                const std::array<real_t, M>& v_data = this->get_data();
                const std::array<real_t, N>& u_data = u.get_data();
                std::array<real_t, M * N>& result_data = result.get_data();

                #ifdef USING_BLAS
                    if constexpr (std::is_same_v<real_t, float>) {
                        cblas_sger(CblasRowMajor, M, N, 1.0f, v_data.data(), 1, u_data.data(), 1, result_data.data(), N);
                    } else if constexpr (std::is_same_v<real_t, double>) {
                        cblas_dger(CblasRowMajor, M, N, 1.0, v_data.data(), 1, u_data.data(), 1, result_data.data(), N);
                    }
                    return result;
                #endif

                #ifdef USING_OPENMP
                    #pragma omp parallel for collapse(2)
                #endif
                for (index_t i = 0; i < M; ++i) {
                    for (index_t j = 0; j < N; ++j) {
                        result_data[i * N + j] = v_data[i] * u_data[j];
                    }
                }
                return result;
            }


    };

} // namespace data_structure


#endif