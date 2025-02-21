#ifndef DATA_STRUCTURE_TENSOR_TENSOR_OMP_HPP
#define DATA_STRUCTURE_TENSOR_TENSOR_OMP_HPP

#include <data_structure/tensor/tensor_base.hpp>

#ifdef USING_OPENMP
#include <omp.h>
#endif

namespace data_structure
{
    
    template<unsigned int Rank, unsigned int Dim, typename real_t = double, typename index_t = unsigned int>
    class tensor_omp : public tensor_base<Rank, Dim, real_t, index_t> {
        public:
            typename tensor_base<Rank, Dim, real_t, index_t>::ts_ptr operator+(const tensor_base<Rank, Dim, real_t>& other) const override {
                auto result = std::make_unique<tensor_omp<Rank, Dim, real_t>>();
                std::vector<real_t> other_data_ = other.get_data();

                #ifdef USING_OPENMP
                #pragma omp parallel for
                #endif
                for (index_t i = 0; i < data_.size(); ++i) {
                    result->data_[i] = this->data_[i] + other_data_[i];
                }
                return result;
            }

            typename tensor_base<Rank, Dim, real_t, index_t>::ts_ptr operator-(const tensor_base<Rank, Dim, real_t>& other) const override {
                auto result = std::make_unique<tensor_omp<Rank, Dim, real_t>>();
                
                std::vector<real_t> other_data_ = other.get_data();

                #ifdef USING_OPENMP
                #pragma omp parallel for
                #endif
                for (index_t i = 0; i < data_.size(); ++i) {
                    result->data_[i] = this->data_[i] - other_data_[i];
                }
                return result;
            }

            typename tensor_base<Rank, Dim, real_t, index_t>::ts_ptr operator*(const tensor_base<Rank, Dim, real_t>& other) const override {
                auto result = std::make_unique<tensor_omp<Rank, Dim, real_t>>();
            
                std::vector<real_t> other_data_ = other.get_data();

                #ifdef USING_OPENMP
                #pragma omp parallel for
                #endif
                for (index_t i = 0; i < data_.size(); ++i) {
                    result->data_[i] = this->data_[i] * other_data_[i];
                }
                return result;
            }
    };


} // namespace data_structure



#endif