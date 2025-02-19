#ifndef DATA_STRUCTURE_TENSOR_OMP_HPP
#define DATA_STRUCTURE_TENSOR_OMP_HPP
#define USING_OPENMP

#include "tensor_base.hpp"


#ifdef USING_OPENMP
#include <omp.h>
#endif

namespace data_structure
{

    template<unsigned int Rank, unsigned int Dim, typename real_t = double, typename index_t = unsigned int>
    class tensor_omp : public tensor_base<Rank, Dim, real_t> {
        private:
            std::vector<real_t> data_;
            std::array<index_t, Rank> shape_;

            index_t _idx_(const std::array<index_t, Rank>& ids) const {
                index_t index = 0, stride = 1;
                for (index_t i = Rank; i-- > 0;) {
                    if (ids[i] >= shape_[i])
                        throw std::out_of_range("index out of bounds");
                    index += ids[i] * stride;
                    stride *= shape_[i];
                }
                return index;
            }

        public:
            using tsomp_ptr = std::unique_ptr<tensor_omp>;
            
            const std::vector<real_t>& get_data() const override { return data_; }
            unsigned int& data_size() const override { return data_.size(); }

            tensor_omp() {
                index_t total_size = 1;
                shape_.fill(Dim);
                for (index_t s : shape_) total_size *= s;
                data_.resize(total_size, 0.0);
            }

            void init_tensor() override {
                #ifdef USING_OPENMP
                #pragma omp parallel for
                #endif
                for (size_t i = 0; i < data_.size(); ++i) {
                    data_[i] = static_cast<real_t>(0.0);
                }
            }

            void init(const real_t* raw_data, size_t size) override {
                if (size != data_.size())
                    throw std::invalid_argument("size mismatch in raw pointer initialization.");
                #ifdef USING_OPENMP
                #pragma omp parallel for
                #endif
                for (size_t i = 0; i < data_.size(); ++i) {
                    data_[i] = raw_data[i];
                }
            }

            void init(const std::vector<real_t>& vec) override {
                if (vec.size() != data_.size())
                    throw std::invalid_argument("size mismatch in vector initialization.");
                #ifdef USING_OPENMP
                #pragma omp parallel for
                #endif
                for (size_t i = 0; i < data_.size(); ++i) {
                    data_[i] = vec[i];
                }
            }

            void init(std::initializer_list<real_t> list) override {
                if (list.size() != data_.size())
                    throw std::invalid_argument("size mismatch in initializer_list initialization.");
                #ifdef USING_OPENMP
                #pragma omp parallel for
                #endif
                size_t i = 0;
                for (const auto& val : list) {
                    data_[i++] = val;
                }
            }

            typename tensor_base<Rank, Dim, real_t>::ts_ptr operator+(const tensor_base<Rank, Dim, real_t>& other) const override {
                auto result = std::make_unique<tensor_omp<Rank, Dim, real_t>>();
                #ifdef USING_OPENMP
                #pragma omp parallel for
                #endif
                
                std::vector<real_t> other_data_ = other.get_data();
                for (size_t i = 0; i < data_.size(); ++i) {
                    result->data_[i] += other_data_[i];
                }
                return result;
            }

            typename tensor_base<Rank, Dim, real_t>::ts_ptr operator-(const tensor_base<Rank, Dim, real_t>& other) const override {
                auto result = std::make_unique<tensor_omp<Rank, Dim, real_t>>();
                #ifdef USING_OPENMP
                #pragma omp parallel for
                #endif
                
                std::vector<real_t> other_data_ = other.get_data();
                for (size_t i = 0; i < data_.size(); ++i) {
                    result->data_[i] -= other_data_[i];
                }
                return result;
            }

            typename tensor_base<Rank, Dim, real_t>::ts_ptr operator*(const tensor_base<Rank, Dim, real_t>& other) const override {
                auto result = std::make_unique<tensor_omp<Rank, Dim, real_t>>();
                #ifdef USING_OPENMP
                #pragma omp parallel for
                #endif
                
                std::vector<real_t> other_data_ = other.get_data();
                for (size_t i = 0; i < data_.size(); ++i) {
                    result->data_[i] *= other_data_[i];
                }
                return result;
            }

            template<unsigned int Rank2, unsigned int Dim2>
            std::unique_ptr<tensor_omp<Rank + Rank2, Dim, real_t>> outer(const tensor_base<Rank2, Dim2, real_t>& other) const override {

                const auto* other_cast = dynamic_cast<const tensor_omp<Rank2, Dim2, real_t>*>(&other);
                if (!other_cast)
                    throw std::invalid_argument("invalid tensor type in outer product");

                auto result = std::make_unique<tensor_omp<Rank + Rank2, Dim, real_t>>();

                const auto& other_data = other_cast->get_data();

                #ifdef USING_OPENMP
                #pragma omp parallel for collapse(2)
                #endif
                for (size_t i = 0; i < data_.size(); ++i) {
                    for (size_t j = 0; j < other_data.size(); ++j) {
                        result->data_[i * other_data.size() + j] = data_[i] * other_data[j];
                    }
                }

                return result;
            }

            template<std::size_t Rank2, std::size_t Dim2>
            std::unique_ptr<tensor_base<Rank + Rank2 - 2, Dim, real_t>> 
            dot(const tensor_base<Rank2, Dim2, real_t>& other) const override 
            {
                if (this->size() != other.size()) {
                    throw std::invalid_argument("Tensor sizes do not match for dot product.");
                }
            
                const std::vector<real_t>& other_data = other.get_data();
                std::size_t total_size = this->size();
            
                real_t global_dot = 0.0;
            
                #ifdef USING_OPENMP
                #pragma omp parallel for reduction(+:global_dot)
                #endif
                for (std::size_t i = 0; i < total_size; ++i) {
                    global_dot += this->data_[i] * other_data[i];
                }
            
                auto result_tensor = std::make_unique<tensor_omp<Rank + Rank2 - 2, Dim, real_t>>();
                result_tensor->init({global_dot});
            
                return result_tensor;
            }
            
    };

} // namespace data_structure

#endif