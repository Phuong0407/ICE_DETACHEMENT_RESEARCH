#ifndef DATA_STRUCTURE_TENSOR_OMP_HPP
#define DATA_STRUCTURE_TENSOR_OMP_HPP

#include <data_structure/tensor/tensor_base.hpp>

#ifdef USING_OPENMP
#include <omp.h>
#endif

namespace data_structure
{

    template<unsigned int Rank, unsigned int Dim, typename real_t = double, typename index_t = unsigned int>
    class tensor_omp : public tensor_base<Rank, Dim, real_t, index_t> {
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
            
            void set_data(index_t idx, real_t val) override {
                if (idx >= data_.size())
                    throw std::invalid_argument("index to set data is out of bound.");
                this->data_[idx] = val;
            }
            const std::vector<real_t>& get_data() const override { return data_; }
            unsigned int data_size() const override { return data_.size(); }


            tensor_omp() {
                index_t total_size = 1;
                shape_.fill(Dim);
                for (index_t s : shape_) total_size *= s;
                data_.resize(total_size, 0.0);
            }

            void init() override {
                #ifdef USING_OPENMP
                #pragma omp parallel for
                #endif
                for (index_t i = 0; i < data_.size(); ++i) {
                    data_[i] = static_cast<real_t>(0.0);
                }
            }

            void init(const real_t* raw_data, index_t size) override {
                if (size != data_.size())
                    throw std::invalid_argument("size mismatch in raw pointer initialization.");
                #ifdef USING_OPENMP
                #pragma omp parallel for
                #endif
                for (index_t i = 0; i < data_.size(); ++i) {
                    data_[i] = raw_data[i];
                }
            }

            void init(const std::vector<real_t>& vec) override {
                if (vec.size() != data_.size())
                    throw std::invalid_argument("size mismatch in vector initialization.");
                #ifdef USING_OPENMP
                #pragma omp parallel for
                #endif
                for (index_t i = 0; i < data_.size(); ++i) {
                    data_[i] = vec[i];
                }
            }

            void init(std::initializer_list<real_t> list) override {
                if (list.size() != data_.size())
                    throw std::invalid_argument("size mismatch in initializer_list initialization.");
                const auto* list_ptr = list.begin();

                #ifdef USING_OPENMP
                #pragma omp parallel for
                #endif

                for (index_t i = 0; i < data_.size(); ++i) {
                    data_[i] = list_ptr[i];
                }
            }

            typename tensor_base<Rank, Dim, real_t>::ts_ptr operator+(const tensor_base<Rank, Dim, real_t>& other) const override {
                auto result = std::make_unique<tensor_omp<Rank, Dim, real_t>>();
                std::vector<real_t> other_data_ = other.get_data();

                #ifdef USING_OPENMP
                #pragma omp parallel for
                #endif
                
                for (index_t i = 0; i < data_.size(); ++i) {
                    result->data_[i] += other_data_[i];
                }
                return result;
            }

            typename tensor_base<Rank, Dim, real_t>::ts_ptr operator-(const tensor_base<Rank, Dim, real_t>& other) const override {
                auto result = std::make_unique<tensor_omp<Rank, Dim, real_t>>();
                
                std::vector<real_t> other_data_ = other.get_data();

                #ifdef USING_OPENMP
                #pragma omp parallel for
                #endif

                for (index_t i = 0; i < data_.size(); ++i) {
                    result->data_[i] -= other_data_[i];
                }
                return result;
            }

            typename tensor_base<Rank, Dim, real_t>::ts_ptr operator*(const tensor_base<Rank, Dim, real_t>& other) const override {
                auto result = std::make_unique<tensor_omp<Rank, Dim, real_t>>();
               
                std::vector<real_t> other_data_ = other.get_data();

                #ifdef USING_OPENMP
                #pragma omp parallel for
                #endif
                
                for (index_t i = 0; i < data_.size(); ++i) {
                    result->data_[i] *= other_data_[i];
                }
                return result;
            }

            template<unsigned int Rank2, unsigned int Dim2>
            std::unique_ptr<tensor_omp<Rank + Rank2, Dim, real_t>> outer(const tensor_base<Rank2, Dim2, real_t>& other) const {

                const auto* other_cast = dynamic_cast<const tensor_omp<Rank2, Dim2, real_t>*>(&other);
                if (!other_cast)
                    throw std::invalid_argument("invalid tensor type in outer product");

                auto result = std::make_unique<tensor_omp<Rank + Rank2, Dim, real_t>>();
                const auto& other_data = other_cast->get_data();

                #ifdef USING_OPENMP
                #pragma omp parallel for collapse(2)
                #endif

                for (index_t i = 0; i < data_.size(); ++i) {
                    for (index_t j = 0; j < other_data.size(); ++j) {
                        result->set_data(i * other_data.size() + j, data_[i] * other_data[j]);

                        // ->data_[i * other_data.size() + j] = data_[i] * other_data[j];
                    }
                }

                return result;
            }

            template<index_t Rank2, index_t Dim2>
            std::unique_ptr<tensor_base<Rank + Rank2 - 2, Dim, real_t>> dot(const tensor_base<Rank2, Dim2, real_t>& other) const {
                if (this->data_size() != other.data_size()) {
                    throw std::invalid_argument("Tensor sizes do not match for dot product.");
                }
            
                const std::vector<real_t>& other_data = other.get_data();
                index_t total_size = this->data_size();
            
                real_t global_dot = 0.0;
            
                #ifdef USING_OPENMP
                #pragma omp parallel for reduction(+:global_dot)
                #endif
                for (index_t i = 0; i < total_size; ++i) {
                    global_dot += this->data_[i] * other_data[i];
                }
            
                auto result_tensor = std::make_unique<tensor_omp<Rank + Rank2 - 2, Dim, real_t>>();
                result_tensor->init({global_dot});
            
                return result_tensor;
            }
            
            void print() const override {
                for (const auto& val : data_) {
                    std::cout << val << " ";
                }
                std::cout << std::endl;
            }
    };

} // namespace data_structure

#endif