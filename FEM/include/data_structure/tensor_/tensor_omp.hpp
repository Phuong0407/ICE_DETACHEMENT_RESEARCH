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

            const std::array<index_t, Rank>& get_shape() const override {
                return shape_;
            }

            void resize(const std::array<index_t, Rank>& new_shape) {
                shape_ = new_shape;
                index_t total_size = 1;
                for (index_t s : shape_) total_size *= s;
                data_.resize(total_size, 0.0);
            }            

        public:
            using tsomp_ptr = std::unique_ptr<tensor_omp>;
            
            tensor_omp(const tensor_omp& other) {
                data_ = other.data_;
            }
    
            tensor_omp& operator=(const tensor_omp& other) {
                if (this != &other) {
                    data_ = other.data_;
                }
                return *this;
            }
    
            tensor_omp(tensor_omp&& other) noexcept {
                data_ = std::move(other.data_);
            }
    
            tensor_omp& operator=(tensor_omp&& other) noexcept {
                if (this != &other) {
                    data_ = std::move(other.data_);
                }
                return *this;
            }
    
            void set_data(const std::vector<real_t>& new_data) override {
                if (new_data.size() != data_.size())
                    throw std::invalid_argument("size mismatch!");
                data_ = new_data;
            }

            void set_data(index_t idx, real_t val) override {
                if (idx >= data_.size()) {
                    throw std::out_of_range("index out of bounds in set_data()");
                }
                data_[idx] = val;
            }            

            void set_data(real_t value) override {
                if (data_.size() != 1)
                    throw std::invalid_argument("set_data(scalar) can only be used for scalar tensors (size=1)!");
                
                data_[0] = value;
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

                std::vector<real_t> temp_data(list);
                #ifdef USING_OPENMP
                #pragma omp parallel for
                #endif
                for (index_t i = 0; i < data_.size(); ++i) {
                    data_[i] = temp_data[i];
                }
            }

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

            template<unsigned int Rank2, unsigned int Dim2>
            std::unique_ptr<tensor_omp<Rank + Rank2, Dim, real_t, index_t>> outer(const tensor_base<Rank2, Dim2, real_t>& other) const {

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
                    }
                }
                return result;
            }

            template<unsigned int Rank2, unsigned int Dim2>
            std::unique_ptr<tensor_omp<Rank + Rank2 - 2, Dim, real_t, index_t>> dot(const tensor_base<Rank2, Dim2, real_t>& other) const {
                auto outer_tensor = this->outer(other);
                return outer_tensor->contract(Rank - 1, Rank);
            }

            typename tensor_base<Rank - 2, Dim, real_t, index_t>::ts_ptr contract(index_t i, index_t j) override {
                static_assert(Rank > 2, "Cannot contract a tensor with Rank <= 2");
            
                if (i >= Rank || j >= Rank || i == j) {
                    throw std::invalid_argument("Invalid contraction indices.");
                }
            
                std::array<index_t, Rank> shape = this->get_shape();
                if (shape[i] != shape[j]) {
                    throw std::invalid_argument("Contracted dimensions must be equal.");
                }
            
                // Compute the new shape after contraction
                std::array<index_t, Rank - 2> new_shape;
                index_t new_idx = 0;
                for (index_t idx = 0; idx < Rank; ++idx) {
                    if (idx != i && idx != j) {
                        new_shape[new_idx++] = shape[idx];
                    }
                }
            
                // Create the result tensor with the new shape
                auto result = std::make_unique<tensor_omp<Rank - 2, Dim, real_t, index_t>>();
                result->resize(new_shape);
            
                const std::vector<real_t>& A = this->get_data();
                std::vector<real_t>& C = result->get_data();
            
                index_t result_size = result->data_size();
            
                // Perform contraction
                #ifdef USING_OPENMP
                #pragma omp parallel for
                #endif
                for (index_t idx = 0; idx < result_size; ++idx) {
                    real_t sum = 0.0;
                    for (index_t k = 0; k < shape[i]; ++k) {
                        std::array<index_t, Rank> full_idx = this->unravel_index(idx, new_shape);
                        full_idx[i] = k;
                        full_idx[j] = k;
                        index_t idx_A = this->compute_index(full_idx);
                        sum += A[idx_A];
                    }
                    C[idx] = sum;
                }
            
                return result;
            }

            template<typename ...indices_t>
            real_t operator()(indices_t ... indices) const {
                static_assert(sizeof...(indices) == Rank, "incorrect number of indices to access the tensor");
                std::array<index_t, Rank> idx_array = {static_cast<index_t>(indices)...};
                return data_[compute_index(idx_array)];
            }
            template<typename ...indices_t>
            real_t& operator()(indices_t ... indices) {
                static_assert(sizeof...(indices) == Rank, "incorrect number of indices to access the tensor");
                std::array<index_t, Rank> idx_array = {static_cast<index_t>(indices)...};
                return data_[compute_index(idx_array)];
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