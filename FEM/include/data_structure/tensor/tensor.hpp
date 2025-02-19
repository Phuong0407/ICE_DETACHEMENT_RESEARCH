#ifndef DATA_STRUCTURE_TENSOR_HPP
#define DATA_STRUCTURE_TENSOR_HPP

#include <data_structure/tensor/tensor_omp.hpp>

#ifdef USE_CUDA
#include <data_structure/tensor/tensor_cu.hpp>
#endif

#include <memory>
#include <iostream>

namespace data_structure {

    enum class running_device_t {
        CPU,
        GPU
    };

    template<unsigned int Rank, unsigned int Dim, typename real_t = double, typename index_t = unsigned int>
    class tensor {
    private:
        running_device_t backend_;
        std::unique_ptr<tensor_base<Rank, Dim, real_t, index_t>> tensor_impl_;

    public:
        tensor(running_device_t backend = running_device_t::CPU) : backend_(backend) {
            if (backend_ == running_device_t::CPU) {
                tensor_impl_ = std::make_unique<tensor_omp<Rank, Dim, real_t, index_t>>();
            }
        #ifdef USE_CUDA
            else if (backend_ == running_device_t::GPU) {
                tensor_impl_ = std::make_unique<tensor_cu<Rank, Dim, real_t, index_t>>();
            }
        #endif
            else
                throw std::invalid_argument("unsupported backend selected.");
        }

        void init() {
            tensor_impl_->init();
        }

        void init(const real_t* raw_data, index_t size) {
            tensor_impl_->init(raw_data, size);
        }

        void init(const std::vector<real_t>& vec) {
            tensor_impl_->init(vec);
        }

        void init(std::initializer_list<real_t> list) {
            tensor_impl_->init(list);
        }

        tensor operator+(const tensor& other) const {
            tensor result(backend_);
            result.tensor_impl_ = (*tensor_impl_ + *other.tensor_impl_);
            return result;
        }

        tensor operator-(const tensor& other) const {
            tensor result(backend_);
            result.tensor_impl_ = (*tensor_impl_ - *other.tensor_impl_);
            return result;
        }

        tensor operator*(const tensor& other) const {
            tensor result(backend_);
            result.tensor_impl_ = (*tensor_impl_ * *other.tensor_impl_);
            return result;
        }

        template<unsigned int Rank2, unsigned int Dim2>
        tensor<Rank + Rank2 - 2, Dim, real_t, index_t> dot(const tensor<Rank2, Dim2, real_t, index_t>& other) const {
            tensor<Rank + Rank2 - 2, Dim, real_t, index_t> result(backend_);
            result.tensor_impl_ = tensor_impl_->dot(*other.tensor_impl_);
            return result;
        }

        template<unsigned int Rank2, unsigned int Dim2>
        tensor<Rank + Rank2, Dim, real_t, index_t> outer(const tensor<Rank2, Dim2, real_t, index_t>& other) const {
            tensor<Rank + Rank2, Dim, real_t, index_t> result(backend_);
            result.tensor_impl_ = tensor_impl_->outer(*other.tensor_impl_);
            return result;
        }

        std::vector<real_t> get_data() const {
            return tensor_impl_->get_data();
        }

        unsigned int data_size() const {
            return tensor_impl_->data_size();
        }

        void print() const {
            tensor_impl_->print();
        }
    };

} // namespace data_structure

#endif