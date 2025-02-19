#ifndef DATA_STRUCTURE_TENSOR_HPP
#define DATA_STRUCTURE_TENSOR_HPP

#include <memory>
#include <iostream>
#include "tensor_omp.hpp"
#include "config.h"

#ifdef USING_CUDA
#include "tensor_cu.hpp"
#endif

namespace data_structure {

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
        #ifdef USING_CUDA
            else if (backend_ == running_device_t::GPU) {
                tensor_impl_ = std::make_unique<tensor_cu<Rank, Dim, real_t, index_t>>();
            }
        #endif
            else
                throw std::invalid_argument("unsupported backend selected.");
        }

        tensor(const tensor& other) : backend_(other.backend_) {
            if (other.tensor_impl_) {
                tensor_impl_ = std::make_unique<tensor_omp<Rank, Dim, real_t, index_t>>(*other.tensor_impl_);
            }
        }

        tensor& operator=(const tensor& other) {
            if (this != &other) {
                backend_ = other.backend_;
                if (other.tensor_impl_) {
                    tensor_impl_ = std::make_unique<tensor_omp<Rank, Dim, real_t, index_t>>(*other.tensor_impl_);
                }
            }
            return *this;
        }

        tensor(tensor&& other) noexcept : backend_(other.backend_), tensor_impl_(std::move(other.tensor_impl_)) {}

        tensor& operator=(tensor&& other) noexcept {
            if (this != &other) {
                backend_ = other.backend_;
                tensor_impl_ = std::move(other.tensor_impl_);
            }
            return *this;
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
            result.tensor_impl_ = std::make_unique<tensor_omp<Rank, Dim, real_t, index_t>>(*tensor_impl_ + *other.tensor_impl_);
            return result;
        }

        tensor operator-(const tensor& other) const {
            tensor result(backend_);
            result.tensor_impl_ = std::make_unique<tensor_omp<Rank, Dim, real_t, index_t>>(*tensor_impl_ - *other.tensor_impl_);
            return result;
        }

        tensor operator*(const tensor& other) const {
            tensor result(backend_);
            result.tensor_impl_ = std::make_unique<tensor_omp<Rank, Dim, real_t, index_t>>(*tensor_impl_ * *other.tensor_impl_);
            return result;
        }

        template<unsigned int Rank2, unsigned int Dim2>
        tensor<Rank + Rank2 - 2, Dim, real_t, index_t> dot(const tensor<Rank2, Dim2, real_t, index_t>& other) const {
            tensor<Rank + Rank2 - 2, Dim, real_t, index_t> result(backend_);
            result.tensor_impl_ = std::make_unique<tensor_omp<Rank + Rank2 - 2, Dim, real_t, index_t>>(*tensor_impl_->dot(*other.tensor_impl_));
            return result;
        }

        template<unsigned int Rank2, unsigned int Dim2>
        tensor<Rank + Rank2, Dim, real_t, index_t> outer(const tensor<Rank2, Dim2, real_t, index_t>& other) const {
            tensor<Rank + Rank2, Dim, real_t, index_t> result(backend_);
            result.tensor_impl_ = std::make_unique<tensor_omp<Rank + Rank2, Dim, real_t, index_t>>(*tensor_impl_->outer(*other.tensor_impl_));
            return result;
        }

        std::vector<real_t> get_data() const {
            return tensor_impl_->get_data();
        }

        void set_data(const tensor_base<Rank, Dim, real_t, index_t>& new_tensor) {
            tensor_impl_ = std::make_unique<tensor_omp<Rank, Dim, real_t, index_t>>(new_tensor);
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