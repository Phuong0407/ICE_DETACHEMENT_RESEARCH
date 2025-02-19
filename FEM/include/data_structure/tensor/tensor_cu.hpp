#ifndef DATA_STRUCTURE_TENSOR_CU_HPP
#define DATA_STRUCTURE_TENSOR_CU_HPP

#include "tensor_base.hpp"

#include <iostream>
#include <vector>
#include <cuda_runtime.h>

namespace data_structure {

    template<unsigned int Rank, unsigned int Dim, typename real_t = double, typename index_t = unsigned int>
    class tensor_cu : public tensor_base<Rank, Dim, real_t, index_t> {
        private:
            index_t size_;
            real_t* device_data_;

        public:
            tensor_cu() {
                size_ = 1;
                for (unsigned int i = 0; i < Rank; ++i) {
                    size_ *= Dim;
                }
                cudaMalloc(&device_data_, size_ * sizeof(real_t));   
            }

            ~tensor_cu() override {
                cudaFree(device_data_);
            }

            void init() override {
                cudaMemset(device_data_, 0, size_ * sizeof(real_t));
            }

            void init(const real_t* raw_data, index_t size) override {
                if (size != size_) throw std::invalid_argument("Size mismatch.");
                cudaMemcpy(device_data_, raw_data, size_ * sizeof(real_t), cudaMemcpyHostToDevice);
            }

            void init(const std::vector<real_t>& vec) override {
                if (vec.size() != size_) throw std::invalid_argument("Size mismatch.");
                cudaMemcpy(device_data_, vec.data(), size_ * sizeof(real_t), cudaMemcpyHostToDevice);
            }

            std::vector<real_t> get_data_host() const {
                std::vector<real_t> host_data(size_);
                cudaMemcpy(host_data.data(), device_data_, size_ * sizeof(real_t), cudaMemcpyDeviceToHost);
                return host_data;
            }

            const std::vector<real_t>& get_data() const override {
                return get_data_host();
            }

            unsigned int data_size() const override { 
                return static_cast<unsigned int>(size_); 
            }

            void set_data(index_t idx, real_t val) override {
                cudaMemcpy(device_data_ + idx, &val, sizeof(real_t), cudaMemcpyHostToDevice);
            }

            void print() const override {
                std::vector<real_t> host_data = get_data_host();
                for (const auto& val : host_data) {
                    std::cout << val << " ";
                }
                std::cout << std::endl;
            }
    };

} // namespace data_structure

#endif