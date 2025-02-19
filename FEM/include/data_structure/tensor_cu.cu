#include <cuda_runtime.h>
#include <iostream>
#include "tensor_cu.hpp"

namespace data_structure {

// ✅ Custom atomicAdd for double (for GPUs < sm_60)
__device__ double atomicAdd_double(double* address, double val) {
    unsigned long long int* address_as_ull = (unsigned long long int*)address;
    unsigned long long int old = *address_as_ull, assumed;

    do {
        assumed = old;
        old = atomicCAS(address_as_ull, assumed,
                        __double_as_longlong(val + __longlong_as_double(assumed)));
    } while (assumed != old);

    return __longlong_as_double(old);
}

// ✅ CUDA kernel for dot product computation
__global__ void dot_kernel(const double* a, const double* b, double* result, std::size_t n) {
    __shared__ double temp[1024];  // Shared memory for block reduction
    int idx = threadIdx.x + blockIdx.x * blockDim.x;
    int local_idx = threadIdx.x;

    temp[local_idx] = (idx < n) ? a[idx] * b[idx] : 0.0;
    __syncthreads();

    // Parallel reduction within block
    for (int stride = blockDim.x / 2; stride > 0; stride >>= 1) {
        if (local_idx < stride) {
            temp[local_idx] += temp[local_idx + stride];
        }
        __syncthreads();
    }

    // Store result of block reduction
    if (local_idx == 0) {
        atomicAdd_double(result, temp[0]);  // ✅ Use custom atomicAdd for double
    }
}

// ✅ Dot product function for CUDA tensors
template<unsigned int Rank, unsigned int Dim, typename real_t>
std::unique_ptr<tensor_cu<Rank + Rank - 2, Dim, real_t>> dot(
    const tensor_cu<Rank, Dim, real_t>& a,
    const tensor_cu<Rank, Dim, real_t>& b)
{
    if (a.data_size() != b.data_size()) {
        throw std::invalid_argument("Size mismatch.");
    }

    std::size_t size = a.data_size();
    real_t* device_result;
    real_t host_result = 0.0;
    
    cudaMalloc(&device_result, sizeof(real_t));
    cudaMemcpy(device_result, &host_result, sizeof(real_t), cudaMemcpyHostToDevice);

    // Launch CUDA kernel
    int threads_per_block = 1024;
    int num_blocks = (size + threads_per_block - 1) / threads_per_block;
    dot_kernel<<<num_blocks, threads_per_block>>>(a.device_data_, b.device_data_, device_result, size);

    // Copy result back to host
    cudaMemcpy(&host_result, device_result, sizeof(real_t), cudaMemcpyDeviceToHost);
    cudaFree(device_result);

    // ✅ Return result as tensor
    auto result_tensor = std::make_unique<tensor_cu<Rank + Rank - 2, Dim, real_t>>();
    result_tensor->init({host_result});
    return result_tensor;
}

// ✅ CUDA kernel for outer product computation
__global__ void outer_kernel(const double* a, const double* b, double* c, std::size_t m, std::size_t n) {
    int row = blockIdx.y * blockDim.y + threadIdx.y;
    int col = blockIdx.x * blockDim.x + threadIdx.x;

    if (row < m && col < n) {
        c[row * n + col] = a[row] * b[col];
    }
}

// ✅ Outer product function for CUDA tensors
template<unsigned int Rank1, unsigned int Dim1, unsigned int Rank2, unsigned int Dim2, typename real_t>
std::unique_ptr<tensor_cu<Rank1 + Rank2, Dim1, real_t>> outer(
    const tensor_cu<Rank1, Dim1, real_t>& a,
    const tensor_cu<Rank2, Dim2, real_t>& b)
{
    std::size_t m = a.data_size();  // Number of rows
    std::size_t n = b.data_size();  // Number of columns
    std::size_t size = m * n;       // Total elements in the outer product result

    // Allocate memory for result tensor on GPU
    real_t* device_c;
    cudaMalloc(&device_c, size * sizeof(real_t));

    // ✅ Define optimized grid and block dimensions
    dim3 threads_per_block(16, 16);
    dim3 num_blocks((n + threads_per_block.x - 1) / threads_per_block.x,
                    (m + threads_per_block.y - 1) / threads_per_block.y);

    // ✅ Launch CUDA kernel
    outer_kernel<<<num_blocks, threads_per_block>>>(a.device_data_, b.device_data_, device_c, m, n);

    // ✅ Create and initialize result tensor
    auto result_tensor = std::make_unique<tensor_cu<Rank1 + Rank2, Dim1, real_t>>();
    cudaMemcpy(result_tensor->device_data_, device_c, size * sizeof(real_t), cudaMemcpyDeviceToDevice);

    cudaFree(device_c);  // ✅ Free temporary GPU memory

    return result_tensor;
}

} // namespace data_structure
