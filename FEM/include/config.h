#ifndef CONFIG_H
#define CONFIG_H

#include <iostream>
#include <thread>

#ifdef USING_CUDA
    #include <cuda_runtime.h>
#endif

#ifdef _OPENMP
    #define USING_OPENMP 1
    #include <omp.h>
#endif

#if defined(__has_include)
    #if __has_include(<cblas.h>)
        #define USING_BLAS 1
        #include <cblas.h>
    #endif
#endif

#if !defined(USING_BLAS) && !defined(USING_OPENMP)
    #warning "Neither BLAS nor OpenMP detected. Using sequential execution."
#endif

inline int get_cpu_threads() {
#ifdef USING_OPENMP
    return omp_get_max_threads();
#else
    return std::thread::hardware_concurrency();
#endif
}

inline int get_gpu_count() {
#ifdef USING_CUDA
    int device_count = 0;
    cudaGetDeviceCount(&device_count);
    return device_count;
#else
    return 0;
#endif
}

enum class running_device_t {
    CPU,
    GPU
};

inline running_device_t get_running_device() {
    return (get_gpu_count() > 0) ? running_device_t::GPU : running_device_t::CPU;
}

inline const int CPU_THREADS = get_cpu_threads();
inline const bool HAS_GPU = (get_gpu_count() > 0);
inline const running_device_t RUNNING_DEVICE = get_running_device();

inline void print_system_info() {
    std::cout << "====================================\n";
    std::cout << "System Hardware Configuration:\n";
    std::cout << "------------------------------------\n";
    std::cout << "CPU Threads Available: " << CPU_THREADS << "\n";

    if (HAS_GPU) {
#ifdef USING_CUDA
        int device;
        cudaGetDevice(&device);

        cudaDeviceProp prop;
        cudaGetDeviceProperties(&prop, device);

        std::cout << "GPU Detected: " << prop.name << "\n";
        std::cout << "GPU Compute Capability: " << prop.major << "." << prop.minor << "\n";
        std::cout << "Total GPU Memory: " << prop.totalGlobalMem / (1024 * 1024) << " MB\n";
#endif
    } else {
        std::cout << "No GPU Detected. Running on CPU.\n";
    }

    std::cout << "====================================\n";
}

template<typename real_t>
constexpr real_t ETOL = 
#ifdef USE_HYPER_PRECISION
    static_cast<real_t>(1e-15);
#elif defined(USE_ULTRA_PRECISION)
    static_cast<real_t>(1e-12);
#elif defined(USE_AVERAGE_PRECISION)
    static_cast<real_t>(1e-6);
#else
    static_cast<real_t>(1e-12);
#endif
;

#endif // CONFIG_H