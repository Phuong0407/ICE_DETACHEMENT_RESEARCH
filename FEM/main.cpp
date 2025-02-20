
// #include <iostream>
// #include <data_structure/tensor/tensor.hpp>

// int main() {
//     using namespace data_structure;

//     std::cout << "====================================" << std::endl;
//     std::cout << "         Tensor API Test            " << std::endl;
//     std::cout << "====================================" << std::endl;

//     // ✅ Create a CPU tensor (Rank=2, Dim=3)
//     tensor<2, 3, double> tensor_cpu1, tensor_cpu2;
    
//     // ✅ Initialize tensors
//     tensor_cpu1.init({1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0});
//     tensor_cpu2.init({7.0, 8.0, 9.0, 10.0, 11.0, 12.0, 13.0, 14.0, 15.0});

//     std::cout << "Tensor 1: ";
//     tensor_cpu1.print();

//     std::cout << "Tensor 2: ";
//     tensor_cpu2.print();

//     // ✅ Perform Addition
//     auto sum_tensor = tensor_cpu1 + tensor_cpu2;
//     std::cout << "Sum Tensor: ";
//     sum_tensor.print();

//     // ✅ Perform Subtraction
//     auto diff_tensor = tensor_cpu1 - tensor_cpu2;
//     std::cout << "Difference Tensor: ";
//     diff_tensor.print();

//     // ✅ Perform Element-wise Multiplication
//     auto prod_tensor = tensor_cpu1 * tensor_cpu2;
//     std::cout << "Product Tensor: ";
//     prod_tensor.print();

//     // ✅ Compute Dot Product
//     auto cpu_dot = tensor_cpu1.dot(tensor_cpu2);
//     std::cout << "Dot Product Result: ";
//     cpu_dot.print();

//     // ✅ Compute Outer Product
//     auto cpu_outer = tensor_cpu1.outer(tensor_cpu2);
//     std::cout << "Outer Product Result: ";
//     cpu_outer.print();

// #ifdef USING_CUDA
//     std::cout << "\n====================================" << std::endl;
//     std::cout << "         GPU Tensor Test            " << std::endl;
//     std::cout << "====================================" << std::endl;

//     // ✅ Create GPU tensors
//     tensor<2, 3, double> tensor_gpu1(running_device_t::GPU), tensor_gpu2(running_device_t::GPU);

//     // ✅ Initialize GPU tensors
//     tensor_gpu1.init({1.0, 2.0, 3.0, 4.0, 5.0, 6.0});
//     tensor_gpu2.init({7.0, 8.0, 9.0, 10.0, 11.0, 12.0});

//     std::cout << "GPU Tensor 1: ";
//     tensor_gpu1.print();

//     std::cout << "GPU Tensor 2: ";
//     tensor_gpu2.print();

//     // ✅ Perform GPU Addition
//     auto sum_gpu = tensor_gpu1 + tensor_gpu2;
//     std::cout << "GPU Sum Tensor: ";
//     sum_gpu.print();

//     // ✅ Compute GPU Dot Product
//     auto gpu_dot = tensor_gpu1.dot(tensor_gpu2);
//     std::cout << "GPU Dot Product Result: ";
//     gpu_dot.print();
// #endif

//     return 0;
// }


#include <iostream>
#include <data_structure/tensor/tensor.hpp>

int main() {
    std::cout << "====================================\n";
    std::cout << "         Tensor API Test            \n";
    std::cout << "====================================\n";

    using namespace data_structure;

    tensor<2, 3, double> tensor1(running_device_t::CPU);
    tensor<2, 3, double> tensor2(running_device_t::CPU);

    tensor1.init({1, 2, 3, 4, 5, 6, 7, 8, 9});
    tensor2.init({9, 8, 7, 6, 5, 4, 3, 2, 1});

    std::cout << "Tensor 1: ";
    tensor1.print();

    std::cout << "Tensor 2: ";
    tensor2.print();

    auto sum_tensor = tensor1 + tensor2;
    std::cout << "Sum Tensor: ";
    sum_tensor.print();

    auto diff_tensor = tensor1 - tensor2;
    std::cout << "Difference Tensor: ";
    diff_tensor.print();

    auto prod_tensor = tensor1 * tensor2;
    std::cout << "Product Tensor: ";
    prod_tensor.print();

    try {
        auto dot_tensor = tensor1.dot(tensor2);
        std::cout << "Dot Product Result: ";
        dot_tensor.print();
    } catch (const std::exception& e) {
        std::cerr << "Error in dot product: " << e.what() << std::endl;
    }

    try {
        auto outer_tensor = tensor1.outer(tensor2);
        std::cout << "Outer Product Result: ";
        outer_tensor.print();
    } catch (const std::exception& e) {
        std::cerr << "Error in outer product: " << e.what() << std::endl;
    }

    std::cout << "====================================\n";
    std::cout << "         Tensor API Test Complete   \n";
    std::cout << "====================================\n";

    return 0;
}
