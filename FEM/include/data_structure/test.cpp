// #include <iostream>
// #include <vector>
// #include <random>
// #include <chrono>  // For measuring execution time
// #include <omp.h>   // OpenMP

// #define USING_OPENMP
// #include "tensor_omp.hpp"

// int main() {
//     using namespace data_structure;

//     // Set OpenMP to use 8 threads
//     omp_set_num_threads(8);

//     constexpr std::size_t SIZE = 1'000'000; // Large dataset size

//     // Start measuring total execution time
//     auto total_start = std::chrono::high_resolution_clock::now();

//     for (unsigned int i = 0; i < 100; ++i) {
//         // Create large 1D tensors (Rank = 1, Dim = SIZE)
//         tensor_omp<1, SIZE, double, unsigned int> tensor1, tensor2;

//         // Generate large random data with parallelization
//         std::vector<double> data1(SIZE);
//         std::vector<double> data2(SIZE);

//         #pragma omp parallel
//         {
//             std::random_device rd;
//             std::mt19937 gen(rd());
//             std::uniform_real_distribution<double> dist(0.0, 100.0);

//             #pragma omp for
//             for (std::size_t i = 0; i < SIZE; ++i) {
//                 data1[i] = dist(gen);
//                 data2[i] = dist(gen);
//             }
//         }

//         // Initialize tensors
//         tensor1.init(data1);
//         tensor2.init(data2);

//         // Compute dot product
//         auto result = tensor1.dot(tensor2);

//             // Print the dot product result
//         if (result) {
//             std::cout << "Dot product result: ";
//             result->print();
//         } else {
//             std::cerr << "Dot product computation failed.\n";
//         }
//     }
//     // Stop measuring total execution time
//     auto total_end = std::chrono::high_resolution_clock::now();
//     std::chrono::duration<double> total_time = total_end - total_start;

//     // Print total execution time
//     std::cout << "Total execution time: " << total_time.count() << " seconds\n";



//     return 0;
// }


#include <iostream>
#include "tensor_cu.hpp"

int main() {
    using namespace data_structure;

    constexpr std::size_t SIZE1 = 3; // Vector size
    constexpr std::size_t SIZE2 = 4;

    tensor_cu<1, SIZE1, double> tensor1;
    tensor_cu<1, SIZE2, double> tensor2;

    std::vector<double> data1 = {1.0, 2.0, 3.0};  // 3x1 vector
    std::vector<double> data2 = {4.0, 5.0, 6.0, 7.0};  // 1x4 vector

    tensor1.init(data1);
    tensor2.init(data2);

    std::cout << "Tensors initialized on GPU.\n";

    // Compute outer product
    auto result = outer(tensor1, tensor2);

    // Print result
    std::cout << "Outer product result:\n";
    result->print();

    return 0;
}
