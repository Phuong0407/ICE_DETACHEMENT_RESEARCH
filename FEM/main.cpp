#include <iostream>
#include <data_structure/matrix/matrix.hpp>

int main() {
    using namespace data_structure;

    // Define two 3x3 matrices for testing
    matrix_base<3, 3> A{{1.0, 2.0, 3.0,
                         4.0, 5.0, 6.0,
                         7.0, 8.0, 9.0}};

    matrix_base<3, 3> B{{9.0, 8.0, 7.0,
                         6.0, 5.0, 4.0,
                         3.0, 2.0, 1.0}};

    // Perform A - B
    matrix_base<3, 3> C = A + B;

    // Print the results
    std::cout << "Matrix A:\n";
    for (unsigned int i = 0; i < A.row_size(); ++i) {
        for (unsigned int j = 0; j < A.col_size(); ++j) {
            std::cout << A(i, j) << " ";
        }
        std::cout << std::endl;
    }

    std::cout << "\nMatrix B:\n";
    for (unsigned int i = 0; i < B.row_size(); ++i) {
        for (unsigned int j = 0; j < B.col_size(); ++j) {
            std::cout << B(i, j) << " ";
        }
        std::cout << std::endl;
    }

    std::cout << "\nMatrix C (A * B):\n";
    for (unsigned int i = 0; i < C.row_size(); ++i) {
        for (unsigned int j = 0; j < C.col_size(); ++j) {
            std::cout << C(i, j) << " ";
        }
        std::cout << std::endl;
    }

    return 0;
}
