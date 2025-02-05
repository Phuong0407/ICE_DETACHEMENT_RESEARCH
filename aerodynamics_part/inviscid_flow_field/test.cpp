#include "geometry.hpp"
#include "axisymmetric_fem.hpp"
#include "csr_sparse_solver.hpp"
#include "post_processing.hpp"

#include <iostream>
#include <vector>
#include <fstream>
#include <cstdlib>

#include <chrono>

int main() {
    auto start = std::chrono::high_resolution_clock::now();

    Geometry G;
    G.initPhysicalDimensions(3,3,3,4,2);
    G.initComputationalDimensions(100, 100, 100, 100);
    G.initCoordsArrSize();
    G.initGridStep();
    G.generateBoundaryNode();
    G.generateAlgebraicGrid();
    // G.generateEllipicGrid();
    G.generateGridConnection();
    G.generateNeumannBoundaryConditions(0.0, 4.0/3.0, 0.0, -1.0);

    std::vector<double> x = G.get_x();
    std::vector<double> y = G.get_y();
    std::vector<std::vector<std::size_t>> c = G.getElementConnectionData();
    std::vector<BoundaryEdge> d = G.getNeumannBoundaryConditions();
    AxisymmetricFiniteElement fem;
    fem.computeGlobalStiffnessMatrix(x, y, c);
    fem.computeNeumannBoundaryConditions(x, y, c, d);
    fem.computeSolution();
    std::vector<double> solution = fem.getSolution();

    // std::size_t number_nodes = x.size();
    // Gradient_Recovery g;
    // g.initNumberOfNode(number_nodes);
    // g.generateNodeRecoveryPatchConnManager(c);
    // g.generateGradientRecovery(x, y, solution);
    // std::vector<double> dev_x = g.get_dev_x();
    // std::vector<double> dev_y = g.get_dev_y();
    
    // // for (std::size_t i = 0; i < x.size(); ++i) {
    // //     std::cout << "(Dev_x, Dev y)("<< x[i] << ", " << y[i] << ") = " << dev_x[i] << ", " << dev_y[i] << std::endl;
    // // }

    // // Open a file for writing
    // std::ofstream outFile("output.txt");

    // // Check if the file is open
    // if (!outFile.is_open()) {
    //     std::cerr << "Error: Unable to open file for writing!" << std::endl;
    //     return 1;
    // }

    // // Write the data to the file
    // for (std::size_t i = 0; i < x.size(); ++i) {
    //     outFile << "(Dev_x, Dev_y)(" << x[i] << ", " << y[i] << ") = "
    //             << dev_x[i] << ", " << dev_y[i] << std::endl;
    // }

    // // Close the file
    // outFile.close();

    // std::ofstream nodeFile("grid_points.dat");
    // std::ofstream connFile("grid_connections.dat");

    // if (!nodeFile.is_open() || !connFile.is_open()) {
    //     std::cerr << "Error: Unable to open output files!" << std::endl;
    //     return 1;
    // }

    // if (!x.empty()) {
    //     for (std::size_t i = 0; i < x.size(); ++i) {
    //         nodeFile << x[i] << " " << y[i] << "\n";
    //     }
    //     nodeFile << "\n";
    // } else {
    //     std::cerr << "Warning: No grid points found!" << std::endl;
    // }
    // nodeFile.close();

    // if (!c.empty()) {
    //     for (const auto &d : c) {
    //         std::size_t idx1 = d[0];
    //         std::size_t idx2 = d[1];
    //         std::size_t idx3 = d[2];

    //         connFile << x[idx1] << " " << y[idx1] << "\n";
    //         connFile << x[idx2] << " " << y[idx2] << "\n\n";
    //         connFile << x[idx2] << " " << y[idx2] << "\n";
    //         connFile << x[idx3] << " " << y[idx3] << "\n\n";

    //         connFile << x[idx3] << " " << y[idx3] << "\n";
    //         connFile << x[idx1] << " " << y[idx1] << "\n\n";
    //     }
    // } else {
    //     std::cerr << "Warning: No connectivity data found!" << std::endl;
    // }
    // connFile.close();

    // std::ofstream scriptFile("plot_grid.gnu");
    // scriptFile << "set title 'Computational Mesh'\n";
    // scriptFile << "set xlabel 'X-axis'\n";
    // scriptFile << "set ylabel 'Y-axis'\n";
    // scriptFile << "set grid\n";
    // scriptFile << "set key off\n";
    // scriptFile << "plot 'grid_points.dat' with points pt 7 ps 1 lc rgb 'blue', \\\n";
    // scriptFile << "     'grid_connections.dat' with lines lc rgb 'black'\n";
    // scriptFile.close();

    // system("gnuplot -p plot_grid.gnu");
    // return 0;

    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = end - start;

    std::cout << "Execution time: " << elapsed.count() << " seconds\n";

}