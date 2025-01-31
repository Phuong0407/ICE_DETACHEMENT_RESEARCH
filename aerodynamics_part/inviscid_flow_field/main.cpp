#include "geometry.h"
#include <iostream>
#include <vector>
#include <fstream>
#include <cstdlib>

int main() {
    Geometry G;
    G.initPhysicalDimensions(3,3,3,4,2);
    G.initComputationalDimensions(5,5,5,5);
    G.initCoordsArrSize();
    G.initGridStep();
    G.generateBoundaryNode();
    G.generateAlgebraicGrid();
    G.generateEllipicGrid();
    G.generateGridConnection();
    G.generateNeumannBoundaryNodeIndices();

    CoordArr x = G.get_x();
    CoordArr y = G.get_y();
    ElemConnArr c = G.getElementConnectionData();

    std::ofstream nodeFile("grid_points.dat");
    std::ofstream connFile("grid_connections.dat");

    if (!nodeFile.is_open() || !connFile.is_open()) {
        std::cerr << "Error: Unable to open output files!" << std::endl;
        return 1;
    }

    if (!x.empty()) {
        for (std::size_t i = 0; i < x.size(); ++i) {
            nodeFile << x[i] << " " << y[i] << "\n";
        }
        nodeFile << "\n";
    } else {
        std::cerr << "Warning: No grid points found!" << std::endl;
    }
    nodeFile.close();

    if (!c.empty()) {
        for (const auto &d : c) {
            std::size_t idx1 = d[0];
            std::size_t idx2 = d[1];
            std::size_t idx3 = d[2];

            connFile << x[idx1] << " " << y[idx1] << "\n";
            connFile << x[idx2] << " " << y[idx2] << "\n\n";

            connFile << x[idx2] << " " << y[idx2] << "\n";
            connFile << x[idx3] << " " << y[idx3] << "\n\n";

            connFile << x[idx3] << " " << y[idx3] << "\n";
            connFile << x[idx1] << " " << y[idx1] << "\n\n";
        }
    } else {
        std::cerr << "Warning: No connectivity data found!" << std::endl;
    }
    connFile.close();

    std::ofstream scriptFile("plot_grid.gnu");
    scriptFile << "set title 'Computational Mesh'\n";
    scriptFile << "set xlabel 'X-axis'\n";
    scriptFile << "set ylabel 'Y-axis'\n";
    scriptFile << "set grid\n";
    scriptFile << "set key off\n";
    scriptFile << "plot 'grid_points.dat' with points pt 7 ps 1 lc rgb 'blue', \\\n";
    scriptFile << "     'grid_connections.dat' with lines lc rgb 'black'\n";
    scriptFile.close();

    system("gnuplot -p plot_grid.gnu");
    return 0;
}