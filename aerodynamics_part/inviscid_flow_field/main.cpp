#include "geometry.h"
#include <iostream>
#include <vector>
#include <fstream>
#include <cstdlib>  // For system() function

int main() {
    Geometry G;
    G.initPhysicalDimensions(3,3,3,4,2);
    G.initComputationalDimensions(5,5,5,5);
    G.initCoordsArrSize();
    G.initGridStep();
    G.generateBoundaryNode();
    // G.transfiniteInterpolation();
    // G.generateStretchingGrid(0.9, 0.8);
    G.generateAlgebraicGrid();
    G.generateEllipicGrid();
    G.generateGridConnection();

    Point2DArr x = G.get_x();
    Point2DArr y = G.get_y();
    ElemConnArr c = G.getElementConnectionData();

    std::size_t M = G.getTotalHorizontalElements();
    std::ofstream nodeFile("grid_points.dat");
    std::ofstream connFile("grid_connections.dat");

    if (!x.empty() && nodeFile.is_open() && connFile.is_open()) {
        for (std::size_t i = 0; i < x.size(); ++i) {
            for (std::size_t j = 0; j < x[i].size(); ++j) {
                nodeFile << x[i][j] << " " << y[i][j] << "\n";
            }
            nodeFile << "\n";
        }
        nodeFile.close();

        if (!c.empty()) {
            for (const auto &d : c) {
                std::size_t i1 = d[0] / (M + 1), j1 = d[0] % (M + 1);
                std::size_t i2 = d[1] / (M + 1), j2 = d[1] % (M + 1);
                std::size_t i3 = d[2] / (M + 1), j3 = d[2] % (M + 1);

                connFile << x[i1][j1] << " " << y[i1][j1] << "\n";
                connFile << x[i2][j2] << " " << y[i2][j2] << "\n\n";
                
                connFile << x[i2][j2] << " " << y[i2][j2] << "\n";
                connFile << x[i3][j3] << " " << y[i3][j3] << "\n\n";
                
                connFile << x[i3][j3] << " " << y[i3][j3] << "\n";
                connFile << x[i1][j1] << " " << y[i1][j1] << "\n\n";
            }
        }
        connFile.close();
    } else {
        std::cerr << "Error: Unable to write grid data files!" << std::endl;
        return 1;
    }

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