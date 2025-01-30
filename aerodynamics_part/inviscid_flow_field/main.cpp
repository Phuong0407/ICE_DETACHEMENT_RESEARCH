#include "geometry.h"
#include <iostream>

int main() {
    Geometry G;
    G.initPhysicalDimensions(3,3,3,4,2);
    G.initComputationalDimensions(10,10,10,10);
    G.initCoordsArrSize();
    G.initGridStep();
    // G.generate_boundary_node();
    // G.transfinite_interpolation();
    // G.generate_stretching_grid(0.5, 0.5);
    G.generateAlgebraicGrid();
    G.generateEllipicGrid();
    G.generateGridConnection();

    Point2DArr x = G.get_x();
    Point2DArr y = G.get_y();
    ElemConnArr c = G.getElementConnectionData();

    if (c.size() != 0) {
        for (const auto & d : c) {
            std::cout << "[" << d[0] << ", " << d[1] << ", " << d[2] << "]" << std::endl;
        }
    }

    if (c.size() != 0) {
        for (const auto & d : c) {
            std::size_t i = d[0] / 31, j = d[0] % 31;
            std::cout << "[(" << x[i][j] << ", " << y[i][j] << "),";
            i = d[1] / 31, j = d[1] % 31;
            std::cout << "(" << x[i][j] << ", " << y[i][j] << "),";
            i = d[2] / 31, j = d[2] % 31;
            std::cout << "(" << x[i][j] << ", " << y[i][j] << ")]," << std::endl;
        }
    }

}