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
    Point2DArr a = G.get_x();
    Point2DArr b = G.get_y();

    if (a.size()!=0){
        for (std::size_t i = 0; i < a.size(); ++i) {
            for (std::size_t j = 0; j < a[i].size(); ++j) {
                std::cout << "(" << a[i][j] << ", " << b[i][j] << ")" << std::endl;
            }
        }
    }
}