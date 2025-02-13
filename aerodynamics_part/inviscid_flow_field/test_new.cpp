#include "./grid_generation/grid_generation.hpp"

#include <stdio.h>
#include <stdlib.h>
#include <stddef.h>
#include <iostream>
#include <fstream>
#include <vector>
#include <array>
#include <cstdlib>

void visualize_mesh_via_gnu_plot(const std::vector<double>& x, 
                                 const std::vector<double>& y, 
                                 const std::vector<std::array<unsigned int, 3>>& elements) 
{
    if (x.size() != y.size()) {
        std::cerr << "Error: x and y coordinate vectors must have the same size!\n";
        return;
    }

    std::ofstream nodeFile("grid_points.dat");
    std::ofstream connFile("grid_connections.dat");

    if (!nodeFile || !connFile) {
        std::cerr << "Error: Unable to open output files!\n";
        return;
    }

    for (size_t i = 0; i < x.size(); i++) {
        nodeFile << x[i] << " " << y[i] << "\n";
    }

    for (const auto& elem : elements) {
        unsigned int idx1 = elem[0];
        unsigned int idx2 = elem[1];
        unsigned int idx3 = elem[2];

        if (idx1 >= x.size() || idx2 >= x.size() || idx3 >= x.size()) {
            std::cerr << "Error: Element index (" << idx1 << ", " << idx2 << ", " << idx3
                      << ") out of range (valid range: 0 - " << x.size() - 1 << ")!\n";
            return;
        }

        connFile << x[idx1] << " " << y[idx1] << "\n" << x[idx2] << " " << y[idx2] << "\n\n";
        connFile << x[idx2] << " " << y[idx2] << "\n" << x[idx3] << " " << y[idx3] << "\n\n";
        connFile << x[idx3] << " " << y[idx3] << "\n" << x[idx1] << " " << y[idx1] << "\n\n";
    }

    std::ofstream scriptFile("plot_grid.gnu");
    if (!scriptFile) {
        std::cerr << "Error: Unable to create Gnuplot script!\n";
        return;
    }

    scriptFile << "set title 'Computational Mesh'\n"
               << "set xlabel 'X-axis'\n"
               << "set ylabel 'Y-axis'\n"
               << "set grid\n"
               << "set key off\n"
               << "plot 'grid_points.dat' with points pt 7 ps 0.1 lc rgb 'blue', \\\n"
               << "     'grid_connections.dat' with lines lc rgb 'black'\n";
    nodeFile.close();
    connFile.close();
    scriptFile.close();

    system("/usr/bin/gnuplot -p plot_grid.gnu");
}

int main() {
    std::vector<double> x, y;
    std::vector<double> x2, y2;
    std::vector<std::array<unsigned int, 3>> connection;
    double L1 = 4.0;
    double L2 = 2.0;
    double H1 = 3.0;
    double H2 = 3.0;
    double H3 = 3.0;
    double a = (L1 - L2) / H2;

    unsigned int M1 = 20;
    unsigned int M2 = 20;
    unsigned int M3 = 10;
    unsigned int N = 20;

    auto xAB = [](double h, double H1, double H2) -> double { return H1 + H2 * h; };
    auto yAB = [](double h, double a, double H2) -> double { return a * H2 * h; };
    auto xCD = [](double h, double H1, double H2) -> double { return H1 + H2 * h; };
    auto yCD = [](double h, double L) -> double { return L; };
    auto dxAB = [](double h, double H1, double H2) -> double { return H2 - H1; };
    auto dyAB = [](double h, double a, double H2) -> double { return a * H2; };
    auto dxCD = [](double h, double H2) -> double { return H2; };
    auto dyCD = [](double h) -> double { return 0.0; };

    auto xAB_wrapped = [xAB, H1, H2](double h) -> double { return xAB(h, H1, H2); };
    auto yAB_wrapped = [yAB, a, H2](double h) { return yAB(h, a, H2); };
    auto xCD_wrapped = [xCD, H1, H2](double h) -> double { return xCD(h, H1, H2); };
    auto yCD_wrapped = [L1, yCD](double h) { return yCD(h, L1); };
    auto dxAB_wrapped = [dxAB, H1, H2](double h) { return dxAB(h, H1, H2); };
    auto dyAB_wrapped = [dyAB, a, H2](double h) { return dyAB(h, a, H2); };
    auto dxCD_wrapped = [dxCD, H2](double h) { return dxCD(h, H2); };

    grid_generation::generate_grid_by_hermite_interpolation(N, M2, 3.0, 3.0, 0.0, 0.0, xAB_wrapped, yAB_wrapped, xCD_wrapped, yCD_wrapped, dxAB_wrapped, dyAB_wrapped, dxCD_wrapped, dyCD, x2, y2);
    grid_generation::generate_full_grid(H1, H2, H3, N, M1, M2, M3, 1.8, 0.1, 2.0, x2, y2, x, y);
    grid_generation::generate_grid_connection(N, M1 + M2 + M3, connection);
    visualize_mesh_via_gnu_plot(x, y, connection);
}