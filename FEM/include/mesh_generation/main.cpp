#include "point.hpp"
#include "edge.hpp"
#include "vector.hpp"
#include "line.hpp"
#include "polygon.hpp"
#include "advancing_front.hpp"
#include "front_advancing.hpp"

#include <iostream>
#include <fstream>
#include <cstdlib>

#define _spdim 2
#define double_t double
#define index_t unsigned int

using namespace mesh_generation;

void save_and_plot_discretized_boundary(const std::vector<edge<_spdim, double_t>>& discretized_boundary) {
    std::string data_filename = "boundary.dat";

    // Step 1: Save the discretized boundary to a file
    std::ofstream data_file(data_filename);
    if (!data_file.is_open()) {
        std::cerr << "Error: Could not open file " << data_filename << " for writing.\n";
        return;
    }

    for (const auto& edge : discretized_boundary) {
        point<_spdim, double_t> A = edge.first();
        point<_spdim, double_t> B = edge.second();

        data_file << A(0) << " " << A(1) << "\n";  // First point
        data_file << B(0) << " " << B(1) << "\n";  // Second point
        data_file << "\n";  // Blank line to separate edges
    }

    data_file.close();
    std::cout << "Discretized boundary saved to " << data_filename << std::endl;

    // Step 2: Generate Gnuplot script
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
               << "plot 'boundary.dat' with points pt 7 ps 1 lc rgb 'blue', \\\n";
    data_file.close();
    scriptFile.close();
    // Step 3: Run Gnuplot
    system("/usr/bin/gnuplot -p plot_grid.gnu");  
}


int main() {
    // mesh_generation::point<_spdim> p(1.0, 2.0);
    // mesh_generation::point<_spdim> q(2.0, 4.0);
    // mesh_generation::point<_spdim> r(0.0, 0.0);
    // mesh_generation::edge<_spdim> e(p, q);
    // mesh_generation::point<_spdim> f = r.project_to_edge(e);
    // f.print();

    mesh_generation::polygon<_spdim, double> polygon;
    
    polygon.add_vertex(0.0, 0.0);
    polygon.add_vertex(3.0, 0.0);
    polygon.add_vertex(6.0, 2.0);
    polygon.add_vertex(9.0, 2.0);
    polygon.add_vertex(9.0, 4.0);
    polygon.add_vertex(0.0, 4.0);

    std::vector<double_t> xc = {0.1, 0.1};
    std::vector<double_t> D = {0.5, 0.5};
    std::vector<double_t> delta = {0.1, 0.1};
    std::vector<mesh_generation::point<2, double>> _S = {
        mesh_generation::point<2, double>(3.5, -1.0),
        mesh_generation::point<2, double>(6.5, 1.0)
    };

    mesh_generation::front_advancing<_spdim, double_t, index_t> mesh(polygon);

    std::vector<index_t> idx_edge_to_stretch = {1, 1, 1, 1, 1, 1};
    std::vector<index_t> n_per_edge = {10, 10, 10, 10, 20, 20};
    double q = 3.0;
    std::vector<double_t> Q = {q, q, q, q, q, q};

    mesh.initialize_background_mesh(xc, D, delta, _S);
    mesh.generate_mesh(idx_edge_to_stretch, n_per_edge, Q);
    mesh.print_discretized_boundary();

    std::vector<edge<_spdim, double_t>> boundary_mesh = mesh.get_discretized_boundary();
    save_and_plot_discretized_boundary(boundary_mesh);
    std::cout << boundary_mesh.size();
    return 0;
}