#ifndef GRID_GENERATION_H
#define GRID_GENERATION_H

#include "grid.hpp"

#include <vector>
#include <array>
#include <cmath>

class grid_generation {
public:
    static void grid_generation::generate_boundary(double L1, double L2, double H1, double H2, double H3, unsigned int M1, unsigned int M2, unsigned int M3, std::vector<double> &bottom_x, std::vector<double> &bottom_y, std::vector<double> &top_y);
    static void generate_streching_grid(double alpha, double eta0, double eta1, double eta2, unsigned int N, const std::vector<double> bottom_x, const std::vector<double> bottom_y, const std::vector<double> top_y, std::vector<double> &x, std::vector<double> &y);
};

void grid_generation::generate_boundary(double L1, double L2, double H1, double H2, double H3, unsigned int M1, unsigned int M2, unsigned int M3, std::vector<double> &bottom_x, std::vector<double> &bottom_y, std::vector<double> &top_y) {
    unsigned int M = M1 + M2 + M3 + 1;
    double dx1 = H1 / M1, dx2 = H2 / M2, dx3 = H3 / M3;

    top_y.resize(M), bottom_x.resize(M), bottom_y.resize(M);
    for (size_t i = 0; i <= M1; ++i) {
        double x = dx1 * i;
        bottom_x[i] = dx1 * x;
        bottom_y[i] = 0.0;
        top_y[i] = L1;
    }
    double a = (L1 - L2) / H2;
    for (size_t i = M1 + 1; i <= M1 + M2; ++i) {
        double x = H1 + dx2 * (i - M1);
        double y = a * dx2 * (i - M1);
        bottom_x[i] = x;
        bottom_y[i] = y;
        top_y[i] = L1;
    }
    for (size_t i = M1 + M2 + 1; i <= M; ++i) {
        double x = H1 + H2 + dx3 * (i - M1 - M2);
        bottom_x[i] = x;
        bottom_y[i] = L1 - L2;
        top_y[i] = L1;
    }
}



void grid_generation::generate_streching_grid(double alpha, double eta0, double eta1, double eta2, unsigned int N, const std::vector<double> bottom_x, const std::vector<double> bottom_y, const std::vector<double> top_y, std::vector<double> &x, std::vector<double> &y) {
    unsigned int M = bottom_x.size();
    x.resize((N + 1) * M), y.resize((N + 1) * M);
    for (unsigned int i = 0; i < N; ++i) {
        double eta = i * 1.0 / N;
        for (unsigned int j = 0; j <= M; ++j) {
            unsigned int ij = i * M + j;
            x[ij] = bottom_x[j];
            if (eta <= eta1)
                y[ij] = (top_y[j] - bottom_y[j]) * eta1 * (std::exp(alpha * eta / eta1) - 1.0) / (std::exp(alpha) - 1) + bottom_y[j];
            else
                y[ij] = (top_y[j] - bottom_y[j]) * (1.0 - (1.0 - eta1) * (std::exp(alpha * (1 - eta) / (1 - eta1)) - 1.0) / (std::exp(alpha) - 1)) + bottom_y[j];
        }
    }
}


#endif