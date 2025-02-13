#ifndef grid_generation_h
#define grid_generation_h

#include <vector>
#include <array>
#include <cmath>
#include <functional>
#include <cstdarg>
#include <iostream>

class grid_generation {
public:
    static void generate_top_and_bottom_boundary(
        double L1,
        double L2,
        double H1,
        double H2,
        double H3,
        unsigned int M1,
        unsigned int M2,
        unsigned int M3,
        std::vector<double> &bottom_x,
        std::vector<double> &bottom_y,
        std::vector<double> &top_y
    );

    static void generate_left_and_right_boundary(
        double L1,
        double L2,
        double H,
        unsigned int N,
        std::vector<double> &left_x,
        std::vector<double> &left_y,
        std::vector<double> &right_x,
        std::vector<double> &right_y
    );

    static void generate_grid_by_transfinite_interpolation(
        const std::vector<double> &bottom_x,
        const std::vector<double> &bottom_y,
        const std::vector<double> &top_x,
        const std::vector<double> &top_y,
        const std::vector<double> &left_x,
        const std::vector<double> &left_y,
        const std::vector<double> &right_x,
        const std::vector<double> &right_y,
        std::vector<double> &x,
        std::vector<double> &y
    );

    static void generate_stretching_grid(
        double alpha,
        double eta1,
        unsigned int N,
        const std::vector<double> bottom_x,
        const std::vector<double> bottom_y,
        const std::vector<double> top_x,
        const std::vector<double> top_y,
        std::vector<double> &x,
        std::vector<double> &y
    );

    static void generate_grid_by_hermite_interpolation(
        unsigned int N,
        unsigned int M,
        double Q1,
        double Q2,
        double T1,
        double T2,
        std::function<double(double)> xAB,
        std::function<double(double)> yAB,
        std::function<double(double)> xCD,
        std::function<double(double)> yCD,
        std::function<double(double)> dxAB,
        std::function<double(double)> dyAB,
        std::function<double(double)> dxCD,
        std::function<double(double)> dyCD,
        std::vector<double> &x,
        std::vector<double> &y
    );

    static void generate_full_grid(
        double H1,
        double H2,
        double H3,
        unsigned int N,
        unsigned int M1,
        unsigned int M2,
        unsigned int M3,
        double Q1,
        double P2,
        double Q2,
        const std::vector<double> &x2,
        const std::vector<double> &y2,
        std::vector<double> &x,
        std::vector<double> &y
    );

    static void generate_grid_connection(
        unsigned int N,
        unsigned int M,
        std::vector<std::array<unsigned int, 3>> &connection
    );
};

void grid_generation::generate_top_and_bottom_boundary(
    double L1,
    double L2,
    double H1,
    double H2,
    double H3,
    unsigned int M1,
    unsigned int M2,
    unsigned int M3,
    std::vector<double> &bottom_x,
    std::vector<double> &bottom_y,
    std::vector<double> &top_y
) {
    unsigned int M = M1 + M2 + M3;
    double dx1 = H1 / M1, dx2 = H2 / M2, dx3 = H3 / M3;

    top_y.resize(M + 1), bottom_x.resize(M + 1), bottom_y.resize(M + 1);
    for (unsigned int i = 0; i <= M1; ++i) {
        bottom_x[i] = dx1 * i;
        bottom_y[i] = 0.0;
        top_y[i] = L1;
    }
    double a = (L1 - L2) / H2;
    for (unsigned int i = M1 + 1; i <= M1 + M2; ++i) {
        double x = H1 + dx2 * (i - M1);
        double y = a * dx2 * (i - M1);
        bottom_x[i] = x;
        bottom_y[i] = y;
        top_y[i] = L1;
    }
    for (unsigned int i = M1 + M2 + 1; i <= M; ++i) {
        double x = H1 + H2 + dx3 * (i - M1 - M2);
        bottom_x[i] = x;
        bottom_y[i] = L1 - L2;
        top_y[i] = L1;
    }
}



void grid_generation::generate_left_and_right_boundary(
    double L1,
    double L2,
    double H,
    unsigned int N,
    std::vector<double> &left_x,
    std::vector<double> &left_y,
    std::vector<double> &right_x,
    std::vector<double> &right_y
) {
    left_x.resize(N + 1), left_y.resize(N + 1);
    right_x.resize(N + 1), right_y.resize(N + 1);
    double dy1 = L1 / N, dy2 = L2 / N, dL = L1 - L2;
    for (unsigned int i = 0; i <= N; ++i) {
        left_x[i] = 0.0;
        right_x[i] = H;
        left_y[i] = dy1 * i;
        right_y[i] = dL + dy2 * i;
    }
}


void grid_generation::generate_grid_by_transfinite_interpolation(
    const std::vector<double> &bottom_x,
    const std::vector<double> &bottom_y,
    const std::vector<double> &top_x,
    const std::vector<double> &top_y,
    const std::vector<double> &left_x,
    const std::vector<double> &left_y,
    const std::vector<double> &right_x,
    const std::vector<double> &right_y,
    std::vector<double> &x,
    std::vector<double> &y
) {
    unsigned int N = left_x.size() - 1;
    unsigned int M = bottom_x.size() - 1;
    double dxi = 1.0 / M, deta = 1.0 / N;
    x.resize((M + 1) * (N + 1)), y.resize((M + 1) * (N + 1));
    for (unsigned int i = 0; i <= N; ++i) {
        double eta = deta * i;
        for (unsigned int j = 0; j <= M; ++j) {
            double xi = dxi * j;
            unsigned int ij = i * (M + 1) + j;
            x[ij] = (1 - xi) * left_x[i] + xi * right_x[i]
                     + (1 - eta) * bottom_x[j] + eta * bottom_x[j]
                     - (1- xi) * (1 - eta) * bottom_x[0] - (1 - xi) * eta * bottom_x[0]
                     - (1 - eta) * xi * bottom_x[M] - xi * eta * bottom_x[M];
            y[ij] = (1 - xi) * left_y[i] + xi * right_y[i] +
                      (1 - eta) * bottom_y[j] + eta * top_y[j] -
                      (1- xi) * (1 - eta) * bottom_y[0] - (1 - xi) * eta * top_y[0] -
                      (1 - eta) * xi * bottom_y[M] - xi * eta * top_y[M];
        }        
    }
}



void grid_generation::generate_stretching_grid(
    double alpha,
    double eta1,
    unsigned int N,
    const std::vector<double> bottom_x,
    const std::vector<double> bottom_y,
    const std::vector<double> top_x,
    const std::vector<double> top_y,
    std::vector<double> &x,
    std::vector<double> &y
) {
    unsigned int M = bottom_x.size() - 1;
    x.resize((N + 1) * (M + 1)), y.resize((N + 1) * (M + 1));
    for (unsigned int i = 0; i <= N; ++i) {
        double eta = i * 1.0 / N;
        for (unsigned int j = 0; j <= M; ++j) {
            unsigned int ij = i * (M + 1) + j;
            x[ij] = bottom_x[j];
            if (eta <= eta1)
                y[ij] = (top_y[j] - bottom_y[j]) * eta1 * (std::exp(alpha * eta / eta1) - 1.0) / (std::exp(alpha) - 1) + bottom_y[j];
            else
                y[ij] = (top_y[j] - bottom_y[j]) * (1.0 - (1.0 - eta1) * (std::exp(alpha * (1 - eta) / (1 - eta1)) - 1.0) / (std::exp(alpha) - 1)) + bottom_y[j];
        }
    }
}



void grid_generation::generate_grid_by_hermite_interpolation(
    unsigned int N, unsigned int M,
    double Q1, double Q2,
    double T1, double T2,
    std::function<double(double)> xAB, std::function<double(double)> yAB,
    std::function<double(double)> xCD, std::function<double(double)> yCD,
    std::function<double(double)> dxAB, std::function<double(double)> dyAB,
    std::function<double(double)> dxCD, std::function<double(double)> dyCD,
    std::vector<double> &x, std::vector<double> &y
) {
    x.resize((N + 1) * (M + 1));
    y.resize((N + 1) * (M + 1));

    auto h = [](double xi, double Q) -> double { return (std::tanh(Q * (xi - 0.5)) - std::tanh(-0.5 * Q)) / (2.0 * std::tanh(0.5 * Q)); };
    
    auto Psi1 = [](double eta) -> double { return (2 * eta + 1) * (eta - 1) * (eta - 1); };
    auto Psi2 = [](double eta) -> double { return eta * eta * (-2 * eta + 3); };
    auto Psi3 = [](double eta) -> double { return eta * (eta - 1) * (eta - 1); };
    auto Psi4 = [](double eta) -> double { return eta * eta * (eta - 1); };

    for (unsigned int j = 0; j <= M; ++j) {
        double xi = static_cast<double>(j) / M;
        double hAB = h(xi, Q1);
        double hCD = h(xi, Q2);
        for (unsigned int i = 0; i <= N; ++i) {
            double eta = static_cast<double>(i) / N;
            unsigned int ij = i * (M + 1) + j;
            x[ij] = Psi1(eta) * xAB(hAB) + Psi2(eta) * xCD(hCD) + T1 * Psi3(eta) * dyAB(hAB) + T2 * Psi4(eta) * dyCD(hCD);
            y[ij] = Psi1(eta) * yAB(hAB) + Psi2(eta) * yCD(hCD) - T1 * Psi3(eta) * dxAB(hAB) - T2 * Psi4(eta) * dxCD(hCD);
        }
    }
}



void grid_generation::generate_full_grid(
    double H1,
    double H2,
    double H3,
    unsigned int N,
    unsigned int M1,
    unsigned int M2,
    unsigned int M3,
    double Q1,
    double P2,
    double Q2,
    const std::vector<double> &x2,
    const std::vector<double> &y2,
    std::vector<double> &x,
    std::vector<double> &y
) {
    
    auto h1 = [](double xi, double Q) -> double { return (std::tanh(Q * (xi - 0.5)) - std::tanh(-0.5 * Q)) / (2.0 * std::tanh(0.5 * Q)); };
    auto h2 = [](double xi, double P, double Q) -> double { return P * xi + (1 - P) * (1 - std::tanh(Q * (1 - xi)) / std::tanh(Q)); };

    unsigned int M = M1 + M2 + M3;
    unsigned int num_nodes = (M + 1) * (N + 1);

    double dx1 = H1 / M1, dx3 = H3 / M3;
    x.resize(num_nodes), y.resize(num_nodes);
    for (unsigned int i = 0; i <= N; ++i) {
        unsigned int ij2l = i * (M2 + 1);
        unsigned int ij2r = i * (M2 + 1) + M2;
        double y_left = y2[ij2l];
        double y_right = y2[ij2r];
        for (unsigned int j = 0; j < M1; ++j) {
            unsigned int ij = i * (M + 1) + j;
            double xi = static_cast<double>(j) / M1;
            x[ij] = H1 * h1(xi, Q1);
            y[ij] = y_left;
        }
        for (unsigned int j = M1; j <= M1 + M2; ++j) {
            unsigned int ij = i * (M + 1) + j;
            unsigned int ij2 = i * (M2 + 1) + (j - M1);
            x[ij] = x2[ij2];
            y[ij] = y2[ij2];
        }
        for (unsigned int j = M1 + M2 + 1; j <= M; ++j) {
            unsigned int ij = i * (M + 1) + j;
            double xi = static_cast<double>(j - M1 - M2) / M3;
            x[ij] = H1 + H2 + H3 * h2(xi, P2, Q2);
            y[ij] = y_right;
        }
    }
}



void grid_generation::generate_grid_connection(
    unsigned int N,
    unsigned int M,
    std::vector<std::array<unsigned int, 3>> &connection
) {
    unsigned int elem_id = 0;
    connection.resize(N * M * 2);
    for (unsigned int i = 0; i < N; ++i) {
        for (unsigned int j = 0; j < M; ++j) {
            unsigned int n1 = i * (M + 1) + j;
            unsigned int n2 = n1 + 1;
            unsigned int n3 = n2 + M + 1;
            unsigned int n4 = n1 + M + 1;

            connection[elem_id][0] = n1;
            connection[elem_id][1] = n2;
            connection[elem_id][2] = n3;
            elem_id++;

            connection[elem_id][0] = n1;
            connection[elem_id][1] = n3;
            connection[elem_id][2] = n4;
            elem_id++;
        }
    }

    for (unsigned int i = 0; i < N * M * 2; ++i) {
            std::cout << "[" << connection[i][0] << ", " << connection[i][1] << ", " << connection[i][2] << "]" << std::endl;
    }
}



#endif