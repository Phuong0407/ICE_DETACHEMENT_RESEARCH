/**
 * The flow domain is define as follow
 * 6------------------------------------5
 * |                                    |
 * |                                    |
 * |                                    |
 * |                                    |
 * |                  3-----------------4
 * |                 /
 * |                /
 * 1---------------2
 * The domain is symmetric about the axis 1-2 which has y = 0
 * The edge 1 - 2 - 3 - 4 has zero flux
 * The edge 5 - 6 has wall condition
 * The edge 1 - 6 has inlet conditon, velocity U_in
 * The edge 4 - 5 has outlet conditon, velocity U_out * 
 */


#ifndef GEOMETRY_H
#define GEOMETRY_H

#include <cmath>
#include <vector>
#include <cstddef>
#include <iostream>
#include <stdexcept>

using Size = std::size_t;
using Ind = std::size_t;
using CoordArr = std::vector<double>;
using Coord2DArr = std::vector<CoordArr>;
using ElemConnArr = std::vector<std::vector<unsigned int>>;

class Geometry {
private:
    double L1, L2, H1, H2, H3;
    double dx1, dx2, dx3, dy1, dy2;
    Size N, M, M1, M2, M3;
    Coord2DArr x, y;
    CoordArr x_leftbound, y_leftbound;
    CoordArr x_rightbound, y_rightbound;
    CoordArr x_bottombound, y_bottombound;
    CoordArr x_topbound, y_topbound;
    ElemConnArr ElemConnData;

public:
    Geometry() = default;
    void initPhysicalDimensions(double H1, double H2, double H3, double L1, double L2);
    void initComputationalDimensions(Size M1, Size M2, Size M3, Size N);
    void initGridStep();
    void initCoordsArrSize();
    void generateAlgebraicGrid();
    void generateEllipicGrid();
    void generateBoundaryNode();
    void transfiniteInterpolation();
    void generateStretchingGrid(double alpha, double eta1);
    void generateGridConnection();

public:
    const Coord2DArr& get_x() const {return x;}
    const Coord2DArr& get_y() const {return y;}
    const ElemConnArr& getElementConnectionData() const {return ElemConnData;}
    Size getTotalHorizontalElements() const {return M;}
    Size getTotalVerticalElements() const {return N;}
};

void Geometry::initPhysicalDimensions(double H1, double H2, double H3, double L1, double L2) {
    this->H1 = H1;
    this->H2 = H2;
    this->H3 = H3;
    this->L1 = L1;
    this->L2 = L2;
}

void Geometry::initComputationalDimensions(Size M1, Size M2, Size M3, Size N) {
    this->M1 = M1;
    this->M2 = M2;
    this->M3 = M3;
    this->N = N;
    this->M = M1 + M2 + M3;
}

void Geometry::initGridStep() {
    dx1 = H1 / M1;
    dx2 = H2 / M2;
    dx3 = H3 / M3;
    dy1 = L1 / N;
    dy2 = L2 / N;
}

void Geometry::initCoordsArrSize() {
    x.resize(N + 1, std::vector<double>(M + 1));
    y.resize(N + 1, std::vector<double>(M + 1));
    x_leftbound.resize(N + 1);
    y_leftbound.resize(N + 1);
    x_rightbound.resize(N + 1);
    y_rightbound.resize(N + 1);
    x_bottombound.resize(M + 1);
    y_bottombound.resize(M + 1);
    x_topbound.resize(M + 1);
    y_topbound.resize(M + 1);
}

void Geometry::generateAlgebraicGrid() {
    for (Ind i = 0; i <= N; ++i) {
        for (Ind j = 0; j <= M1; ++j)
            x[i][j] = dx1 * j, y[i][j] = dy1 * i;

        double a = ((dy2 - dy1) * i + L1 - L2) / H2;
        double b = dy1 * i - a * H1;
        for (Ind j = M1 + 1; j <= M1 + M2; ++j) {
            double x_1 = H1 + dx2 * (j - M1);
            double y_1 = a * x_1 + b;
            x[i][j] = x_1;
            y[i][j] = y_1;
        }

        for (Ind j = M1 + M2 + 1; j <= M; ++j) {
            x[i][j] = H1 + H2 + dx3 * (j - M1 - M2);
            y[i][j] = L1 - L2 + dy2 * i;
        }
    }
}

void Geometry::generateEllipicGrid() {
    if (x.empty() || y.empty())
        throw std::runtime_error("ERROR: ALGEBRAIC-GENERATING COORDINATES ARRAY IS EMPTY. GENERATE IT FIRST BY USING generateAlgebraicGrid()!");

    Coord2DArr new_x(N + 1, std::vector<double>(M + 1, 0.0));
    Coord2DArr new_y(N + 1, std::vector<double>(M + 1, 0.0));
    Coord2DArr alpha(N + 1, std::vector<double>(M + 1, 0.0));
    Coord2DArr beta(N + 1, std::vector<double>(M + 1, 0.0));
    Coord2DArr gamma(N + 1, std::vector<double>(M + 1, 0.0));
    new_x = x, new_y = y;

    double min_err = 1e-6;
    std::size_t max_iter = 10000;
    std::vector<double> error_1(max_iter, 0.0);
    std::vector<double> error_2(max_iter, 0.0);

    for (std::size_t iter = 0; iter < max_iter; ++iter) {
        for (int i = 1; i < N; ++i) {
            for (int j = 1; j < M; ++j) {
                double Dxi_x = (x[i][j + 1] - x[i][j - 1]);
                double Dxi_y = (y[i][j + 1] - y[i][j - 1]);
                double Deta_x = (x[i + 1][j] - x[i - 1][j]);
                double Deta_y = (y[i + 1][j] - y[i - 1][j]);
                alpha[i][j] = 0.25 * (std::pow(Deta_x, 2.0) + std::pow(Deta_y, 2.0));
                gamma[i][j] = 0.25 * (std::pow(Dxi_x, 2.0) + std::pow(Dxi_y, 2.0));
                beta[i][j] = 0.0625 * (Deta_x * Dxi_x + Deta_y * Dxi_y);
                double denominator = alpha[i][j] + gamma[i][j] + 1e-9;
                new_x[i][j] = 0.5 / denominator * (-2.0 * beta[i][j] * (x[i + 1][j + 1] - x[i + 1][j - 1] - x[i - 1][j + 1] + x[i - 1][j - 1]) +
                              alpha[i][j] * (x[i][j + 1] + x[i][j - 1]) + gamma[i][j] * (x[i + 1][j] + x[i - 1][j]));
                new_y[i][j] = 0.5 / denominator * (-2.0 * beta[i][j] * (y[i + 1][j + 1] - y[i + 1][j - 1] - y[i - 1][j + 1] + y[i - 1][j - 1]) +
                              alpha[i][j] * (y[i][j + 1] + y[i][j - 1]) + gamma[i][j] * (y[i + 1][j] + y[i - 1][j]));
            }
        }
        double max_error_x = 0.0, max_error_y = 0.0;
        for (std::size_t i = 1; i < N; ++i) {
            for (std::size_t j = 1; j < M; ++j) {
                max_error_x = std::max(max_error_x, std::abs(new_x[i][j] - x[i][j]));
                max_error_y = std::max(max_error_y, std::abs(new_y[i][j] - y[i][j]));
            }
        }
        error_1[iter] = max_error_x;
        error_2[iter] = max_error_y;
        x = new_x;
        y = new_y;
        if (error_1[iter] < min_err && error_2[iter] < min_err)
            break;
    }
}

void Geometry::transfiniteInterpolation() {
    double deta = 1.0 / N;
    double dxi = 1.0 / M;
    for (std::size_t i = 0; i < N + 1; ++i) {
        double eta = deta * i;
        for (std::size_t j = 0; j < M + 1; ++j) {
            double xi = dxi * j;
            x[i][j] = (1 - xi) * x_leftbound[i] + xi * x_rightbound[i]
                     + (1 - eta) * x_bottombound[j] + eta * x_topbound[j]
                     - (1- xi) * (1 - eta) * x_bottombound[0] - (1 - xi) * eta * x_topbound[0]
                     - (1 - eta) * xi * x_bottombound[M] - xi * eta * x_topbound[M];
            y[i][j] = (1 - xi) * y_leftbound[i] + xi * y_rightbound[i] +
                      (1 - eta) * y_bottombound[j] + eta * y_topbound[j] -
                      (1- xi) * (1 - eta) * y_bottombound[0] - (1 - xi) * eta * y_topbound[0] -
                      (1 - eta) * xi * y_bottombound[M] - xi * eta * y_topbound[M];
        }        
    }
}

void Geometry::generateBoundaryNode() {
    for (std::size_t i = 0; i <= M1; ++i) {
        double x = H1/M1 * i;
        x_bottombound[i] = x;
        y_bottombound[i] = 0.0;
        x_topbound[i] = x;
        y_topbound[i] = L1;
    }
    double a = (L1 - L2) / H2;
    for (std::size_t i = M1 + 1; i <= M1 + M2; ++i) {
        double x = H1 + dx2 * (i - M1);
        double y = a * dx2 * (i - M1);
        x_bottombound[i] = x;
        y_bottombound[i] = y;
        x_topbound[i] = x;
        y_topbound[i] = L1;
    }
    for (std::size_t i = M1 + M2 + 1; i <= M; ++i) {
        double x = H1 + H2 + dx3 * (i - M1 - M2);
        x_bottombound[i] = x;
        y_bottombound[i] = L1 - L2;
        x_topbound[i] = x;
        y_topbound[i] = L1;
    }
    for (std::size_t i = 0; i <= N; ++i) {
        x_leftbound[i] = 0.0;
        x_rightbound[i] = H1 + H2 + H3;
        y_leftbound[i] = L1/N * i;
        y_rightbound[i] = L1 - L2 + L2/N * i;
    }
}

void Geometry::generateStretchingGrid(double alpha, double eta1) {
    double deta = 1.0 / N;
    double dxi = 1.0 / M;
    for (std::size_t i = 0; i < N + 1; ++i) {
        double eta = i * deta;
        for (std::size_t j = 0; j <= M; ++j) {
            double xi = j * dxi;
            x[i][j] = x_bottombound[j];
            if (eta <= eta1) {
                y[i][j] = (y_topbound[j] - y_bottombound[j]) * eta1 * (std::exp(alpha * eta / eta1) - 1.0) / (std::exp(alpha) - 1) + y_bottombound[j];
            }
            else {
                y[i][j] = (y_topbound[j] - y_bottombound[j]) * (1.0 - (1.0 - eta1) * (std::exp(alpha * (1 - eta) / (1 - eta1)) - 1.0) / (std::exp(alpha) - 1)) + y_bottombound[j];
            }
        }
    }
}


void Geometry::generateGridConnection() {
    for (std::size_t i = 0; i < N; ++i) {
        for (std::size_t j = 0; j < M; ++j) {
            unsigned int ind_elem = i * (M + 1) + j;
            unsigned int first_node_ind = ind_elem;
            unsigned int second_node_ind = first_node_ind + 1;
            unsigned int third_node_ind = second_node_ind + M + 1;
            unsigned int fourth_node_ind = first_node_ind + M + 1;
            ElemConnData.push_back({first_node_ind, second_node_ind, third_node_ind});
            ElemConnData.push_back({first_node_ind, third_node_ind, fourth_node_ind});
        }
    }
}

#endif