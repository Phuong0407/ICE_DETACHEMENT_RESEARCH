#ifndef AXISYMMETRIC_FINITE_ELEMENT_METHOD_H
#define AXISYMMETRIC_FINITE_ELEMENT_METHOD_H

// #include "geometry.h"
#include "csr_sparse_solver.hpp"
#include <cmath>
#include <vector>
#include <array>
#include <tuple>

using Weight = double;
using GaussPoint = double;
using GaussQuadraturePoint = std::pair<GaussPoint, Weight>;
struct GaussQuadrature8Point {
    static constexpr std::array<GaussQuadraturePoint, 8> Point = {{
        { -0.96028985649753623168, 0.10122853629037625915 },
        { -0.79666647741362673959, 0.22238103445337447054 },
        { -0.52553240991632898582, 0.31370664587788728734 },
        { -0.18343464249564980494, 0.36268378337836209114 },
        {  0.18343464249564980494, 0.36268378337836209114 },
        {  0.52553240991632898582, 0.31370664587788728734 },
        {  0.79666647741362673959, 0.22238103445337447054 },
        {  0.96028985649753623168, 0.10122853629037625915 }
    }};
};
const GaussQuadrature8Point Gauss8Point;

inline double triangleArea(double x1, double y1, double x2, double y2, double x3, double y3) {
    return 0.5 * std::abs(x1 * (y2 - y3) + x2 * (y3 - y1) + x3 * (y1 - y2));
}

class AxisymmetricFiniteElement {
private:
    std::vector<std::vector<double>> GloMat;
    std::vector<double> ForceMat;
    std::vector<double> Solution;
    std::tuple<double, double, double> computeShapeFunctionContributions(double x1, double y1, double x2, double y2, double x3, double y3, const GaussQuadrature8Point& Gauss8Point, double Flux);

public:
    AxisymmetricFiniteElement() = default;
    void computeGlobalStiffnessMatrix(const std::vector<double> &x, const std::vector<double> &y, const std::vector<std::vector<std::size_t>> &ElemConnData);
    std::vector<std::vector<double>> computeLocalStiffnessMatrix(double x1, double y1, double x2, double y2, double x3, double y3);
    void computeNeumannBoundaryConditions(const std::vector<double> &x, const std::vector<double> &y, const std::vector<std::vector<std::size_t>>& ElemConnData, const std::vector<BoundaryEdge> &NeumannBoundaryConditions);
    void computeSolution();
    std::vector<double> getSolution() const {return Solution;}
};

void AxisymmetricFiniteElement::computeGlobalStiffnessMatrix(const std::vector<double> &x, const std::vector<double> &y, const std::vector<std::vector<std::size_t>> &ElemConnData) {
    std::size_t N = x.size();
    GloMat.clear(), GloMat.resize(N, std::vector<double>(N, 0.0));
    std::size_t NumOfElem = ElemConnData.size();
    for (std::size_t IndElem = 0; IndElem < NumOfElem; ++IndElem) {
        const auto &ElemConn = ElemConnData[IndElem];
        double x1 = x[ElemConn[0]], y1 = y[ElemConn[0]];
        double x2 = x[ElemConn[1]], y2 = y[ElemConn[1]];
        double x3 = x[ElemConn[2]], y3 = y[ElemConn[2]];
        std::vector<std::vector<double>> LocMat = computeLocalStiffnessMatrix(x1, y1, x2, y2, x3, y3);

        for (std::size_t LocRowInd = 0; LocRowInd < 3; ++LocRowInd) {
            std::size_t GloRowInd = ElemConnData[IndElem][LocRowInd];
            for (std::size_t LocColInd = 0; LocColInd < 3; ++LocColInd) {
                std::size_t GloColInd = ElemConnData[IndElem][LocColInd];
                GloMat[GloRowInd][GloColInd] += LocMat[LocRowInd][LocColInd];
            }
        }
    }
    // for (const auto& elem : GloMat) {
    //     for (const auto & a : elem) {
    //         std::cout << a << " ";
    //     } std::cout << std::endl;
    // }
}

std::vector<std::vector<double>> AxisymmetricFiniteElement::computeLocalStiffnessMatrix(double x1, double y1, double x2, double y2, double x3, double y3) {
    std::vector<double> c(3), d(3);
    std::vector<std::vector<double>> LocMat(3, std::vector<double>(3, 0.0));

    c[0] = y2 - y3, c[1] = y3 - y1, c[2] = y1 - y2;
    d[0] = x3 - x2, d[1] = x1 - x3, d[2] = x2 - x1;

    double r_centroid = (y1 + y2 + y3) / 3.0;
    double A = triangleArea(x1, y1, x2, y2, x3, y3);
    double Coeff = M_PI * r_centroid / (2 * A);

    for (std::size_t i = 0; i < 3; ++i) {
        for (std::size_t j = 0; j < 3; ++j) {
            LocMat[i][j] = Coeff * (c[i] * c[j] + d[i] * d[j]);
        }
    }
    return LocMat;
}

void AxisymmetricFiniteElement::computeNeumannBoundaryConditions(const std::vector<double> &x, const std::vector<double> &y, const std::vector<std::vector<std::size_t>> &ElemConnData, const std::vector<BoundaryEdge> &NeumannBoundaryConditions) {
    ForceMat.clear();
    ForceMat.resize(x.size(), 0.0);

    std::size_t n = NeumannBoundaryConditions.size();

    for (std::size_t i = 0; i < n; ++i) {
        double Flux = NeumannBoundaryConditions[i].get_flux();
        if (std::abs(Flux) < 1e-12) continue;

        std::size_t node1 = NeumannBoundaryConditions[i].get_first_node_ind();
        std::size_t node2 = NeumannBoundaryConditions[i].get_second_node_ind();
        std::size_t node3 = 0;

        double x1 = 0.0, y1 = 0.0, x2 = 0.0, y2 = 0.0, x3 = 0.0, y3 = 0.0;

        for (const auto& element : ElemConnData) {
            if (element[0] == node1 && element[1] == node2) {
                node3 = element[2];
                x1 = x[node1], y1 = y[node1];
                x2 = x[node2], y2 = y[node2];
                x3 = x[node3], y3 = y[node3];
            } 
            else if (element[1] == node1 && element[2] == node2) {
                node3 = element[0];
                x1 = x[node1], y1 = y[node1];
                x2 = x[node2], y2 = y[node2];
                x3 = x[node3], y3 = y[node3];
            } 
            else if (element[2] == node1 && element[0] == node2) {
                node3 = element[1];
                x1 = x[node1], y1 = y[node1];
                x2 = x[node2], y2 = y[node2];
                x3 = x[node3], y3 = y[node3];
            }
        }
        auto [N1, N2, N3] = computeShapeFunctionContributions(x1, y1, x2, y2, x3, y3, Gauss8Point, Flux);
        ForceMat[node1] += N1;
        ForceMat[node2] += N2;
        ForceMat[node3] += N3;
    }
    // for (std::size_t i = 0; i < ForceMat.size(); ++i) {
    //     std::cout << ForceMat[i] << std::endl;
    // }
}

std::tuple<double, double, double> AxisymmetricFiniteElement::computeShapeFunctionContributions(double x1, double y1, double x2, double y2, double x3, double y3, const GaussQuadrature8Point& Gauss8Point, double Flux) {
    double b1 = x2 * y3 - x3 * y2, b2 = x3 * y1 - x1 * y3, b3 = x1 * y2 - x2 * y1;
    double c1 = y2 - y3, c2 = y3 - y1, c3 = y1 - y2;
    double d1 = x3 - x2, d2 = x1 - x3, d3 = x2 - x1;

    double Area = triangleArea(x1, y1, x2, y2, x3, y3);
    double edge_length = std::sqrt((x1 - x2) * (x1 - x2) + (y1 - y2) * (y1 - y2));
    double N1 = 0.0, N2 = 0.0, N3 = 0.0;

    for (const auto& GaussPoint : Gauss8Point.Point) {
        double xi = GaussPoint.first;
        double weight = GaussPoint.second;
        double x_ = (1 - xi) / 2.0 * x1 + (1 + xi) / 2.0 * x2;
        double y_ = (1 - xi) / 2.0 * y1 + (1 + xi) / 2.0 * y2;
        N1 += (b1 + c1 * x_ + d1 * y_) * y_ * weight;
        N2 += (b2 + c2 * x_ + d2 * y_) * y_ * weight;
        N3 += (b3 + c3 * x_ + d3 * y_) * y_ * weight;
    }
    double factor = M_PI * Flux / Area;
    N1 *= factor * edge_length / 2.0;
    N2 *= factor * edge_length / 2.0;
    N3 *= factor * edge_length / 2.0;
    return {N1, N2, N3};
}

void AxisymmetricFiniteElement::computeSolution() {
    Solution.clear();
    CSRSparseSolver Solver;
    Solver.initSparseSolver(GloMat);
    Solver.solveSparseSystem(ForceMat, Solution);
}

#endif