#ifndef AXISYMMETRIC_FINITE_ELEMENT_METHOD_H
#define AXISYMMETRIC_FINITE_ELEMENT_METHOD_H

#include "geometry.h"
#include <cmath>
#include <vector>
#include <array>

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

using ColMat = std::vector<double>;
using RowMat = std::vector<double>;
using Mat = std::vector<RowMat>;

inline double triangleArea(double x1, double y1, double x2, double y2, double x3, double y3) {
    return 0.5 * std::abs(x1 * (y2 - y3) + x2 * (y3 - y1) + x3 * (y1 - y2));
}

class AxisymmetricFiniteElement {
private:
    Mat GloMat;
    ColMat ForceMat;
    ColMat Solution;

public:
    AxisymmetricFiniteElement() = default;
    void computeGlobalStiffnessMatrix(const CoordArr &x, const CoordArr &y, const ElemConnArr &ElemConnData);
    Mat computeLocalStiffnessMatrix(double x1, double y1, double x2, double y2, double x3, double y3);
    void computeNeumannBoundaryConditions(const CoordArr &x, const CoordArr &y, const ElemConnArr& ElemConnData, const std::vector<BoundaryEdge> &NeumannBoundaryConditions);
};

void AxisymmetricFiniteElement::computeGlobalStiffnessMatrix(const CoordArr &x, const CoordArr &y, const ElemConnArr &ElemConnData) {
    std::size_t N = x.size();
    GloMat.clear(), GloMat.resize(N, RowMat(N, 0.0));

    for (std::size_t IndElem = 0; IndElem < ElemConnData.size(); ++IndElem) {
        std::vector<Ind> ind_node;
        for (const auto ElemConn : ElemConnData[IndElem])
            ind_node.push_back(ElemConn);

        double x1 = x[ind_node[0]], y1 = y[ind_node[0]];
        double x2 = x[ind_node[0]], y2 = y[ind_node[0]];
        double x3 = x[ind_node[0]], y3 = y[ind_node[0]];
        Mat LocMat = computeLocalStiffnessMatrix(x1, y1, x2, y2, x3, y3);

        for (Ind LocRowInd = 0; LocRowInd < 3; ++LocRowInd) {
            Ind GloRowInd = ElemConnData[IndElem][LocRowInd];
            for (Ind LocColInd = 0; LocColInd < 3; ++LocColInd) {
                Ind GloColInd = ElemConnData[IndElem][LocColInd];
                GloMat[GloRowInd][GloColInd] += LocMat[LocRowInd][LocColInd];
            }
        }
    }
}

Mat AxisymmetricFiniteElement::computeLocalStiffnessMatrix(double x1, double y1, double x2, double y2, double x3, double y3) {
    RowMat c(3), d(3);
    Mat LocMat(3, RowMat(3, 0.0));

    c[0] = x2 - x3, c[1] = x3 - x1, c[2] = x1 - x2;
    d[0] = y3 - y1, d[1] = y1 - y3, d[2] = y2 - y1;

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

void AxisymmetricFiniteElement::computeNeumannBoundaryConditions(const CoordArr &x, const CoordArr &y, const ElemConnArr& ElemConnData, const std::vector<BoundaryEdge> &NeumannBoundaryConditions) {
    ForceMat.clear(), ForceMat.resize(x.size(), 0.0);
    Size n = NeumannBoundaryConditions.size();
    for (std::size_t i = 0; i < n; ++i) {
        double Flux = NeumannBoundaryConditions[i].get_flux();
        if (std::abs(Flux) < 1e-12) continue;
        std::size_t node1 = NeumannBoundaryConditions[i].get_first_node_ind();
        std::size_t node2 = NeumannBoundaryConditions[i].get_second_node_ind();
        for (const auto& element : ElemConnData) {
            double x1, y1, x2, y2, x3, y3;
            if (element[0] == node1 && element[1] == node2) {
                x1 = x[node1], y1 = x[node1];
                x2 = x[node2], y2 = y[node2];
                x3 = x[element[2]], y3 = y[element[2]];
            } else if (element[1] == node1 && element[2] == node2) {
                x2 = x[node1], y2 = x[node1];
                x3 = x[node2], y3 = y[node2];
                x1 = x[element[0]], y2 = y[element[0]];

            } else if (element[2] == node1 && element[0] == node2) {
                x3 = x[node1], y3 = x[node1];
                x1 = x[node2], y1 = y[node2];
                x2 = x[element[1]], y2 = y[element[1]];
            }
            double b1 = y2 * x3 - y3 * x3, b2 = y3 * x1 - y1 * x3, b3 = y1 * x2 - y2 * x1;
            double c1 = x2 - x3, c2 = x3 - x1, c3 = x1 - x2;
            double d1 = y3 - y1, d2 = y1 - y3, d3 = y2 - y1;
            double Area = triangleArea(x1, y1, x2, y2, x3, y3);
            double edge_length = std::sqrt((x1 - x2) * (x1 - x2) + (y1 - y2) * (y1 - y2));
            double N1 = 0.0, N2 = 0.0, N3 = 0.0;
            for (const auto& GaussPoint : Gauss8Point.Point) {
                double xi = GaussPoint.first;
                double weight = GaussPoint.second;
                double x_ = (1 - xi) / 2.0 * x1 +  xi / 2.0 * x2;
                double y_ = (1 - xi) / 2.0 * x2 +  xi / 2.0 * y2;
                N1 += (b1 + d1 * x_ + c1 * y_) * y_ * weight;
                N2 += (b2 + d2 * x_ + c2 * y_) * y_ * weight;
                N3 += (b3 + d3 * x_ + c3 * y_) * y_ * weight;
            }
            N1 *= - M_PI * Flux / Area * edge_length / 2.0;
            N2 *= - M_PI * Flux / Area * edge_length / 2.0;
            N3 *= - M_PI * Flux / Area * edge_length / 2.0;
            if (element[0] == node1 && element[1] == node2) {
                ForceMat[node1] += N1;
                ForceMat[node2] += N2;
                ForceMat[element[2]] += N3;
            } else if (element[1] == node1 && element[2] == node2) {
                ForceMat[node1] += N2;
                ForceMat[node2] += N3;
                ForceMat[element[0]] += N1;
            } else if (element[2] == node1 && element[0] == node2) {
                ForceMat[node1] += N3;
                ForceMat[node2] += N1;
                ForceMat[element[1]] += N2;
            }
        }
    }
}

#endif