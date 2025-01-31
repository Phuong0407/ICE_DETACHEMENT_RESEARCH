#ifndef AXISYMMETRIC_FINITE_ELEMENT_METHOD_H
#define AXISYMMETRIC_FINITE_ELEMENT_METHOD_H

#include <cmath>
#include <vector>

#ifndef GEOMETRY_H

using Size = std::size_t;
using Ind = std::size_t;
using CoordArr = std::vector<double>;
using Coord2DArr = std::vector<CoordArr>;
using ElemConnArr = std::vector<std::vector<unsigned int>>;

#endif

using ColMat = std::vector<double>;
using RowMat = std::vector<double>;
using Mat = std::vector<RowMat>;

inline double triangleArea(double x1, double y1, double x2, double y2, double x3, double y3) {
    return 0.5 * std::abs(x1 * (y2 - y3) + x2 * (y3 - y1) + x3 * (y1 - y2));
}

class AxisymmetricFiniteElement {
private:
    Mat GloMat;
    ColMat NeumannBoundaryCondition;
    ColMat Solution;

public:
    AxisymmetricFiniteElement() = default;
    void computeGlobalStiffnessMatrix(CoordArr x, CoordArr y, ElemConnArr ElemConnData);
    Mat computeLocalStiffnessMatrix(double x1, double y1, double x2, double y2, double x3, double y3);
};

void AxisymmetricFiniteElement::computeGlobalStiffnessMatrix(CoordArr x, CoordArr y, ElemConnArr ElemConnData) {
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


#endif