#ifndef POST_PROCESSING_H
#define POST_PROCESSING_H

#include "dense_system_solver.hpp"

#include <array>
#include <vector>
#include <algorithm>
#include <unordered_set>
#include <iostream>
#include <vector>

std::vector<std::vector<double>> computeTransposeProduct(const std::vector<std::vector<double>> &A) {
    size_t m = A.size();
    size_t n = A[0].size();
    std::vector<std::vector<double>> ATA(n, std::vector<double>(n, 0.0));
    for (size_t i = 0; i < n; ++i) {
        for (size_t j = 0; j < n; ++j) {
            for (size_t k = 0; k < m; ++k) {
                ATA[i][j] += A[k][i] * A[k][j];
            }
        }
    }
    return ATA;
}

std::vector<double> computeTransposeVectorProduct(const std::vector<std::vector<double>> &A, const std::vector<double> &B) {
    size_t m = A.size();
    size_t n = A[0].size();

    std::vector<double> ATB(n, 0.0);

    for (size_t i = 0; i < n; ++i) {
        for (size_t k = 0; k < m; ++k) {
            ATB[i] += A[k][i] * B[k];
        }
    }
    return ATB;
}

class RecoveryPatchConnManager {
protected:
    const std::size_t OrderOfPatch = 6;
    std::size_t number_nodes;
    std::vector<std::unordered_set<std::size_t>> FirstLevelNodePatch, FullLevelNodePatch;

private:
    std::vector<std::unordered_set<std::size_t>> genFirstLevelPatchNode(const std::vector<std::vector<std::size_t>> &ElemConnData, std::size_t number_nodes);
    std::unordered_set<std::size_t> findAllInternalNode(const std::vector<std::unordered_set<std::size_t>> &FirstLevelNodePatch, const std::unordered_set<std::size_t> &FirstOrderPatch);

public:
    RecoveryPatchConnManager() = default;
    void initNumberOfNode(std::size_t number_nodes) {this->number_nodes = number_nodes;}
    void generateNodeRecoveryPatchConnManager(const std::vector<std::vector<std::size_t>> &ElemConnData);
};

void RecoveryPatchConnManager::generateNodeRecoveryPatchConnManager(const std::vector<std::vector<std::size_t>> &ElemConnData) {
    FirstLevelNodePatch = genFirstLevelPatchNode(ElemConnData, number_nodes);
    FullLevelNodePatch = FirstLevelNodePatch;

    for (std::size_t i = 0; i < number_nodes; ++i) {
        const auto &NodePatch = FirstLevelNodePatch[i];
        if (NodePatch.size() >= OrderOfPatch) continue;
        std::unordered_set<std::size_t> InternalNodeForPatch = findAllInternalNode(FirstLevelNodePatch, NodePatch);
        for (const auto &InternalNode : InternalNodeForPatch) {
            FullLevelNodePatch[i].insert(FirstLevelNodePatch[InternalNode].begin(), FirstLevelNodePatch[InternalNode].end());
        } FullLevelNodePatch[i].erase(i);
    }

    for (std::size_t i = 0; i < number_nodes; ++i) {
        const auto &NodePatch = FullLevelNodePatch[i];
        if (NodePatch.size() < OrderOfPatch) {
            for (const auto &neighbor : NodePatch) {
                FullLevelNodePatch[i].insert(FullLevelNodePatch[neighbor].begin(), FullLevelNodePatch[neighbor].end());
                FullLevelNodePatch[i].erase(i);
            }
        }
    }
}

std::vector<std::unordered_set<std::size_t>> RecoveryPatchConnManager::genFirstLevelPatchNode(const std::vector<std::vector<std::size_t>> &ElemConnData, std::size_t number_nodes) {
    std::vector<std::unordered_set<std::size_t>> FirstLevelNodePatch(number_nodes);
    for (std::size_t i = 0; i < number_nodes; ++i) {
        for (const auto& element : ElemConnData) {
            if (element[0] == i || element[1] == i || element[2] == i) {
                FirstLevelNodePatch[i].insert(element[0]);
                FirstLevelNodePatch[i].insert(element[1]);
                FirstLevelNodePatch[i].insert(element[2]);
            }
        }
    }
    return FirstLevelNodePatch;
}

std::unordered_set<std::size_t> RecoveryPatchConnManager::findAllInternalNode(const std::vector<std::unordered_set<std::size_t>> &FirstLevelNodePatch, const std::unordered_set<std::size_t> &NodePatch) {
    std::unordered_set<std::size_t> InternalNode;
    for (const auto& ConnectedNodeInPatch : NodePatch) {
        if (FirstLevelNodePatch[ConnectedNodeInPatch].size() >= OrderOfPatch)
            InternalNode.insert(ConnectedNodeInPatch);
    }
    return InternalNode;
}

class Gradient_Recovery : public RecoveryPatchConnManager {
private:
    std::vector<double> Dev_x, Dev_y, Dev_xx, Dev_xy, Dev_yy;
    std::vector<double> longest_length;

private:
    void computeLongestLengthByFirstLevelPatch(const std::vector<double> &x, const std::vector<double> &y);

public:
    void generateGradientRecovery(const std::vector<double> &x, const std::vector<double> &y, const std::vector<double> &solution);
    std::vector<double> get_dev_x() const {return Dev_x;}
    std::vector<double> get_dev_y() const {return Dev_y;}
    std::vector<double> get_dev_xx() const {return Dev_xx;}
    std::vector<double> get_dev_xy() const {return Dev_xy;}
    std::vector<double> get_dev_yy() const {return Dev_yy;}

};

void Gradient_Recovery::generateGradientRecovery(const std::vector<double> &x, const std::vector<double> &y, const std::vector<double> &solution) {
    computeLongestLengthByFirstLevelPatch(x, y);
    for (std::size_t i = 0; i < number_nodes; ++i) {
        double max_length = longest_length[i];
        std::size_t SizeOfFullLevelPatch = FullLevelNodePatch[i].size() + 1;
        std::vector<double> xi, eta, local_solution;
        xi.push_back(0.0), eta.push_back(0.0), local_solution.push_back(solution[i]);
        for (std::size_t j : FullLevelNodePatch[i]) {
            xi.push_back((x[j] - x[i])/max_length);
            eta.push_back((y[j] - y[i])/max_length);
            local_solution.push_back(solution[j]);
        }
        std::vector<std::vector<double>> A(SizeOfFullLevelPatch, std::vector<double>(OrderOfPatch, 0.0));
        for (std::size_t j = 0; j < SizeOfFullLevelPatch; ++j) {
            A[j][0] = 1.0;
            A[j][1] = xi[j], A[j][2] = eta[j];
            A[j][3] = xi[j] * xi[j], A[j][4] = xi[j] * eta[j], A[j][5] = eta[j] * eta[j];
        }
        std::vector<std::vector<double>> A_squared = computeTransposeProduct(A);
        std::vector<double> B = computeTransposeVectorProduct(A, local_solution);
        std::vector<double> local_derivative;
        DenseSystemSolver ::solveDenseMatrixSystem(A_squared, B, local_derivative);
        Dev_x.push_back(local_derivative[1]/max_length);
        Dev_y.push_back(local_derivative[2]/max_length);
        Dev_xx.push_back(2 * local_derivative[3]/(max_length*max_length));
        Dev_xy.push_back(local_derivative[4]/(max_length*max_length));
        Dev_yy.push_back(2 * local_derivative[5]/(max_length*max_length));
    }
}

void Gradient_Recovery::computeLongestLengthByFirstLevelPatch(const std::vector<double> &x, const std::vector<double> &y) {
    for (std::size_t i = 0; i < number_nodes; ++i) {
        std::vector<double> x_, y_;
        x_.push_back(x[i]), y_.push_back(y[i]);
        for (std::size_t j : FirstLevelNodePatch[i])
            x_.push_back(x[j]), y_.push_back(y[j]);
        
        double max_length = std::sqrt((x_[0] - x_[1]) * (x_[0] - x_[1]) + (y_[0] - y_[1]) * (y_[0] - y_[1]));
        for (std::size_t j = 2; j < x_.size(); ++j) {
            double length = std::sqrt((x_[0] - x_[j]) * (x_[0] - x_[j]) + (y_[0] - y_[j]) * (y_[0] - y_[j]));
            if (length > max_length) max_length = length;
        }
        longest_length.push_back(max_length);
    }
}

// class Gradient_Interpolation : public Gradient_Recovery {
// private:
//     std::size_t getIndexOfBelongingElem(double x, double y);
//     bool isPointInTriangle(double x, double y);

// public:
//     Gradient_Interpolation() = default;
//     std::tuple<double, double> interpolateDev(double x, double y);
// };

// bool Gradient_Interpolation::isPointInTriangle(double z, double r) {
//     double z1 = Coord[TriElem[0]].first, r1 = Coord[TriElem[0]].second;
//     double z2 = Coord[TriElem[1]].first, r2 = Coord[TriElem[1]].second;
//     double z3 = Coord[TriElem[2]].first, r3 = Coord[TriElem[2]].second;

//     double Denominator = (r2 - r3) * (z1 - z3) + (z3 - z2) * (r1 - r3);

//     double a = ((r2 - r3) * (z - z3) + (z3 - z2) * (r - r3)) / Denominator;
//     double b = ((r3 - r1) * (z - z3) + (z1 - z3) * (r - r3)) / Denominator;
//     double c = 1.0 - a - b;

//     if (a >= 0 && a <= 1 && b >= 0 && b <= 1 && c >= 0 && c <= 1) {
//         return true;
//     } else {
//         return false;
//     }
// }

// void Gradient_Interpolation::interpolateDev(double z, double r) {
//     Ind IndElem = getIndexOfBelongingElem(z, r);
//     Ind IndFirstNode = ElemConnData[IndElem][0]; 
//     Ind IndSecondNode = ElemConnData[IndElem][1]; 
//     Ind IndThirdNode = ElemConnData[IndElem][2]; 
    
//     double z1 = Coord[IndFirstNode].first;
//     double r1 = Coord[IndFirstNode].second;
//     double z2 = Coord[IndSecondNode].first;
//     double r2 = Coord[IndSecondNode].second;
//     double z3 = Coord[IndThirdNode].first;
//     double r3 = Coord[IndThirdNode].second;

//     double Area = AreaPerElem[IndElem];

//     double FirstCoefficient = ((z2 * r3 - z3 * r2) + (r2 - r3) * z + (z3 - z2) * r) / (2 * Area);
//     double SecondCoefficient = ((z3 * r1 - z1 * r3) + (r3 - r1) * z + (z1 - z3) * r) / (2 * Area);
//     double ThirdCoefficient = 1.0 - FirstCoefficient - SecondCoefficient;

//     double zDev = FirstCoefficient * RecoveredDev[IndFirstNode][0]
//                     + SecondCoefficient * RecoveredDev[IndSecondNode][0]
//                     + ThirdCoefficient * RecoveredDev[IndThirdNode][0];
//     double rDev = FirstCoefficient * RecoveredDev[IndFirstNode][1]
//                     + SecondCoefficient * RecoveredDev[IndSecondNode][1]
//                     + ThirdCoefficient * RecoveredDev[IndThirdNode][1];
//     double zzDev = FirstCoefficient * RecoveredDev[IndFirstNode][2]
//                     + SecondCoefficient * RecoveredDev[IndSecondNode][2]
//                     + ThirdCoefficient * RecoveredDev[IndThirdNode][2];
//     double zrDev = FirstCoefficient * RecoveredDev[IndFirstNode][3]
//                     + SecondCoefficient * RecoveredDev[IndSecondNode][3]
//                     + ThirdCoefficient * RecoveredDev[IndThirdNode][3];
//     double rrDev = FirstCoefficient * RecoveredDev[IndFirstNode][4]
//                     + SecondCoefficient * RecoveredDev[IndSecondNode][4]
//                     + ThirdCoefficient * RecoveredDev[IndThirdNode][4];
//     InterpolatedDev = {zDev, rDev, zzDev, zrDev, rrDev};
// }

// void Gradient_Interpolation::calcAreaElem() {
//     for (const auto & Elem : ElemConnData) {
//         Ind IndFirstNode = Elem[0]; 
//         Ind IndSecondNode = Elem[1];
//         Ind IndThirdNode = Elem[2];
//         double z1 = Coord[IndFirstNode].first;
//         double r1 = Coord[IndFirstNode].second;
//         double z2 = Coord[IndSecondNode].first;
//         double r2 = Coord[IndSecondNode].second;
//         double z3 = Coord[IndThirdNode].first;
//         double r3 = Coord[IndThirdNode].second;
//         double Area = calcTriangularArea(z1, r1, z2, r2, z3, r3);
//         AreaPerElem.push_back(Area);
//     }
// }

// std::size_t Gradient_Interpolation::getIndexOfBelongingElem(double x, double y) {
//     for (Ind IndElem = 0; IndElem < ElemConnData.size(); ++IndElem) {
//         TriConn TriElem = ElemConnData[IndElem];
//         if (isPointInTriangle(x, y, TriElem)) {
//             return IndElem;
//         } else {
//             continue;
//         }
//     }
//     return INVALID_ELEMENT_INDEX;
// }

#endif