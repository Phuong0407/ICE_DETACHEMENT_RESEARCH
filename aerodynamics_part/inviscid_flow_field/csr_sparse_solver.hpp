#ifndef CSR_SPARSE_SOLVER_H
#define CSR_SPARSE_SOLVER_H

#include <cmath>
#include <vector>
#include <numeric>
#include <iostream>
#include <algorithm>
#include <stdexcept>
#include <limits>

class CSRSparseSolver {
private:
    double ZeroComparedTol;
    double MatNorm;
    long unsigned int SizeMat;
    long unsigned int NonZeroElem;

    std::vector<double> Val;
    std::vector<long unsigned int> ColInd;
    std::vector<long unsigned int> RowPtr;

    const double CONV_TOL = 1.0E-6;
    const int MAX_ITER = 100000;
    void multSparseMat(const std::vector<double>& X, std::vector<double>& Y) const;

public:
    CSRSparseSolver() = default;
    void initSparseSolver(const std::vector<std::vector<double>>& SPARSE_MAT);
    void solveSparseSystem(const std::vector<double>& RHS_VECT, std::vector<double>& X) const;
};

void CSRSparseSolver::initSparseSolver(const std::vector<std::vector<double>>& SPARSE_MAT) {
    SizeMat = SPARSE_MAT.size();
    MatNorm = 0.0;

    for (const auto& row : SPARSE_MAT) {
        MatNorm = std::max(MatNorm, std::accumulate(row.begin(), row.end(), 0.0, [](double sum, double Val) {
            return sum + std::abs(Val);
        }));
    }

    ZeroComparedTol = std::numeric_limits<double>::epsilon() * MatNorm;

    NonZeroElem = 0;
    for (const auto& row : SPARSE_MAT) {
        for (double Val : row) {
            if (std::abs(Val) > ZeroComparedTol) {
                ++NonZeroElem;
            }
        }
    }

    Val.resize(NonZeroElem);
    ColInd.resize(NonZeroElem);
    RowPtr.resize(SizeMat + 1);

    int index = 0;
    RowPtr[0] = 0;

    for (std::size_t i = 0; i < SizeMat; ++i) {
        for (std::size_t j = 0; j < SizeMat; ++j) {
            if (std::abs(SPARSE_MAT[i][j]) > ZeroComparedTol) {
                Val[index] = SPARSE_MAT[i][j];
                ColInd[index] = j;
                ++index;
            }
        }
        RowPtr[i + 1] = index;
    }
}

void CSRSparseSolver::multSparseMat(const std::vector<double>& X, std::vector<double>& Y) const {
    if (X.size() != SizeMat) {
        throw std::invalid_argument("VECTOR SIZE DOES NOT MATCH MATRIX SIZE.");
    }

    if (Y.size() != SizeMat) {
        Y.resize(SizeMat, 0.0);
    } else {
        std::fill(Y.begin(), Y.end(), 0.0);
    }

    for (std::size_t i = 0; i < SizeMat; ++i) {
        for (std::size_t j = RowPtr[i]; j < RowPtr[i + 1]; ++j) {
            Y[i] += Val[j] * X[ColInd[j]];
        }
    }
}

void CSRSparseSolver::solveSparseSystem(const std::vector<double>& RHS_VECT, std::vector<double>& X) const {
    if (RHS_VECT.size() != SizeMat) {
        throw std::invalid_argument("VECTOR SIZE DOES NOT MATCH MATRIX SIZE.");
    }

    if (X.size() != SizeMat) {
        X.resize(SizeMat, 0.0);
    }
    else {
        std::fill(X.begin(), X.end(), 0.0);
    }

    std::vector<double> r(SizeMat, 0.0), p(SizeMat, 0.0), Ap(SizeMat, 0.0);
    double alpha, beta, rnorm, bnorm, tol = CONV_TOL;

    bnorm = std::sqrt(std::inner_product(RHS_VECT.begin(), RHS_VECT.end(), RHS_VECT.begin(), 0.0));
    r = RHS_VECT;
    p = r;
    rnorm = std::inner_product(r.begin(), r.end(), r.begin(), 0.0);

    int ITER = 0;
    while (std::sqrt(rnorm) / bnorm > tol && ITER < MAX_ITER) {
        ++ITER;
        multSparseMat(p, Ap);
        alpha = rnorm / std::inner_product(p.begin(), p.end(), Ap.begin(), 0.0);

        for (std::size_t i = 0; i < SizeMat; ++i) {
            X[i] += alpha * p[i];
            r[i] -= alpha * Ap[i];
        }

        if (std::sqrt(std::inner_product(r.begin(), r.end(), r.begin(), 0.0)) / bnorm < tol) {
            break;
        }

        double new_rnorm = std::inner_product(r.begin(), r.end(), r.begin(), 0.0);
        beta = new_rnorm / rnorm;
        rnorm = new_rnorm;

        for (std::size_t i = 0; i < SizeMat; ++i) {
            p[i] = r[i] + beta * p[i];
        }
    }

    if (ITER == MAX_ITER) {
        std::cerr << "FAILED TO CONVERGE IN " << MAX_ITER << " ITERATIONS.\n";
    }
}

#endif