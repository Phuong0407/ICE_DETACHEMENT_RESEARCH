#ifndef MESH_GENERATION_BACKGROUND_MESH_DAT_HPP
#define MESH_GENERATION_BACKGROUND_MESH_DAT_HPP

#include <geometry/geometry.hpp>
#include <data_structure/matrix/matrix.hpp>
#include <data_structure/matrix/vector.hpp>

#include <Eigen/Dense>

#include <array>
#include <cmath>
#include <vector>
#include <limits>
#include <iostream>
#include <algorithm>

namespace mesh_generation
{

    template<unsigned int _spdim, typename real_t = double>
    class control_point_t {
        public:
            std::array<real_t, _spdim> delta;
            std::array<geometry::vector<_spdim, real_t>, _spdim> alpha;

            control_point_t(
                const std::array<real_t, _spdim>& delta,
                const std::array<geometry::vector<_spdim, real_t>, _spdim>& alpha
            ) : delta(delta), alpha(alpha) {}

    };

    template<typename real_t = double>
    class distribution_point_t {
        public:
            real_t delta;
            real_t xc;
            real_t D;

            distribution_point_t(real_t delta, real_t xc, real_t D) : delta(delta), xc(xc), D(D) {}
            ~distribution_point_t() = default;
    };


    template<unsigned int _spdim, typename real_t = double, typename index_t = unsigned int>
    class background_mesh_dat {
        private:
            std::vector<geometry::point<_spdim, real_t>> vertices;
            std::vector<std::array<unsigned int, _spdim + 1>> connections;
            std::vector<control_point_t<_spdim, real_t>> control_datas;
            std::vector<data_structure::matrix<_spdim, _spdim, real_t, index_t>> transformations;

            std::vector<geometry::point<_spdim, real_t>> distribution_points;
            std::vector<distribution_point_t<real_t>> distribution_dats;

            void compute_control_point_transformation() {
                transformations.clear();
                index_t num_vertices = vertices.size();
                for (index_t i = 0; i < num_vertices; ++i) {
                    data_structure::matrix<_spdim, _spdim, real_t, index_t> transformation;
                    for (unsigned int j = 0; j < _spdim; ++j) {
                        real_t delta = control_datas[i].delta[j];
                        geometry::vector<_spdim, real_t> alpha = control_datas[i].alpha[j];
                        data_structure::vector<_spdim, real_t, index_t> alpha_vec(alpha.operator()());
                        transformation += alpha_vec.outer(alpha_vec) * (1 / delta);
                    }
                    transformations.push_back(transformation);
                }
            }

            const index_t id_element(const geometry::point<_spdim, real_t> P, real_t *beta1, real_t *beta2, real_t *beta3) const {
                if constexpr(_spdim == 2) {
                    index_t num_elements = connections.size();
                    for (index_t i = 0; i < num_elements; ++i) {
                        const auto& element = connections[i];

                        index_t id1 = element[0];
                        index_t id2 = element[1];
                        index_t id3 = element[2];

                        real_t x = P(0), y = P(1);

                        real_t x1 = vertices[id1](0), y1 = vertices[id1](1);
                        real_t x2 = vertices[id2](0), y2 = vertices[id2](1);
                        real_t x3 = vertices[id3](0), y3 = vertices[id3](1);
                        
                        real_t det_T = (x2 - x1) * (y3 - y1) - (x3 - x1) * (y2 - y1);

                        (*beta1) = ((y2 - y3) * (x - x3) + (x3 - x2) * (y - y3)) / det_T;
                        (*beta2) = ((y3 - y1) * (x - x3) + (x1 - x3) * (y - y3)) / det_T;
                        (*beta3) = 1.0 - beta1 - beta2;
                        
                        if (beta1 >= 0 && beta1 <= 1 && beta2 >= 0 && beta2 <= 1 && beta3 >= 0 && beta3 <= 1)
                            return i;
                    }
                    return static_cast<index_t>(-1);
                }
            }

            const real_t compute_local_spacing(const geometry::point<_spdim, real_t> P) const {
                real_t local_spacing = 0.0;
                index_t num_distri_points = distribution_points.size();
                for (unsigned int i = 0; i < num_distri_points; ++i) {
                    const geometry::point<_spdim, real_t> & Q = distribution_points[i];
                    real_t x = P.dist(Q);
                    real_t delta_x = 0.0;
                    const distribution_point_t<real_t> dat = distribution_dats[i];
                    real_t xc = dat.xc;
                    if (x - xc <= ETOL<real_t>)
                        delta_x = dat.delta;
                    else {
                        const real_t D = dat.D;
                        delta_x = dat.delta * std::exp(std::abs((x - xc)/(D - xc)) * std::log(2.0));
                    }
                    local_spacing = std::min(local_spacing, delta_x);
                }
                return local_spacing;
            }

        public:
            background_mesh_dat() = default;
            ~background_mesh_dat() = default;

            background_mesh_dat(
                std::vector<geometry::point<_spdim, real_t>> vertices,
                std::vector<std::array<unsigned int, _spdim + 1>> connections,
                std::vector<control_point_t<_spdim, real_t>> control_datas
            ) : vertices(vertices), connections(connections), control_datas(control_datas)
            {
                compute_control_point_transformation();
            }

            const data_structure::matrix<_spdim, _spdim, real_t, index_t>& compute_local_transformation(const geometry::point<_spdim, real_t> &P) const {
                real_t beta1 = 0.0, beta2 = 0.0, beta3 = 0.0;
                index_t id_elem = id_element(P, &beta1, &beta2, &beta3);
                if (id_elem == static_cast<index_t>(-1))
                    throw std::runtime_error("point is outside the background mesh.");

                const auto& element = connections[id_elem];
                index_t id1 = element[0];
                index_t id2 = element[1];
                index_t id3 = element[2];
                
                data_structure::matrix<_spdim, _spdim, real_t, index_t> Tb;

                #ifdef USING_OPENMP
                    #pragma omp parallel for collapse(2)
                #endif
                for (unsigned int i = 0; i < _spdim; ++i) {
                    for (unsigned int j = 0; j < _spdim; ++j) {
                        Tb(i, j) = beta1 * transformations[id1](i, j) +
                                   beta2 * transformations[id2](i, j) +
                                   beta3 * transformations[id3](i, j);
                    }
                }

                Eigen::MatrixXd Tb_eigen(_spdim, _spdim);
                for (unsigned int i = 0; i < _spdim; ++i) {
                    for (unsigned int j = 0; j < _spdim; ++j) {
                        Tb_eigen(i, j) = Tb(i, j);
                    }
                }

                Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> solver(Tb_eigen);
                if (solver.info() != Eigen::Success) {
                    throw std::runtime_error("eigenvalue decomposition failed.");
                }
                
                data_structure::matrix<_spdim, _spdim, real_t, index_t> T;
                Eigen::VectorXd eigenvalues = solver.eigenvalues();
                Eigen::MatrixXd eigenvectors = solver.eigenvectors();
                Eigen::VectorXd delta_star(_spdim);
                for (unsigned int i = 0; i < _spdim; ++i) {
                    real_t local_spacing = compute_local_spacing(P);
                    delta_star(i) = std::min(eigenvalues(i), local_spacing);
                }
                                
                #ifdef USING_OPENMP
                    #pragma omp parallel for collapse(2)
                #endif
                for (unsigned int i = 0; i < _spdim; ++i) {
                    for (unsigned int j = 0; j < _spdim; ++j) {
                        for (unsigned int k = 0; k < _spdim; ++k) {
                            if (eigenvalues(k) > 0) {
                                T(i, j) += (1.0 / eigenvalues(k)) * eigenvectors(i, k) * eigenvectors(j, k);
                            }
                        }
                    }
                }
                return T;
            }
            
    };


} // namespace mesh_generation

#endif