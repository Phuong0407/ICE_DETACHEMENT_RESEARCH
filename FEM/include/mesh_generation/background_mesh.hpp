#ifndef MESH_GENERATION_BACKGROUND_MESH_HPP
#define MESH_GENERATION_BACKGROUND_MESH_HPP

#include "geometry.hpp"
#include "point.hpp"
#include "edge.hpp"
#include "polygon.hpp"
#include "triangle.hpp"
#include "matrix.hpp"

namespace mesh_generation
{
    
    template<unsigned int _spdim, typename double_t = double>
    class background_mesh {
    private:
        std::vector<point<_spdim, double_t>> vertices;
        std::vector<std::array<unsigned int, _spdim>> connections;
        std::vector<vector<_spdim, double_t>> area;
        std::vector<double_t> delta;
        std::vector<double_t> xc;
        std::vector<double_t> D;
        std::vector<vector<_spdim, double_t>> alpha;
        std::vector<matrix<_spdim, double_t>> T;

    public:
        background_mesh() = default;
        ~background_mesh() = default;

        void _compute_transformation_for_control_point() {
            T.clear();
            unsigned int num_control_points = vertices.size();
            for (unsigned int i = 0; i < num_control_points; ++i) {
                matrix<_spdim, double_t> _T;
                for (unsigned int j = 0; j < _spdim; ++j) {
                    for (unsigned int k = 0; k < _spdim; ++k) {
                        for (unsigned int l = 0; l < _spdim; ++l) {
                            _T(j, k) += alpha[j][l] * alpha[l][j] / delta[i];
                        }
                    }
                }
                T.push_back(_T);
            }
        }

        void _compute_cell_areas() {
            unsigned int num_cells = connections.size();
            area.clear();
            area.resize(num_cells);
            for (unsigned int i = 0; i < num_cells; ++i) {
                unsigned int id1 = connections[i][0];
                unsigned int id2 = connections[i][1];
                unsigned int id3 = connections[i][2];

                point<_spdim, double_t> p1 = vertices[id1];
                point<_spdim, double_t> p2 = vertices[id2];
                point<_spdim, double_t> p3 = vertices[id3];

                if constexpr (_spdim == 2) {
                    double_t cell_area = 0.5 * std::abs(
                        (p2(0) - p1(0)) * (p3(1) - p1(1)) - 
                        (p3(0) - p1(0)) * (p2(1) - p1(1))
                    );
                    area[i] = vector<_spdim, double_t>({cell_area});
                }
            }
        }

        [[nodiscard]] unsigned int _idx_cell(point<_spdim, double_t> P) {
            unsigned int num_cells = connections.size();
            for (unsigned int i = 0; i < num_cells; ++i) {
                unsigned int id1 = connections[i][0];
                unsigned int id2 = connections[i][1];
                unsigned int id3 = connections[i][2];
                const point<_spdim, double_t>& p1 = vertices[id1];
                const point<_spdim, double_t>& p2 = vertices[id2];
                const point<_spdim, double_t>& p3 = vertices[id3];

                double_t denominator = (p2.x[1] - p3.x[1]) * (p1.x[0] - p3.x[0]) +
                                    (p3.x[0] - p2.x[0]) * (p1.x[1] - p3.x[1]);

                double_t w1 = ((p2.x[1] - p3.x[1]) * (P.x[0] - p3.x[0]) +
                            (p3.x[0] - p2.x[0]) * (P.x[1] - p3.x[1])) / denominator;
                
                double_t w2 = ((p3.x[1] - p1.x[1]) * (P.x[0] - p3.x[0]) +
                            (p1.x[0] - p3.x[0]) * (P.x[1] - p3.x[1])) / denominator;
                
                double_t w3 = 1.0 - w1 - w2;

                if (w1 >= 0 && w2 >= 0 && w3 >= 0) {
                    return i;
                }
            }
            return std::numeric_limits<unsigned int>::max();
        }

        [[nodiscard]] matrix<2, double_t> interpolate_transformation(const point<2, double_t>& P) const {
            // Step 1: Find the element that contains P
            unsigned int element_idx = _idx_cell(P);
            if (element_idx == std::numeric_limits<unsigned int>::max()) {
                throw std::runtime_error("Point is outside the background mesh.");
            }
        
            // Step 2: Retrieve the element's nodes and transformation matrices
            std::array<unsigned int, 3> node_ids = connections[element_idx]; // Triangle: 3 nodes
            const point<2, double_t>& A = vertices[node_ids[0]];
            const point<2, double_t>& B = vertices[node_ids[1]];
            const point<2, double_t>& C = vertices[node_ids[2]];
        
            const matrix<2, double_t>& T_A = T[node_ids[0]];
            const matrix<2, double_t>& T_B = T[node_ids[1]];
            const matrix<2, double_t>& T_C = T[node_ids[2]];
        
            // Step 3: Compute Shape Function Values at P (Barycentric Coordinates)
            double_t area = 0.5 * ((B.x[0] - A.x[0]) * (C.x[1] - A.x[1]) - (C.x[0] - A.x[0]) * (B.x[1] - A.x[1]));
        
            if (std::abs(area) < std::numeric_limits<double_t>::epsilon()) {
                throw std::runtime_error("Degenerate triangle encountered.");
            }
        
            double_t N1 = ((B.x[1] - C.x[1]) * (P.x[0] - C.x[0]) + (C.x[0] - B.x[0]) * (P.x[1] - C.x[1])) / (2.0 * area);
            double_t N2 = ((C.x[1] - A.x[1]) * (P.x[0] - C.x[0]) + (A.x[0] - C.x[0]) * (P.x[1] - C.x[1])) / (2.0 * area);
            double_t N3 = 1.0 - N1 - N2;
        
            // Step 4: Interpolate T_b at P
            matrix<2, double_t> Tb;
            Tb.set_zero();
            Tb += N1 * T_A;
            Tb += N2 * T_B;
            Tb += N3 * T_C;
        
            return Tb;
        }
        
        std::pair<std::array<double_t, 2>, vector<_spdim, double_t>> compute_eigen_2x2(const matrix<2, double_t>& A) {
            T a = A[0][0], b = A[0][1];
            T c = A[1][0], d = A[1][1];
        
            T trace = a + d;
            T determinant = a * d - b * c;
            T discriminant = std::sqrt(trace * trace - 4 * determinant);
        
            T lambda1 = (trace + discriminant) / 2.0;
            T lambda2 = (trace - discriminant) / 2.0;
        
            std::array<T, 2> eigenvalues = {lambda1, lambda2};
        
            std::array<std::array<T, 2>, 2> eigenvectors;
            for (int i = 0; i < 2; ++i) {
                T lambda = eigenvalues[i];
        
                // Solve (A - lambda I) v = 0
                if (std::abs(b) > std::numeric_limits<T>::epsilon()) {
                    eigenvectors[i] = {b, lambda - a};  // v = (b, lambda-a)
                } else if (std::abs(c) > std::numeric_limits<T>::epsilon()) {
                    eigenvectors[i] = {lambda - d, c};  // v = (lambda-d, c)
                } else {
                    eigenvectors[i] = {1.0, 0.0};  // Default eigenvector if matrix is diagonal
                }
        
                // Normalize eigenvector
                T norm = std::sqrt(eigenvectors[i][0] * eigenvectors[i][0] + eigenvectors[i][1] * eigenvectors[i][1]);
                eigenvectors[i][0] /= norm;
                eigenvectors[i][1] /= norm;
            }
        
            return {eigenvalues, eigenvectors};
        }

        // [[nodiscard]] matrix<_spdim, double_t> interpolate_transformation(point<_spdim, double_t> P) const {
        // }

        // void initialize_background_mesh(const std::vector<double_t>& xc,
        //                                 const std::vector<double_t>& D,
        //                                 const std::vector<double_t>& delta,
        //                                 const std::vector<point<_spdim, double_t>>& S,
        //                                 const std::vector<vector<_spdim, double_t>>& alpha,
        //                                 const std::vector<vector<_spdim, double_t>>& scaling)
        // {
        //     _S = S;
        //     _D = D;
        //     _xc = xc;
        //     target_element_size = delta;
        //     _alpha = alpha;
        //     directional_scaling = scaling;
        // }
    };


} // namespace mesh_generatio


#endif