#ifndef MESH_GENERATION_FRONT_ADVANCING_HPP
#define MESH_GENERATION_FRONT_ADVANCING_HPP

#include "geometry.hpp"
#include "point.hpp"
#include "edge.hpp"
#include "triangle.hpp"
#include "polygon.hpp"
#include "matrix.hpp"

#include <vector>
#include <algorithm>
#include <unordered_map>

namespace mesh_generation
{

    template<unsigned int _spdim, typename double_t = double, typename index_t = unsigned int>
    class front_advancing {
        private:
            background_mesh_controller<_spdim, double_t> background_mesh;
            
            std::vector<edge<_spdim, double_t>> raw_boundary;
            std::vector<edge<_spdim, double_t>> discretized_boundary;
            std::vector<edge<_spdim, double_t>> active_front;

            void _discretize_boundary(std::vector<index_t> idx_edge_to_stretch, std::vector<index_t> n_per_edge, std::vector<double_t> Q) {
                if (idx_edge_to_stretch.size() != raw_boundary.size())
                    throw std::invalid_argument("The number of stretching edges exceeds the number of raw edges to discretize.");
                if (n_per_edge.size() != raw_boundary.size())
                    throw std::invalid_argument("The number of discretization parameters does not match the number of raw edges.");
                if (n_per_edge.size() != Q.size())
                    throw std::invalid_argument("The number of discretization parameters does not match the number of raw edges.");
            
                index_t n_e = raw_boundary.size();
                
                for (index_t i = 0; i < n_e; ++i) {
                    const edge<_spdim, double_t>& _edge = raw_boundary[i];
                    point<_spdim, double_t> A = _edge.first();
                    point<_spdim, double_t> B = _edge.second();                   
                    
                    point<_spdim, double_t> current_point = A;
                    index_t n_e = n_per_edge[i];

                    for (index_t j = 0; j < n_e; ++j) {
                        double_t t = 0.0;
                        double_t xi = static_cast<double_t>(j) / static_cast<double_t>(n_e);
                        if (idx_edge_to_stretch[i] != 0)
                            t = (std::tanh(Q[i] * (xi - 0.5)) - std::tanh(-0.5 * Q[i])) / (2.0 * std::tanh(0.5 * Q[i]));
                        else
                            t = xi;

                        point<_spdim, double_t> new_point = A + (B - A) * t;
                        discretized_boundary.push_back(edge<_spdim, double_t>(current_point, new_point));
                    }
                }
            }
            
            [[nodiscard]] edge<_spdim, double_t> _find_shortest_edge() {
                if (active_front.empty())
                    throw std::runtime_error("Active front is empty. No edges to search.");
            
                edge<_spdim, double_t> min_edge = active_front.begin();
                double_t min_length = min_edge.length();
            
                for (auto it = active_front.begin(); it != active_front.end(); ++it) {
                    double_t length = it.length();
                    if (length < min_length) {
                        min_length = length;
                        min_edge = it;
                    }
                }
                return *min_edge;
            }
            

        public:
            front_advancing() = default;
            ~front_advancing() = default;

            front_advancing(const polygon<_spdim, double_t> &domain) : raw_boundary(domain.get_boundary_edge()) {}

            void initialize_background_mesh(const std::vector<double_t> &_xc,
                                            const std::vector<double_t> &_D,
                                            const std::vector<double_t> &_delta,
                                            const std::vector<point<_spdim, double_t>> &_S)
            {
                background_mesh.initialize_background_mesh(_xc, _D, _delta, _S);
            }

            void generate_mesh(std::vector<index_t> idx_edge_to_stretch, std::vector<index_t> n_per_edge, std::vector<double_t> Q) {
                _discretize_boundary(idx_edge_to_stretch, n_per_edge, Q);
                active_front = discretized_boundary;

                // step 1, search the shortest edge of the active front
                edge<_spdim, double_t> AB = _find_shortest_edge();

                // step 2
                // center of the side M
                point<_spdim, double_t> M = point<_spdim, double_t>::mid_point(AB.first(), AB.second());
                // interpolate T

                // step 3, Determine, in the normalized space, the ideal position P̂ 1 for the vertex of the triangular element.
                // The point P̂ 1 is located on the line perpendicular to the side that passes through the point M̂
                // and at a distance δ1 from the points Â and B̂ . The direction in which P̂ 1 is generated is determined
                // by the orientation of the side. The value δ 1 is chosen according to the criterion
                matrix<_spdim, double_t> T = background_mesh.compute_transformation(M);
                vector<_spdim, double_t> vA = T * vector<_spim, double_t>(A);
                vector<_spdim, double_t> vB = T * vector<_spim, double_t>(B);
                vector<_spdim, double_t> vM = T * vector<_spim, double_t>(M);


            }

            [[nodiscard]] const std::vector<edge<_spdim, double_t>>& get_discretized_boundary() const {
                return discretized_boundary;
            }

            void print_discretized_boundary () const {
                for (const auto& _edge : discretized_boundary) { 
                    _edge.print();
                    std::cout << std::endl;
                }
            }
    };

}

#endif