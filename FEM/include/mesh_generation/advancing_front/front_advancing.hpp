#ifndef MESH_GENERATION_FRONT_ADVANCING_HPP
#define MESH_GENERATION_FRONT_ADVANCING_HPP

#include <geometry/geometry.hpp>
#include <mesh_generation/background_mesh.hpp>

#include <vector>
#include <algorithm>
#include <unordered_map>
#include <map>

namespace mesh_generation
{

    template<unsigned int _spdim, typename real_t = double, typename index_t = unsigned int>
    class front_advancing {
        private:
            background_mesh_dat<_spdim, real_t> background_mesh;
            
            std::vector<geometry::edge<_spdim, real_t>> raw_boundary;
            std::vector<geometry::edge<_spdim, real_t>> discretized_boundary;
            std::vector<geometry::edge<_spdim, real_t>> active_front;

            void discretize_boundary(std::vector<index_t> idx_edge_to_stretch, std::vector<index_t> n_per_edge, std::vector<real_t> Q) {
                if (idx_edge_to_stretch.size() != raw_boundary.size())
                    throw std::invalid_argument("the number of stretching edges exceeds the number of raw edges to discretize.");
                if (n_per_edge.size() != raw_boundary.size())
                    throw std::invalid_argument("the number of discretization parameters does not match the number of raw edges.");
                if (n_per_edge.size() != Q.size())
                    throw std::invalid_argument("the number of discretization parameters does not match the number of raw edges.");
            
                index_t n_e = raw_boundary.size();
                
                for (index_t i = 0; i < n_e; ++i) {
                    const geometry::edge<_spdim, real_t>& _edge = raw_boundary[i];
                    geometry::point<_spdim, real_t> A = _edge.first();
                    geometry::point<_spdim, real_t> B = _edge.second();                   
                    
                    geometry::point<_spdim, real_t> current_point = A;
                    index_t n_e = n_per_edge[i];

                    for (index_t j = 0; j < n_e; ++j) {
                        real_t t = 0.0;
                        real_t xi = static_cast<real_t>(j) / static_cast<real_t>(n_e);
                        if (idx_edge_to_stretch[i] != 0) {
                            t = (std::tanh(Q[i] * (xi - 0.5)) - std::tanh(-0.5 * Q[i])) / (2.0 * std::tanh(0.5 * Q[i]));
                        } else {
                            t = xi;
                        }

                        geometry::point<_spdim, real_t> new_point = A + (B - A) * t;
                        discretized_boundary.push_back(geometry::edge<_spdim, real_t>(current_point, new_point));
                    }
                }
            }
            
            [[nodiscard]] geometry::edge<_spdim, real_t> find_shortest_edge(index_t& idx) {
                if (active_front.empty())
                    throw std::runtime_error("Active front is empty. No edges to search.");
            
                auto it_min = active_front.begin();
                real_t min_length = it_min->length();
                
                index_t current_idx = 0;
                idx = 0;
            
                for (auto it = active_front.begin(); it != active_front.end(); ++it, ++current_idx) {
                    real_t length = it->length();
                    if (length < min_length) {
                        min_length = length;
                        it_min = it;
                        idx = current_idx;
                    }
                }
                return *it_min;
            }
            
            [[nodiscard]] geometry::point<_spdim, real_t> transform_point_to_normalized_space(
                const geometry::point<_spdim, real_t>& P,
                const data_structure::matrix<_spdim, _spdim, real_t, index_t>& T
            )
            {
                data_structure::vector<_spdim, real_t> vP = T * data_structure::vector<_spdim, real_t>(A.operator());
                geometry::point<_spdim, real_t>& result(vP.get_data());
                return result;
            }


        public:
            front_advancing() = default;
            ~front_advancing() = default;

            // front_advancing(const geometry::polygon<_spdim, real_t> &domain) : raw_boundary(domain.get_boundary_edge()) {}

            void initializebackground_mesh(const std::vector<real_t> &_xc,
                                            const std::vector<real_t> &_D,
                                            const std::vector<real_t> &_delta,
                                            const std::vector<geometry::point<_spdim, real_t>> &_S)
            {
                background_mesh.initialize_background_mesh(_xc, _D, _delta, _S);
            }

            void generate_mesh(std::vector<index_t> idx_edge_to_stretch, std::vector<index_t> n_per_edge, std::vector<real_t> Q) {
                discretize_boundary(idx_edge_to_stretch, n_per_edge, Q);
                active_front = discretized_boundary;
                index_t num_active_edges = active_front.size();
                // step 1, search the shortest edge of the active front
                index_t idx_shortest_edge;
                geometry::edge<_spdim, real_t> AB = find_shortest_edge(idx_shortest_edge);
                geometry::point<_spdim, real_t> A = AB.first();
                geometry::point<_spdim, real_t> B = AB.second();
                // step 2 center of the side M interpolate T and A_hat, B_hat, M_hat
                geometry::point<_spdim, real_t> M = geometry::point<_spdim, real_t>::mid_point(AB.first(), AB.second());
                data_structure::matrix<_spdim, _spdim, real_t, index_t> T = background_mesh.compute_local_transformation(M);
                data_structure::vector<_spdim, real_t> vA = T * data_structure::vector<_spdim, real_t>(A);
                data_structure::vector<_spdim, real_t> vB = T * data_structure::vector<_spdim, real_t>(A);
                data_structure::vector<_spdim, real_t> vM = T * data_structure::vector<_spdim, real_t>(M);
                geometry::point<_spdim, real_t> A_hat = transform_point_to_normalized_space(A, T);
                geometry::point<_spdim, real_t> B_hat = transform_point_to_normalized_space(B, T);
                geometry::point<_spdim, real_t> M_hat = transform_point_to_normalized_space(M, T);

                // select all relevant edge
                std::vector<geometry::edge<_spdim, real_t>> relevant_edges;
                real_t LAB = AB.length();
                for (index_t i = idx_shortest_edge; i < num_active_edges; ++i) {
                    const geometry::point<_spdim, real_t> P1 = active_front[i].first();
                    const geometry::point<_spdim, real_t> P2 = active_front[i].second();
                    if (M.dist(P1) - 3 * LAB < ETOL<real_t> && M.dist(P2) - 3 * LAB < ETOL<real_t>)
                        relevant_edges.push_back(active_front[i]);
                }

                //search backward
                for (index_t i = idx_shortest_edge; i < num_active_edges; ++i) {
                    const geometry::point<_spdim, real_t> P1 = active_front[i].first();
                    const geometry::point<_spdim, real_t> P2 = active_front[i].second();
                    if (M.dist(P1) - 3 * LAB < ETOL<real_t> && M.dist(P2) - 3 * LAB < ETOL<real_t>)
                        relevant_edges.push_back(active_front[i]);
                }

                // step 3, Determine, in the normalized space, the ideal position P̂ 1 for the vertex of the triangular element.
                // The point P̂ 1 is located on the line perpendicular to the side that passes through the point M̂
                // and at a distance δ1 from the points Â and B̂ . The direction in which P̂ 1 is generated is determined
                // by the orientation of the side. The value δ 1 is chosen according to the criterion where L is the distance between points Â and B̂ . Only in situations where the side AB happens to have characteristics very different from those specified by the background mesh will the value of δ 1 be different from unity. However, the above inequalities must be taken into account to ensure geometrical compatibility. Expression (3) is pure empirical, and different inequalities could be devised to serve the same purpose.

                real_t L = A_hat.dist(B_hat);
                real_t delta_1 = 0.0;
                if (0.55 * L - 1 <= ETOL<real_t> && 1 - 2.00 * L < ETOL<real_t>)
                    delta_1 = 1.0;
                else if (0.55 * L - 1 <= ETOL<real_t>)
                    delta_1 = 0.55 * L;
                else
                    delta_1 = 2.00 * L;
                
                geometry::edge<_spdim, real_t> AB_hat(A_hat, B_hat);
                geometry::vector<_spdim, real_t> vP1_hat = AB_hat.perpendicular();
                geometry::point<_spdim, real_t> P1_hat = M_hat + (vP1_hat * delta_1);

                // step 4
                std::map<real_t, geometry::point<_spdim, real_t>, std::greater<real_t>> ideal_point_maps;
                geometry::point<_spdim, real_t> P2_hat = M_hat + (vP1_hat * (delta_1/5 * 4));
                geometry::point<_spdim, real_t> P3_hat = M_hat + (vP1_hat * (delta_1/5 * 3));
                geometry::point<_spdim, real_t> P4_hat = M_hat + (vP1_hat * (delta_1/5 * 2));
                geometry::point<_spdim, real_t> P5_hat = M_hat + (vP1_hat * (delta_1/5 * 4));

                for (const geometry::edge<_spdim, real_t> & front : active_front) {
                    const geometry::point<_spdim, real_t>& Q1 = front.first();
                    const geometry::point<_spdim, real_t>& Q2 = front.second();
                    geometry::point<_spdim, real_t> Q1_hat = transform_point_to_normalized_space(Q1, T);
                    geometry::point<_spdim, real_t> Q2_hat = transform_point_to_normalized_space(Q2, T);

                    if (Q1_hat != A_hat && Q1_hat != B_hat && Q1_hat.dist(P1_hat) < ETOL<real_t>) {
                        point<_spdim, real_t> C1_hat = geometry::circle::circumcircle(Q1_hat, A_hat, B_hat).get_center();
                        real_t dist_C1_to_P1 = C1_hat.dist(P1_hat);
                        ideal_point_maps.insert(dist_C1_to_P1, Q1_hat);
                    }
                    if (Q2_hat != A_hat && Q2_hat != B_hat && Q2_hat.dist(P2_hat) < ETOL<real_t>) {
                        point<_spdim, real_t> C2_hat = geometry::circle::circumcircle(Q2_hat, A_hat, B_hat).get_center();
                        real_t dist_C2_to_P2 = C2_hat.dist(P1_hat);
                        ideal_point_maps.insert(dist_C2_to_P2, Q1_hat);
                    }
                }
                std::vector<geometry::point<_spdim, real_t>> ideal_point_arr(ideal_point_maps.begin(), ideal_point_maps.end());
                ideal_point_arr.push_back(P1_hat);
                ideal_point_arr.push_back(P2_hat);
                ideal_point_arr.push_back(P3_hat);
                ideal_point_arr.push_back(P4_hat);
                ideal_point_arr.push_back(P5_hat);

                // step 5 and step 6
                geometry::point<_spdim, real_t> new_point;
                for (const geometry::point<_spdim, real_t>& I : ideal_point_arr) {
                    geometry::edge<_spdim, real_t> AI_hat(A_hat, I_hat);
                    geometry::edge<_spdim, real_t> BI_hat(B_hat, I_hat);
                    
                    for (const geometry::edge<_spdim, real_t> relevant_edge : relevant_edges) {
                        geometry::edge<_spdim, real_t> relevant_edge_hat;
                        if (AI_hat.is_intersect_with(relevant_edge_hat) && BI_hat.is_intersect_with(relevant_edge_hat))
                            continue;
                        else {
                            data_structure::matrix<_spdim, _spdim, real_t, index_t> Tinverse = T.inverse();
                            data_structure::vector<_spdim, real_t> vI(I.operator());
                            data_structure::vector<_spdim, real_t> vnewpoint = Tinverse * vI;
                            new_point(vnewpoint.operator());
                        }
                    }
                }

                // step 7: Store the new triangle and update the front by adding/removing the relevant sides.
                

            }

            [[nodiscard]] const std::vector<geometry::edge<_spdim, real_t>>& get_discretized_boundary() const {
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