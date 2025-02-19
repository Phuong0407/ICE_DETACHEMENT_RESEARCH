#ifndef MESH_GENERATION_ADVANCING_FRONT_HPP
#define MESH_GENERATION_ADVANCING_FRONT_HPP

#include "point.hpp"
#include "triangle.hpp"
#include "polygon.hpp"

#include <vector>
#include <cmath>
#include <algorithm>
#include <limits>

namespace mesh_generation
{

    template<unsigned int _spdim, typename double_t = double>
    class advancing_front {
    private:
        std::vector<point<_spdim, double_t>> active_front;
        std::vector<point<_spdim, double_t>> control_points;
        std::vector<triangle<_spdim, double_t>> mesh;
        polygon<_spdim, double_t> polygon_boundary;
        double_t node_spacing = 1.0;
        double_t stretching_param = 1.0;

        bool is_point_inside_polygon(const point<_spdim, double_t>& P) const {
            return polygon_boundary.contains(P);
        }

        point<_spdim, double_t> find_new_triangle_point(const point<_spdim, double_t>& A, 
                                                     const point<_spdim, double_t>& B) 
        {
            point<_spdim, double_t> M = (A + B) * 0.5;
            double_t d = std::max(0.55 * A.dist(B), std::min(node_spacing, 2 * A.dist(B)));

            double_t dx = B[1] - A[1];
            double_t dy = A[0] - B[0];
            double_t norm = std::sqrt(dx * dx + dy * dy);
            dx = (dx / norm) * d;
            dy = (dy / norm) * d;

            point<_spdim, double_t> C1 = {M[0] + dx, M[1] + dy};
            point<_spdim, double_t> C2 = {M[0] - dx, M[1] - dy};

            if (is_point_inside_polygon(C1)) return C1;
            if (is_point_inside_polygon(C2)) return C2;

            return polygon_boundary.nearest_point_on_boundary(M);
        }

        bool __valid_triangle__(const point<_spdim, double_t>& A, 
                                const point<_spdim, double_t>& B, 
                                const point<_spdim, double_t>& C) const 
        {
            double_t area = 0.5 * std::abs(A[0] * (B[1] - C[1]) + B[0] * (C[1] - A[1]) + C[0] * (A[1] - B[1]));
            return area > std::numeric_limits<double_t>::epsilon();
        }

        void update_front(const point<_spdim, double_t>& A, 
                         const point<_spdim, double_t>& B, 
                         const point<_spdim, double_t>& C) 
        {
            remove_edge_from_front(A, B);
            active_front.push_back(C);
        }

        void remove_edge_from_front(const point<_spdim, double_t>& A, 
                                 const point<_spdim, double_t>& B) 
        {
            auto it = std::find_if(active_front.begin(), active_front.end(),
                [&](const edge<_spdim, double_t>& e) {
                    return (e.first() == A && e.second() == B) || (e.first() == B && e.second() == A);
                });
            
            if (it != active_front.end()) {
                active_front.erase(it);
            }
        }

    public:
        advancing_front() = default;
        advancing_front(const polygon<_spdim, double_t>& domain) : polygon_boundary(domain) {}

        void initialize_mesh() {
            active_front = polygon_boundary.get_boundary_edge();
        }

        void generate_mesh() {
            while (!active_front.empty()) {
                auto edge = active_front.back();
                active_front.pop_back();
                point<_spdim, double_t> A = edge[0];
                point<_spdim, double_t> B = edge[1];

                point<_spdim, double_t> C = find_new_triangle_point(A, B);

                if (__valid_triangle__(A, B, C)) {
                    mesh.emplace_back(A, B, C);
                    update_front(A, B, C);
                }
            }
        }

        const std::vector<triangle<_spdim, double_t>>& getMesh() const {
            return mesh;
        }
    };

} // namespace mesh_generation

#endif