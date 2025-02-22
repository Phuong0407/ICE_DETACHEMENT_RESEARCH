#include <mesh_generation/background_mesh.hpp>
#include <mesh_generation/front_advancing.hpp>

#include <geometry/point.hpp>
#include <geometry/vector.hpp>

#include <vector>
#include <array>

#define _spdim 2
#define real_t double
#define index_t unsigned int

int main() {
    
    geometry::point<_spdim, real_t> p1(0.0, 0.0);
    geometry::point<_spdim, real_t> p2(0.0, 1.0);
    geometry::point<_spdim, real_t> p3(1.0, 0.0);
    geometry::point<_spdim, real_t> p4(1.0, 1.0);

    std::vector<geometry::point<_spdim, real_t>> vertices;
    vertices.push_back(p1);
    vertices.push_back(p2);
    vertices.push_back(p3);
    vertices.push_back(p4);

    std::vector<std::array<index_t, _spdim + 1>> connections;
    connections.push_back({0, 1, 2});
    connections.push_back({0, 1, 3});

    std::array<real_t, _spdim> delta_1 = {0.5, 1.0};
    std::array<real_t, _spdim> delta_2 = {0.2, 0.8};
    std::array<real_t, _spdim> delta_3 = {0.1, 1.3};
    std::array<real_t, _spdim> delta_4 = {1.5, 0.15};

    geometry::vector<_spdim, real_t> v11(1.0, 0.0);
    geometry::vector<_spdim, real_t> v12(0.0, 1.0);
    geometry::vector<_spdim, real_t> v21(1.0, 0.0);
    geometry::vector<_spdim, real_t> v22(0.0, 1.0);
    geometry::vector<_spdim, real_t> v31(1.3, 2.0);
    geometry::vector<_spdim, real_t> v32(2.9, 1.0);
    geometry::vector<_spdim, real_t> v41(1.0, 0.0);
    geometry::vector<_spdim, real_t> v42(0.0, 1.0);

    std::array<geometry::vector<_spdim, real_t>, _spdim> alpha1 = {v11, v12};
    std::array<geometry::vector<_spdim, real_t>, _spdim> alpha2 = {v21, v22};
    std::array<geometry::vector<_spdim, real_t>, _spdim> alpha3 = {v31, v32};
    std::array<geometry::vector<_spdim, real_t>, _spdim> alpha4 = {v41, v42};

    mesh_generation::control_point_t<_spdim, real_t> c1(delta_1, alpha1);
    mesh_generation::control_point_t<_spdim, real_t> c2(delta_2, alpha2);
    mesh_generation::control_point_t<_spdim, real_t> c3(delta_3, alpha3);
    mesh_generation::control_point_t<_spdim, real_t> c4(delta_4, alpha4);

    std::vector<mesh_generation::control_point_t<_spdim, real_t>> control_point_datas;
    control_point_datas.push_back(c1);
    control_point_datas.push_back(c2);
    control_point_datas.push_back(c3);
    control_point_datas.push_back(c4);

    // mesh_generation::front_advancing a;
    // (vertices, connections, control_point_datas);


}