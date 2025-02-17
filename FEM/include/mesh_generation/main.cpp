#include "point.hpp"
#include "edge.hpp"
#include "vector.hpp"
#include "line.hpp"

#include <iostream>

#define _spdim 2

int main() {
    mesh_generation::point<_spdim> p(1.0, 2.0);
    mesh_generation::point<_spdim> q(2.0, 4.0);
    mesh_generation::point<_spdim> r(0.0, 0.0);
    mesh_generation::edge<_spdim> e(p, q);
    mesh_generation::point<_spdim> f = r.project_to_edge(e);
    f.print();
    return 0;
}