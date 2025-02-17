#ifndef EDGE_HPP
#define EDGE_HPP

#include "mesh_entity.hpp"

namespace mesh_entity {

    template<unsigned int _spdim, typename index_t>
    class edge : public mesh_entity<_spdim, 1, index_t> {
    public:
        edge() = default;

        explicit edge(index_t id) : mesh_entity<_spdim, 1, index_t>(id, 2) {}

        unsigned int nv() override {
            return 2;
        }

        void print() const override {
            std::cout << "edge, id = " << this->id() << "." << std::endl;
        }
};

} // namespace mesh_entity

#endif