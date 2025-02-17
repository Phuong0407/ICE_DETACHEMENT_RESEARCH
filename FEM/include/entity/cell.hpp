#ifndef CELL_HPP
#define CELL_HPP

#include "mesh_entity.hpp"
#include "cell_helper.hpp"

#include <algorithm>
#include <stdexcept>

namespace mesh_entity {

    template<unsigned int _spdim, typename index_t>
    class cell : public mesh_entity<_spdim, _spdim, index_t> {
    private:
        celltype_t _t = celltype_t::undefined;
        cellorder_t _o = cellorder_t::zero;

    public:
        cell() = default;

        explicit cell(index_t id, celltype_t _t, cellorder_t _o) :
            _t(_t), _o(_o), mesh_entity<_spdim, _spdim, index_t>(id, cell_helper::__nse(_spdim, _t))
        {
            if (!cell_helper::__valid_celltype(_spdim, _t)) {
                throw std::invalid_argument("error: cannot initialize a "
                        + std::to_string(cell_helper::__spdim(_t)) + "D cell"
                        " in " + std::to_string(_spdim) + "D space.");
            }
        }

    public:
        celltype_t cell_type() const { return _t; }

        unsigned int order() const override {
            return cellorder_helper::nth(_o);
        }

        unsigned int nv() override {
            return cell_helper::__nv(_spdim, _t, _o);
        }

        void print() const override {
            std::cout << "cell, id = " << this->id() << ", ";
            std::cout << "type = " << celltype_helper::to_string(_t) << ", ";
            std::cout << "order = " << this->order() << "." << std::endl;
        }

    };

} // namespace mesh_entity

#endif