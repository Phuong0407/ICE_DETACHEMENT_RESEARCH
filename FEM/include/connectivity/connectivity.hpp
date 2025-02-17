/**
 * @file connectivity.hpp
 * @author Phuong Diep (diepthanhphuong0407@gmail.com)
 * @brief The connectivity between mesh entities. MULTI-ELEMENT SUPPORTED.
 * @version 0.1
 * @date 2025-02-16
 * 
 * @copyright Copyright (c) 2025
 * 
 */

#ifndef CONNECTIVITY_HPP
#define CONNECTIVITY_HPP

#include <vector>
#include <stdexcept>
#include <type_traits>

namespace mesh_connectivity
{
    template<template<typename ...> class mesh_entity1,
             template<typename ...> class mesh_entity2,
             typename index_t>
    class mesh_connectivity {
        private:
            static_assert(std::is_integral_v<index_t>, "index_t must be an integer type.");
            std::vector<mesh_entity1<_spdim, _index_t>> _h;
            std::vector<mesh_entity2<_spdim, _index_t>> _g;
            std::vector<index_t> indices;
            std::vector<index_t> offsets;

        public:
            mesh_connectivity() = default;

            void init_mesh_entity(
                std::vector<mesh_entity1<_spdim, _index_t>> &_h,
                std::vector<mesh_entity2<_spdim, _index_t>> &_g)
            {
                this->_h = std::move(_h);
                this->_g = std::move(_g);    
            }

            void init_indices(std::vector<unsigned int>& indices) {
                this->indices = std::move(indices);
                // TO DO
                // init the set of offsets used by _h._nse for the offsets (multi-element), use for loops
                _h[i]._nse();
            }
            
            void print() const {}
};

} // namespace mesh_connectivity


#endif