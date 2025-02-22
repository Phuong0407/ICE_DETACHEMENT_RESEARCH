/**
 * @file vertex.hpp
 * @author Phuong Diep (diepthanhphuong0407@gmail.com)
 * @brief derived class of mesh_entity, which represents a vertex in mesh topology.
 * @version 0.1
 * @date 2025-02-16
 * 
 * @copyright Copyright (c) 2025
 * 
 */

#ifndef VERTEX_HPP
#define VERTEX_HPP

#include "mesh_entity.hpp"

/**
 * @brief 
 * 
 */
namespace mesh_entity
{

    /**
     * @class vertex
     * @brief A derived class of mesh_entity that represents a vertex in the mesh topology.
     * 
     * @tparam spdim space dimension (must be 2 or 3).
     * @tparam index_t integer type used for indexing the vertex (e.g., unsigned int, size_t).
     * 
     * @details
     *  Since a vertex has no sub-entities, its geometric dimension is always 0.
     *  It inherits from `mesh_entity<spdim, 0, index_t>`, where:
     *      - `spdim` specifies the spatial dimension of the mesh.
     *      - `0` represents the geometric dimension of a vertex.
     *      - `index_t` defines the type of the vertex index.
     * 
     * @remark Key Features:
     * - A vertex has no sub-entities (`nse = 0`).
     * - Each vertex is uniquely identified by its index (`id`).
     * - The `print()` function provides basic information about the vertex.
     * 
     * @remark Public Member Functions:
     * - `vertex() = default;` → Default constructor.
     * - `explicit vertex(index_t id);` → Constructs a vertex with a given unique index.
     * - `void print() const override;` → Prints vertex information.
     * 
     * @remark Example Usage:
     * @code
     * vertex<2, unsigned int> v(5);
     * v.print(); // Output: vertex, id = 5.
     * @endcode
     * 
     */

    template<unsigned int spdim, typename index_t>
    class vertex : public mesh_entity<spdim, 0, index_t> {
        public:
            vertex() = default;
            
            explicit vertex(index_t id) : mesh_entity<spdim, 0, index_t>(id, 0) {}

            void print() const override {
                std::cout << "vertex, id = " << this->id() << "." << std::endl;
            }
    };

} // namespace mesh_entity

#endif