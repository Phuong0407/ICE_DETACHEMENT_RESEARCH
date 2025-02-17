/**
 * @file mesh_entity.hpp
 * @author Phuong Diep (diepthanhphuong0407@gmail.com)
 * @brief base class, which represents entity in mesh topology, deploy a hierarchical structure.
 * @version 0.1
 * @date 2025-02-16
 * 
 * @copyright Copyright (c) 2025
 * 
 */

#ifndef MESH_ENTITY_HPP
#define MESH_ENTITY_HPP

#include <iostream>
#include <type_traits>

/**
 * @brief 
 * 
 */
namespace mesh_entity
{

    /**
     * @class mesh_entity
     * @brief base class representing a generic mesh entity with basic properties and functionality.
     * 
     * @tparam _spdim space dimension (must be 2 or 3).
     * @tparam _gedim geometric dimension of the entity:
     *              - 0 for vertices
     *              - 1 for edges
     *              - _spdim - 1 for faces, acepted only for 3D mesh
     *                      - supported only triangles and quadrilaterals
     *              - _spdim for cells
     *                      - triangles and quadrilaterals in 2D
     *                      - tetrahedra and hexahedra in 3D
     * @tparam index_t integer type used for indexing mesh entities
     *              (unsigned int, size_t, or custom integer type, etc.).
     * 
     * @details
     *  This class provides a foundation for representing mesh entities in a hierarchical structure.
     *  Each entity is uniquely identified by its index and stores the number of its subentities.
     * 
     * @remark Compile-Time constraints:
     *  - `_spdim` must be 2 or 3 (only 2D and 3D meshes are supported).
     *  - `_gedim` must be less than or equal to `_spdim`
     *              (an entity cannot have a higher geometric dimension than the space it resides in).
     * - `index_t` must be an integral type.
     * 
     * @remark Memory and Performance Considerations:
     *  - Copy operations are explicitly disabled to prevent unintended deep copies.
     *  - Move operations are enabled to allow efficient resource transfer.
     * 
     * @remark Required Implementations for Derived Classes:
     * - `nv()`: Returns the number of vertices in the entity (not required for vertices).
     * - `print()`: Outputs entity-specific information.
     * - `order()`: Computes the number of vertices for faces and cells.
     * 
     * @remark Public Member Functions:
     *  - `id()`: Returns the unique identifier of the mesh entity.
     *  - `nse()`: Returns the number of subentities.
     *  - `sp_dim()`: Returns the spatial dimension (constant for all entities in the same mesh).
     *  - `ge_dim()`: Returns the geometric dimension of the entity.
     * 
     * @remarks Usage Notes:
     *  - This class serves as a base for more specific mesh entity types such as vertices, edges, faces, and cells.
     *  - The `order()` function is implemented only for faces and cells to compute the number of vertices.
     *  - The `nv()` function is implemented for edges, faces, and cells but always returns 0 for vertices. 
    */

    template<unsigned int _spdim, unsigned int _gedim, typename index_t>
    class mesh_entity {
        private:
            static_assert(_spdim == 2 || _spdim == 3, "error: only 2D and 3D are supported.");
            static_assert(_gedim <= _spdim, "error: a mesh entity dimension is less than or equal to space dimension.");            
            static_assert(std::is_integral_v<index_t>, "index_t must be an integral type");

            index_t _id;
            unsigned int _nse;

        public:
            mesh_entity() = default;
            explicit mesh_entity(index_t _id, unsigned int _nse) : _id(_id) , _nse(_nse) {}

            mesh_entity(const mesh_entity&) = delete;
            mesh_entity& operator=(const mesh_entity&) = delete;
            mesh_entity(mesh_entity&&) noexcept = default;
            mesh_entity& operator=(mesh_entity&&) noexcept = default;

            index_t id() const { return _id; }
            unsigned int nse() const { return _nse; }
            static constexpr unsigned int sp_dim() { return _spdim; }
            static constexpr unsigned int ge_dim() { return _gedim; }

            virtual unsigned int order() const { return 0; }
            virtual unsigned int nv() { return 0; }
            virtual void print() const = 0;

            virtual ~mesh_entity() = default;
    };

} // namespace mesh_entity

#endif