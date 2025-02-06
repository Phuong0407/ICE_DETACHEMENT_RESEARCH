#ifndef GRADIENT_RECOVERY_H
#define GRADIENT_RECOVERY_H

#include "data_structures_fem.h"

#include <unordered_set>
#include <vector>

constexpr std::size_t order_of_patch = 6;





inline std::vector<std::unordered_set<std::size_t>> generate_first_level_patch(const triangle_element* element_connection, std::size_t num_elements, std::size_t num_nodes) {
    std::vector<std::unordered_set<std::size_t>> first_level_patch(num_nodes);

    for (std::size_t i = 0; i < num_elements; ++i) {
        std::size_t idx_elem_1 = element_connection[i].node_inds[0];
        std::size_t idx_elem_2 = element_connection[i].node_inds[1];
        std::size_t idx_elem_3 = element_connection[i].node_inds[2];

        first_level_patch[idx_elem_1].insert({idx_elem_2, idx_elem_3});
        first_level_patch[idx_elem_2].insert({idx_elem_3, idx_elem_1});
        first_level_patch[idx_elem_3].insert({idx_elem_1, idx_elem_2});
    }
    return first_level_patch;
}



inline std::unordered_set<std::size_t> find_all_internal_nodes(const std::vector<std::unordered_set<std::size_t>> &first_level_patch, const std::unordered_set<std::size_t> &patch_of_node) {
    std::unordered_set<std::size_t> internal_nodes;
    for (const auto& connected_node : patch_of_node)
        if (first_level_patch[connected_node].size() >= order_of_patch)
            internal_nodes.insert(connected_node);
    return internal_nodes;
}





std::vector<std::unordered_set<std::size_t>> generate_full_level_patch(const triangle_element* element_connection, std::size_t num_elements, std::size_t num_nodes, std::vector<std::unordered_set<std::size_t>> &first_level_patch) {
    first_level_patch.clear();
    first_level_patch = generate_first_level_patch(element_connection, num_elements, num_nodes);
    std::vector<std::unordered_set<std::size_t>> full_level_patch = first_level_patch;

    for (std::size_t i = 0; i < num_nodes; ++i) {
        const auto &node_patch = first_level_patch[i];
        if (node_patch.size() >= order_of_patch) continue;

        std::unordered_set<std::size_t> internal_nodes = find_all_internal_nodes(first_level_patch, node_patch);
        for (const auto &internal_node : internal_nodes) {
            full_level_patch[i].insert(first_level_patch[internal_node].begin(), first_level_patch[internal_node].end());
        } full_level_patch[i].erase(i);
    }

    for (std::size_t i = 0; i < num_nodes; ++i) {
        while (full_level_patch[i].size() < order_of_patch) {
            for (const auto &neighbor : full_level_patch[i]) {
                full_level_patch[i].insert(full_level_patch[neighbor].begin(), full_level_patch[neighbor].end());
            }
        } full_level_patch[i].erase(i);
    }
    return full_level_patch;
}



#endif