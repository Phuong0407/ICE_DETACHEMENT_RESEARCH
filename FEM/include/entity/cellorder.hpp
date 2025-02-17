#ifndef CELLORDER_HPP
#define CELLORDER_HPP

#include <string>

namespace mesh_entity {
    enum class cellorder_t {
        zero,
        linear,
        quadratic,
        cubic,
        quartic,
        quintic,
    };

    class cellorder_helper {
    public:
        static constexpr unsigned int nth(cellorder_t order) {
            if (order == cellorder_t::linear)
                return 1;
            else if (order == cellorder_t::quadratic)
                return 2;
            else if (order == cellorder_t::cubic)
                return 3;
            else if (order == cellorder_t::quartic)
                return 4;
            else if (order == cellorder_t::quintic)
                return 5;
            else
                return 0;
        }

        static std::string to_string(cellorder_t order) {
            if (order == cellorder_t::linear)
                return "linear";
            else if (order == cellorder_t::quadratic)
                return "quadratic";
            else if (order == cellorder_t::cubic)
                return "cubic";
            else if (order == cellorder_t::quartic)
                return "quartic";
            else if (order == cellorder_t::quintic)
                return "quintic";
            else
                return "zero";
        }
    };

} // namespace mesh_entity

#endif