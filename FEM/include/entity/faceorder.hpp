#ifndef FACEORDER_HPP
#define FACEORDER_HPP

#include <string>

namespace mesh_entity
{
    enum class faceorder_t {
        zero,
        linear,
        quadratic,
        cubic,
        quartic,
        quintic,
    };

    class faceorder_helper {
        public:
            static constexpr unsigned int nth(faceorder_t order) {
                if (order == faceorder_t::linear)
                    return 1;
                else if (order == faceorder_t::quadratic)
                    return 2;
                else if (order == faceorder_t::cubic)
                    return 3;
                else if (order == faceorder_t::quartic)
                    return 4;
                else if (order == faceorder_t::quintic)
                    return 5;
                else
                    return 0;
            }

            static std::string to_string(faceorder_t order) {
                if (order == faceorder_t::linear)
                    return "linear";
                else if (order == faceorder_t::quadratic)
                    return "quadratic";
                else if (order == faceorder_t::cubic)
                    return "cubic";
                else if (order == faceorder_t::quartic)
                    return "quartic";
                else if (order == faceorder_t::quintic)
                    return "quintic";
                else
                    return "zero";
            }
    };

} // namespace mesh_entity

#endif