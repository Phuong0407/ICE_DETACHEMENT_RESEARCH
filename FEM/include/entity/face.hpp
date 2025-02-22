#ifndef FACE_HPP
#define FACE_HPP

#include "mesh_entity.hpp"
#include "face_helper.hpp"

namespace mesh_entity {

    template<unsigned int spdim, typename index_t>
    class face : public mesh_entity<spdim, spdim - 1, index_t> {
        private:
            facetype_t _t;
            faceorder_t _o;
        public:
            face() = default;

            explicit face(index_t id, facetype_t _t, faceorder_t _o) :
                mesh_entity<spdim, spdim - 1, index_t>(id, face_helper::__nse(_t)), _t(_t), _o(_o) {}

            unsigned int order() const override {
                return faceorder_helper::nth(_o);
            }

            unsigned int nv() override {
                return face_helper::__nv(spdim - 1, _t, _o);
            }

            void print() const override {
                std::cout << "face, id = " << this->id() << ", ";
                std::cout << "type = " << facetype_helper::to_string(_t) << ", ";
                std::cout << "order = " << this->order() << "." << std::endl;
            }
    };

} // namespace mesh_entity

#endif
