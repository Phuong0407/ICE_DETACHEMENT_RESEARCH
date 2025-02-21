#ifndef DATA_STRUCTURE_TENSOR_TENSOR_BASE_HPP
#define DATA_STRUCTURE_TENSOR_TENSOR_BASE_HPP

#include <iostream>
#include <vector>
#include <array>
#include <memory>
#include <stdexcept>
#include <initializer_list>

namespace data_structure {

    template<unsigned int Rank, typename real_t = double, typename index_t = unsigned int>
    class tensor_base {
    protected:
        std::array<index_t, Rank> _shape;
        std::vector<real_t> _data;

        index_t idx(const std::array<index_t, Rank>& ids) const {
            index_t index = 0, stride = 1;
            for (index_t i = Rank; i-- > 0;) {
                if (ids[i] >= _shape[i])
                    throw std::out_of_range("index out of bounds in inx().");
                index += ids[i] * stride;
                stride *= _shape[i];
            }
            return index;
        }

    public:
        using ts_ptr = std::unique_ptr<tensor_base<Rank, real_t, index_t>>;

        tensor_base() = default;
        virtual ~tensor_base() = default;

        virtual const std::array<index_t, Rank>& get_shape() const = 0;
        virtual void init() = 0;
        virtual void init(const real_t* raw_data, index_t size) = 0;
        virtual void init(const std::vector<real_t>& vec) = 0;
        virtual void init(std::initializer_list<real_t> list) = 0;

        virtual ts_ptr operator+(const tensor_base& other) const = 0;
        virtual ts_ptr operator-(const tensor_base& other) const = 0;
        virtual ts_ptr operator*(const tensor_base& other) const = 0;

        virtual ts_ptr contract(index_t i, index_t j) = 0;

        unsigned int data_size() const { return _data.size(); }

        const std::vector<real_t>& get_data() const { return _data; }

        void set_data(const std::vector<real_t>& new_data) {
            if (new_data.size() != _data.size()) {
                throw std::invalid_argument("mismatched data size in set_data().");
            }
            _data = new_data;
        }

        void set_data(index_t idx, real_t val) {
            if (idx >= _data.size())
                throw std::out_of_range("index out of bounds in set_data()");
            _data[idx] = val;
        }

        virtual void print() const override {
            for (const auto& val : _data) {
                std::cout << val << " ";
            }
            std::cout << std::endl;
        }
    };

}

#endif