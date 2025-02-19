// #ifndef DATA_STRUCTURE_TENSOR_BASE_HPP
// #define DATA_STRUCTURE_TENSOR_BASE_HPP

// #include <iostream>
// #include <vector>
// #include <array>
// #include <memory>
// #include <numeric>
// #include <stdexcept>
// #include <exception>

// namespace data_structure
// {

//     template<unsigned int Rank, unsigned int Dim, typename real_t = double, typename index_t = unsigned int>
//     class tensor_base {
//     public:
//         using ts_ptr = std::unique_ptr<tensor_base<Rank, Dim, real_t>>;

//         tensor_base() = default;
//         virtual ~tensor_base() = default;

//         virtual void init() = 0;
//         virtual void init(const real_t* raw_data, index_t size) = 0;
//         virtual void init(const std::vector<real_t>& vec) = 0;
//         virtual void init(std::initializer_list<real_t> list) = 0;

//         virtual ts_ptr operator+(const tensor_base& other) const = 0;
//         virtual ts_ptr operator-(const tensor_base& other) const = 0;
//         virtual ts_ptr operator*(const tensor_base& other) const = 0;

//         virtual const std::vector<real_t>& get_data() const = 0;
//         virtual void set_data(index_t idx, real_t val) = 0;
//         virtual unsigned int data_size() const = 0;

//         virtual void print() const = 0;
//     };

// } // namespace data_structure

// #endif

#ifndef DATA_STRUCTURE_TENSOR_BASE_HPP
#define DATA_STRUCTURE_TENSOR_BASE_HPP

#include <iostream>
#include <vector>
#include <array>
#include <memory>
#include <numeric>
#include <stdexcept>
#include <exception>

namespace data_structure {

    template<unsigned int Rank, unsigned int Dim, typename real_t = double, typename index_t = unsigned int>
    class tensor_base {
        public:
            using ts_ptr = std::unique_ptr<tensor_base<Rank, Dim, real_t, index_t>>;

            tensor_base() = default;
            virtual ~tensor_base() = default;
            virtual const std::vector<real_t>& get_data() const = 0;
            virtual void set_data(const std::vector<real_t>& new_data) = 0;
    };

} // namespace data_structure

#endif
