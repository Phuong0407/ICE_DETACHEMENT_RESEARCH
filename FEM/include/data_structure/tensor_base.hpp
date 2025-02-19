#ifndef DATA_STRUCTURE_TENSOR_BASE_HPP
#define DATA_STRUCTURE_TENSOR_BASE_HPP

#include <iostream>
#include <vector>
#include <array>
#include <memory>
#include <numeric>
#include <stdexcept>
#include <exception>

namespace data_structure
{

template<unsigned int Rank, unsigned int Dim, typename double_t = double>
class tensor_base {
public:
    using ts_ptr = std::unique_ptr<tensor_base<Rank, Dim, double_t>>;

    tensor_base() = default;
    virtual ~tensor_base() = default;

    virtual void init() = 0;
    virtual void init(const double_t* raw_data, size_t size) = 0;
    virtual void init(const std::vector<double_t>& vec) = 0;

    template<size_t N>
    virtual void init(const std::array<double_t, N>& arr) = 0;

    virtual void init(std::initializer_list<double_t> list) = 0;

    virtual ts_ptr operator+(const tensor_base& other) const = 0;
    virtual ts_ptr operator-(const tensor_base& other) const = 0;
    virtual ts_ptr operator*(const tensor_base& other) const = 0;

    virtual const std::vector<double_t>& get_data() const = 0;

    template<unsigned int Rank2, unsigned int Dim2>
    virtual std::unique_ptr<tensor_base<Rank + Rank2, Dim, double_t>> outer(const tensor_base<Rank2, Dim2, double_t>& other) const = 0;

    template<unsigned int Rank2, unsigned int Dim2>
    virtual std::unique_ptr<tensor_base<Rank + Rank2 - 2, Dim, double_t>> dot(const tensor_base<Rank2, Dim2, double_t>& other) const = 0;

    virtual void print() const = 0;
};

} // namespace data_structure

#endif