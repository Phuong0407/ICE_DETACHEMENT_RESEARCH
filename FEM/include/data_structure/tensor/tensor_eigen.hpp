#pragma one
#include "tensor_base.hpp"

#define USE_EIGEN
#ifdef USE_EIGEN

#include <Eigen/Dense>

namespace data_structure
{
    
    template<unsigned int rank, unsigned int dim, typename double_t = double>
    class tensor : public tensor_base<rank, dim> {
        private:
            Eigen::Matrix<double_t, Eigen::Dynamic, Eigen::Dynamic> tensor;

        public:

        [[nodiscard]] tensor<rank, dim, double_t> & void add(const tensor<rank, dim, double_t> &other) {
            return tensor + other.tensor;
        }
    };

} // namespace data_structure


#endif