#ifndef CONFIG_H
#define CONFIG_H

// DEFINE THE PRECISION TO BE USED
template<typename double_t>
constexpr double_t ETOL =
#ifdef USE_HYPER_PRECISION
    static_cast<double_t>(1e-15);
#elif defined(USE_ULTRA_PRECISION)
    static_cast<double_t>(1e-12);
#elif defined(USE_HIGH_PRECISION)
    static_cast<double_t>(1e-6);
#else
    static_cast<double_t>(1e-12);
#endif // DEFINE THE PRECISION TO BE USED

#endif