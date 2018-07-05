// template <typename saidx_t>
// struct libdivsufsort_config;
//
// template <>
// struct libdivsufsort_config<int32_t>
// {
//     static constexpr uint64_t TR_STACKSIZE = 64;
//     static constexpr uint64_t TR_INSERTIONSORT_THRESHOLD = 8;
//     static constexpr uint64_t SS_SMERGE_STACKSIZE = 32;
// };
//
// template <>
// struct libdivsufsort_config<int64_t>
// {
//     static constexpr uint64_t TR_STACKSIZE = 96;
//     static constexpr uint64_t TR_INSERTIONSORT_THRESHOLD = 8;
//     static constexpr uint64_t SS_SMERGE_STACKSIZE = 64;
// };

#define SS_SMERGE_STACKSIZE (32)
#define TR_STACKSIZE (64)

#include <sdsl/divsufsort_impl.hpp>
