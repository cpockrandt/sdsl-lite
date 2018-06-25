#include <sdsl/divsufsort_config.hpp>

/* minstacksize = log(SS_BLOCKSIZE) / log(3) * 2 */
#if SS_BLOCKSIZE == 0
// # if defined(BUILD_DIVSUFSORT64)
// #  define SS_MISORT_STACKSIZE (96)
// # else
#  define SS_MISORT_STACKSIZE (64)
// # endif
#elif SS_BLOCKSIZE <= 4096
# define SS_MISORT_STACKSIZE (16)
#else
# define SS_MISORT_STACKSIZE (24)
#endif
// #if defined(BUILD_DIVSUFSORT64)
// # define SS_SMERGE_STACKSIZE (64)
// #else
# define SS_SMERGE_STACKSIZE (32)
// #endif

#define TR_INSERTIONSORT_THRESHOLD (8)
// #if defined(BUILD_DIVSUFSORT64)
// # define TR_STACKSIZE (96)
// #else
# define TR_STACKSIZE (64)
// #endif

#include <sdsl/divsufsort_impl.hpp>
