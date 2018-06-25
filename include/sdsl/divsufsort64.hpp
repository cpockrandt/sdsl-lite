// #define BUILD_DIVSUFSORT64

#include <sdsl/divsufsort_config.hpp>

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#ifndef _DIVSUFSORT64_H
#define _DIVSUFSORT64_H 1

#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */

#include <inttypes.h>

#ifndef DIVSUFSORT_API
# ifdef DIVSUFSORT_BUILD_DLL
#  define DIVSUFSORT_API
# else
#  define DIVSUFSORT_API
# endif
#endif

/*- Datatypes -*/
#ifndef SAUCHAR_T
#define SAUCHAR_T
typedef uint8_t sauchar_t;
#endif /* SAUCHAR_T */
#ifndef SAINT_T
#define SAINT_T
typedef int32_t saint_t;
#endif /* SAINT_T */
#ifndef SAIDX64_T
#define SAIDX64_T
typedef int64_t saidx64_t;
#endif /* SAIDX64_T */
#ifndef PRIdSAINT_T
#define PRIdSAINT_T PRId32
#endif /* PRIdSAINT_T */
#ifndef PRIdSAIDX64_T
#define PRIdSAIDX64_T PRId64
#endif /* PRIdSAIDX64_T */


/*- Prototypes -*/

/**
 * Constructs the suffix array of a given string.
 * @param T[0..n-1] The input string.
 * @param SA[0..n-1] The output array of suffixes.
 * @param n The length of the given string.
 * @return 0 if no error occurred, -1 or -2 otherwise.
 */
DIVSUFSORT_API
saint_t
divsufsort64(const sauchar_t *T, saidx64_t *SA, saidx64_t n);

/**
 * Constructs the burrows-wheeler transformed string of a given string.
 * @param T[0..n-1] The input string.
 * @param U[0..n-1] The output string. (can be T)
 * @param A[0..n-1] The temporary array. (can be NULL)
 * @param n The length of the given string.
 * @return The primary index if no error occurred, -1 or -2 otherwise.
 */
DIVSUFSORT_API
saidx64_t
divbwt64(const sauchar_t *T, sauchar_t *U, saidx64_t *A, saidx64_t n);

/**
 * Returns the version of the divsufsort library.
 * @return The version number string.
 */
DIVSUFSORT_API
const char *
divsufsort64_version(void);


/**
 * Constructs the burrows-wheeler transformed string of a given string and suffix array.
 * @param T[0..n-1] The input string.
 * @param U[0..n-1] The output string. (can be T)
 * @param SA[0..n-1] The suffix array. (can be NULL)
 * @param n The length of the given string.
 * @param idx The output primary index.
 * @return 0 if no error occurred, -1 or -2 otherwise.
 */
DIVSUFSORT_API
saint_t
bw_transform64(const sauchar_t *T, sauchar_t *U,
             saidx64_t *SA /* can NULL */,
             saidx64_t n, saidx64_t *idx);

/**
 * Inverse BW-transforms a given BWTed string.
 * @param T[0..n-1] The input string.
 * @param U[0..n-1] The output string. (can be T)
 * @param A[0..n-1] The temporary array. (can be NULL)
 * @param n The length of the given string.
 * @param idx The primary index.
 * @return 0 if no error occurred, -1 or -2 otherwise.
 */
DIVSUFSORT_API
saint_t
inverse_bw_transform64(const sauchar_t *T, sauchar_t *U,
                     saidx64_t *A /* can NULL */,
                     saidx64_t n, saidx64_t idx);

/**
 * Checks the correctness of a given suffix array.
 * @param T[0..n-1] The input string.
 * @param SA[0..n-1] The input suffix array.
 * @param n The length of the given string.
 * @param verbose The verbose mode.
 * @return 0 if no error occurred.
 */
DIVSUFSORT_API
saint_t
sufcheck64(const sauchar_t *T, const saidx64_t *SA, saidx64_t n, saint_t verbose);

/**
 * Search for the pattern P in the string T.
 * @param T[0..Tsize-1] The input string.
 * @param Tsize The length of the given string.
 * @param P[0..Psize-1] The input pattern string.
 * @param Psize The length of the given pattern string.
 * @param SA[0..SAsize-1] The input suffix array.
 * @param SAsize The length of the given suffix array.
 * @param idx The output index.
 * @return The count of matches if no error occurred, -1 otherwise.
 */
DIVSUFSORT_API
saidx64_t
sa_search64(const sauchar_t *T, saidx64_t Tsize,
          const sauchar_t *P, saidx64_t Psize,
          const saidx64_t *SA, saidx64_t SAsize,
          saidx64_t *left);

/**
 * Search for the character c in the string T.
 * @param T[0..Tsize-1] The input string.
 * @param Tsize The length of the given string.
 * @param SA[0..SAsize-1] The input suffix array.
 * @param SAsize The length of the given suffix array.
 * @param c The input character.
 * @param idx The output index.
 * @return The count of matches if no error occurred, -1 otherwise.
 */
DIVSUFSORT_API
saidx64_t
sa_simplesearch64(const sauchar_t *T, saidx64_t Tsize,
                const saidx64_t *SA, saidx64_t SAsize,
                saint_t c, saidx64_t *left);


#ifdef __cplusplus
} /* extern "C" */
#endif /* __cplusplus */

#endif /* _DIVSUFSORT64_H */

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////

// #ifndef _DIVSUFSORT_PRIVATE_H
// #define _DIVSUFSORT_PRIVATE_H 1

#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */

// #if HAVE_CONFIG_H
// # include "config.h"
// #endif
#include <assert.h>
#include <stdio.h>
#if HAVE_STRING_H
# include <string.h>
#endif
#if HAVE_STDLIB_H
# include <stdlib.h>
#endif
#if HAVE_MEMORY_H
# include <memory.h>
#endif
#if HAVE_STDDEF_H
# include <stddef.h>
#endif
#if HAVE_STRINGS_H
# include <strings.h>
#endif
#if HAVE_INTTYPES_H
# include <inttypes.h>
#else
# if HAVE_STDINT_H
#  include <stdint.h>
# endif
#endif
// #if defined(BUILD_DIVSUFSORT64)
// # include "divsufsort64.h"
# ifndef SAIDX_T
#  define SAIDX_T
#  define saidx_t saidx64_t
# endif /* SAIDX_T */
# ifndef PRIdSAIDX_T
#  define PRIdSAIDX_T PRIdSAIDX64_T
# endif /* PRIdSAIDX_T */
// # define divsufsort divsufsort64
// # define divbwt divbwt64
// # define divsufsort_version divsufsort64_version
// # define bw_transform bw_transform64
// # define inverse_bw_transform inverse_bw_transform64
// # define sufcheck sufcheck64
// # define sa_search sa_search64
// # define sa_simplesearch sa_simplesearch64
// # define sssort sssort64
// # define trsort trsort64
// #else
// # include "divsufsort.h"
// #endif


/*- Constants -*/
#if !defined(UINT8_MAX)
# define UINT8_MAX (255)
#endif /* UINT8_MAX */
#if defined(ALPHABET_SIZE) && (ALPHABET_SIZE < 1)
# undef ALPHABET_SIZE
#endif
#if !defined(ALPHABET_SIZE)
# define ALPHABET_SIZE (UINT8_MAX + 1)
#endif
/* for divsufsort.c */
#define BUCKET_A_SIZE (ALPHABET_SIZE)
#define BUCKET_B_SIZE (ALPHABET_SIZE * ALPHABET_SIZE)
/* for sssort.c */
#if defined(SS_INSERTIONSORT_THRESHOLD)
# if SS_INSERTIONSORT_THRESHOLD < 1
#  undef SS_INSERTIONSORT_THRESHOLD
#  define SS_INSERTIONSORT_THRESHOLD (1)
# endif
#else
# define SS_INSERTIONSORT_THRESHOLD (8)
#endif
#if defined(SS_BLOCKSIZE)
# if SS_BLOCKSIZE < 0
#  undef SS_BLOCKSIZE
#  define SS_BLOCKSIZE (0)
# elif 32768 <= SS_BLOCKSIZE
#  undef SS_BLOCKSIZE
#  define SS_BLOCKSIZE (32767)
# endif
#else
# define SS_BLOCKSIZE (1024)
#endif
/* minstacksize = log(SS_BLOCKSIZE) / log(3) * 2 */
#if SS_BLOCKSIZE == 0
// # if defined(BUILD_DIVSUFSORT64)
#  define SS_MISORT_STACKSIZE (96)
// # else
// #  define SS_MISORT_STACKSIZE (64)
// # endif
#elif SS_BLOCKSIZE <= 4096
# define SS_MISORT_STACKSIZE (16)
#else
# define SS_MISORT_STACKSIZE (24)
#endif
// #if defined(BUILD_DIVSUFSORT64)
# define SS_SMERGE_STACKSIZE (64)
// #else
// # define SS_SMERGE_STACKSIZE (32)
// #endif
/* for trsort.c */
#define TR_INSERTIONSORT_THRESHOLD (8)
// #if defined(BUILD_DIVSUFSORT64)
# define TR_STACKSIZE (96)
// #else
// # define TR_STACKSIZE (64)
// #endif


/*- Macros -*/
#ifndef SWAP
# define SWAP(_a, _b) do { t = (_a); (_a) = (_b); (_b) = t; } while(0)
#endif /* SWAP */
#ifndef MIN
# define MIN(_a, _b) (((_a) < (_b)) ? (_a) : (_b))
#endif /* MIN */
#ifndef MAX
# define MAX(_a, _b) (((_a) > (_b)) ? (_a) : (_b))
#endif /* MAX */
#define STACK_PUSH(_a, _b, _c, _d)\
  do {\
    assert(ssize < STACK_SIZE);\
    stack[ssize].a = (_a), stack[ssize].b = (_b),\
    stack[ssize].c = (_c), stack[ssize++].d = (_d);\
  } while(0)
#define STACK_PUSH5(_a, _b, _c, _d, _e)\
  do {\
    assert(ssize < STACK_SIZE);\
    stack[ssize].a = (_a), stack[ssize].b = (_b),\
    stack[ssize].c = (_c), stack[ssize].d = (_d), stack[ssize++].e = (_e);\
  } while(0)
#define STACK_POP(_a, _b, _c, _d)\
  do {\
    assert(0 <= ssize);\
    if(ssize == 0) { return; }\
    (_a) = stack[--ssize].a, (_b) = stack[ssize].b,\
    (_c) = stack[ssize].c, (_d) = stack[ssize].d;\
  } while(0)
#define STACK_POP5(_a, _b, _c, _d, _e)\
  do {\
    assert(0 <= ssize);\
    if(ssize == 0) { return; }\
    (_a) = stack[--ssize].a, (_b) = stack[ssize].b,\
    (_c) = stack[ssize].c, (_d) = stack[ssize].d, (_e) = stack[ssize].e;\
  } while(0)
/* for divsufsort.c */
#define BUCKET_A(_c0) bucket_A[(_c0)]
#if ALPHABET_SIZE == 256
#define BUCKET_B(_c0, _c1) (bucket_B[((_c1) << 8) | (_c0)])
#define BUCKET_BSTAR(_c0, _c1) (bucket_B[((_c0) << 8) | (_c1)])
#else
#define BUCKET_B(_c0, _c1) (bucket_B[(_c1) * ALPHABET_SIZE + (_c0)])
#define BUCKET_BSTAR(_c0, _c1) (bucket_B[(_c0) * ALPHABET_SIZE + (_c1)])
#endif


/*- Private Prototypes -*/
/* sssort.c */
void
sssort(const sauchar_t *Td, const saidx_t *PA,
       saidx_t *first, saidx_t *last,
       saidx_t *buf, saidx_t bufsize,
       saidx_t depth, saidx_t n, saint_t lastsuffix);
/* trsort.c */
void
trsort(saidx_t *ISA, saidx_t *SA, saidx_t n, saidx_t depth);


#ifdef __cplusplus
} /* extern "C" */
#endif /* __cplusplus */

// #endif /* _DIVSUFSORT_PRIVATE_H */

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#include <sdsl/divsufsort_impl.hpp>
