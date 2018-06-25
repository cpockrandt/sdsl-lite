#ifndef _CONFIG_H
#define _CONFIG_H 1

#define PROJECT_VERSION_FULL "2.0.1-15-g22e6b23"

#include <assert.h>
#include <stdio.h>
#include <string.h>

// optional
#include <stdlib.h>
#include <memory.h>
#include <stddef.h>
#include <strings.h>

#include <inttypes.h> // alternative: #include <stdint.h>

#if !defined(UINT8_MAX)
# define UINT8_MAX (255)
#endif

#if defined(ALPHABET_SIZE) && (ALPHABET_SIZE < 1)
# undef ALPHABET_SIZE
#endif
#if !defined(ALPHABET_SIZE)
# define ALPHABET_SIZE (UINT8_MAX + 1)
#endif

#define BUCKET_A_SIZE (ALPHABET_SIZE)
#define BUCKET_B_SIZE (ALPHABET_SIZE * ALPHABET_SIZE)

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

#endif /* _CONFIG_H */
