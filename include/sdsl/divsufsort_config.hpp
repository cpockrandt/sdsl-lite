#ifndef _CONFIG_H
#define _CONFIG_H 1

/** Define to the version of this package. **/
#define PROJECT_VERSION_FULL "2.0.1-15-g22e6b23"

/** Define to 1 if you have the header files. **/
#define HAVE_INTTYPES_H 1
#define HAVE_STDDEF_H 1
#define HAVE_STDINT_H 1
#define HAVE_STDLIB_H 1
#define HAVE_STRING_H 1
#define HAVE_STRINGS_H 1
#define HAVE_MEMORY_H 1
#define HAVE_SYS_TYPES_H 1

#ifndef HAVE__SETMODE
# if HAVE_SETMODE
#  define _setmode setmode
#  define HAVE__SETMODE 1
# endif
# if HAVE__SETMODE && !HAVE__O_BINARY
#  define _O_BINARY 0
#  define HAVE__O_BINARY 1
# endif
#endif

#endif /* _CONFIG_H */
