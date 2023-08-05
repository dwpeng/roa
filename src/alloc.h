#pragma once

#include "log.h"

#include <stdint.h>
#include <stdlib.h>

#ifdef __SIZE_TYPE__
typedef __SIZE_TYPE__ size_t;
#else
typedef unsigned long size_t;
#endif

extern unsigned long int useMemory;
extern unsigned long int getUsedMemory();

static inline unsigned long int
getUsedMemoryKB()
{
  return getUsedMemory() / 1024;
}

static inline unsigned long int
getUsedMemoryMB()
{
  return getUsedMemory() / 1024 / 1024;
}

static inline unsigned long int
getUsedMemoryGB()
{
  return getUsedMemory() / 1024 / 1024 / 1024;
}

// void *dmalloc(size_t size);
#define dmalloc(size)                                                         \
  ({                                                                          \
    void* ptr = malloc(size);                                                 \
    if (ptr == NULL) {                                                        \
      fprintf(stderr, "malloc failed. %s:%d\n", __FILE__, __LINE__);          \
      exit(1);                                                                \
    }                                                                         \
    useMemory += (size);                                                      \
    ptr;                                                                      \
  })

// void *drealloc(void *ptr, size_t old_size, size_t size);
#define drealloc(ptr, old_size, size)                                         \
  ({                                                                          \
    void* new_ptr = realloc(ptr, size);                                       \
    if (new_ptr == NULL) {                                                    \
      fprintf(stderr, "realloc failed. %s:%d\n", __FILE__, __LINE__);         \
      exit(1);                                                                \
    }                                                                         \
    useMemory += (size) - (old_size);                                         \
    new_ptr;                                                                  \
  })

// void dfree(void *ptr, size_t size);

#define dfree(ptr, size)                                                      \
  if ((ptr) != NULL) {                                                        \
    useMemory -= (size);                                                      \
    free((ptr));                                                              \
  }

// void *dcalloc(size_t nmemb, size_t size);
#define dcalloc(nmemb, size)                                                  \
  ({                                                                          \
    void* ptr = calloc(nmemb, size);                                          \
    if (ptr == NULL) {                                                        \
      fprintf(stderr, "calloc failed. %s:%d\n", __FILE__, __LINE__);          \
      exit(1);                                                                \
    }                                                                         \
    useMemory += (size);                                                      \
    ptr;                                                                      \
  })

#define roundup(x)                                                            \
  ({                                                                          \
    uint64_t __x = (x);                                                       \
    __x--;                                                                    \
    __x |= __x >> 1;  /* handle  2 bit numbers */                             \
    __x |= __x >> 2;  /* handle  4 bit numbers */                             \
    __x |= __x >> 4;  /* handle  8 bit numbers */                             \
    __x |= __x >> 8;  /* handle 16 bit numbers */                             \
    __x |= __x >> 16; /* handle 32 bit numbers */                             \
    __x |= __x >> 32; /* handle 64 bit numbers */                             \
    __x++;                                                                    \
    __x;                                                                      \
  })
