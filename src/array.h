#pragma once

#include <stdio.h>
#include <stdlib.h>

#include "alloc.h"
#include "log.h"

typedef struct {
  void** data;
  size_t size;
  size_t capacity;
} Array;

static inline Array*
arrayNew(size_t capacity)
{
  Array* array = (Array*)dmalloc(sizeof(Array));
  array->data = (void**)dmalloc(sizeof(void*) * capacity);
  array->size = 0;
  array->capacity = capacity;
  return array;
}

#define arrayFree(__array)                                                    \
  do {                                                                        \
    if ((__array) != NULL) {                                                  \
      if ((__array)->data) {                                                  \
        dfree((__array)->data, sizeof(void*) * (__array)->capacity);          \
      }                                                                       \
      dfree((__array), sizeof(Array));                                        \
    }                                                                         \
  } while (0)

#define arrayPush(__array__, __data)                                          \
  do {                                                                        \
    size_t old_size = (__array__)->capacity;                                  \
    if ((__array__)->size == (__array__)->capacity) {                         \
      (__array__)->capacity = roundup((__array__)->capacity + 1);             \
      (__array__)->data =                                                     \
          (void**)drealloc((__array__)->data, old_size * sizeof(void*),       \
                           sizeof(void*) * (__array__)->capacity);            \
    }                                                                         \
    (__array__)->data[(__array__)->size] = (__data);                          \
    (__array__)->size++;                                                      \
  } while (0)

static inline void*
arrayGet(Array* array, size_t index)
{
  if (index >= array->size) {
    return NULL;
  }
  return array->data[index];
}

static inline void
arraySet(Array* array, size_t index, void* data)
{
  if (index >= array->size) {
    return;
  }
  array->data[index] = data;
}

static inline void
arrayClear(Array* array)
{
  array->size = 0;
}

#define arrayShrink(array)                                                    \
  do {                                                                        \
    if ((array) != NULL) {                                                    \
      if ((array)->size < (array)->capacity) {                                \
        if ((array)->size == 0 && (array)->capacity > 0) {                    \
          dfree((array)->data, sizeof(void*) * (array)->capacity);            \
          (array)->data = NULL;                                               \
          (array)->capacity = 0;                                              \
          (array)->size = 0;                                                  \
        } else {                                                              \
          size_t old_size = (array)->capacity;                                \
          (array)->capacity = (array)->size;                                  \
          (array)->data =                                                     \
              (void**)drealloc((array)->data, old_size * sizeof(void*),       \
                               sizeof(void*) * (array)->capacity);            \
        }                                                                     \
      }                                                                       \
    }                                                                         \
  } while (0)

#define arrayExtend(arr1, arr2)                                               \
  do {                                                                        \
    for (size_t i = 0; i < (arr2)->size; i++) {                               \
      arrayPush((arr1), arrayGet((arr2), i));                                 \
    }                                                                         \
  } while (0)

#define arrayPop(__array__)                                                   \
  do {                                                                        \
    if ((__array__)->size > 0) {                                              \
      (__array__)->size--;                                                    \
    }                                                                         \
  } while (0)

#define arrayLast(__array__)                                                  \
  ((__array__)->size > 0 ? (__array__)->data[(__array__)->size - 1] : NULL)

