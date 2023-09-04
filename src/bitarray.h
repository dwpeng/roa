#pragma once
#include "alloc.h"
#include <stdio.h>

// clang-format off
typedef struct {
  size_t size;
  size_t __realCols;
  char mask;
  int nbit;
  unsigned char* data;
} BitArray;

const unsigned char bitMask[8] = {
  0b00000001, 0b00000011, 0b00000111, 0b00001111,
  0b00011111, 0b00111111, 0b01111111, 0b11111111
};

// clang-format on

static inline BitArray*
bitarrayNew(size_t size, int nbit)
{
  if (nbit > 8 || nbit < 0) {
    return NULL;
  }

  if (nbit != 1) {
    if (nbit % 2) {
      printf("BitArray only supports nbit = 2^k\n");
      return NULL;
    }
  }

  BitArray* array = (BitArray*)dmalloc(sizeof(BitArray));
  array->size = size;
  array->mask = bitMask[nbit - 1];
  array->nbit = nbit;
  array->__realCols = (size * nbit + 7) / 8;
  array->data = (unsigned char*)dmalloc(array->__realCols);
  memset(array->data, 0, array->__realCols);
  return array;
}

static inline BitArray*
bitarrayClone(BitArray* array)
{
  if (array == NULL) {
    return NULL;
  }
  BitArray* clone = (BitArray*)dmalloc(sizeof(BitArray));
  clone->size = array->size;
  clone->__realCols = array->__realCols;
  clone->mask = array->mask;
  clone->nbit = array->nbit;
  clone->data = (unsigned char*)dmalloc(clone->__realCols);
  memcpy(clone->data, array->data, clone->__realCols);
  return clone;
}

static inline void
bitarrayFree(BitArray* array)
{
  if (array == NULL) {
    return;
  }
  dfree(array->data, sizeof(char) * array->__realCols);
  dfree(array, sizeof(BitArray));
}

static inline void
bitarraySet(BitArray* arr, size_t i, unsigned char value)
{
  if (arr == NULL) {
    return;
  }
  if (i >= arr->size) {
    printf("[Set] index out of bounds. size: %zu, i: %zu\n", arr->size, i);
    return;
  }
  arr->data[(i * arr->nbit) / 8] &= ~(arr->mask << (i * arr->nbit) % 8);
  arr->data[(i * arr->nbit) / 8] |=
      ((value & arr->mask) << (i * arr->nbit) % 8);
}

static inline unsigned char
bitarrayGet(BitArray* arr, size_t i)
{
  if (arr == NULL) {
    return 0;
  }
  if (i >= arr->size) {
    printf("[Get] index out of bounds. size: %zu, i: %zu\n", arr->size, i);
    return 0;
  }
  return (arr->data[(i * arr->nbit) / 8] >> ((i * arr->nbit) % 8)) & arr->mask;
}

static inline void
bitarrayOr(BitArray* arr, BitArray* other)
{
  if (arr == NULL || other == NULL) {
    return;
  }
  if (arr->size != other->size) {
    printf("[Or] size mismatch. arr->size: %zu, other->size: %zu\n", arr->size,
           other->size);
    return;
  }
#ifdef parallel
#pragma omp parallel for
#endif
  for (size_t i = 0; i < arr->__realCols; i++) {
    arr->data[i] |= other->data[i];
  }
}
