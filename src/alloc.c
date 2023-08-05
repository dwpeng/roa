#include "alloc.h"

unsigned long int useMemory = 0;

unsigned long int
getUsedMemory()
{
  return useMemory;
}
