#include "log.h"
#include <stdio.h>

int __pglog_default_level = PGLOG_LEVEL_INFO;

void
log_set_level(int level)
{
  __pglog_default_level = level;
}

int
log_get_level()
{
  return __pglog_default_level;
}
