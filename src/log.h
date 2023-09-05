#pragma once

#include <stdio.h>
#include <string.h>
#include <time.h>
#include <unistd.h>

#define PGLOG_LEVEL_DEBUG 10
#define PGLOG_LEVEL_INFO 20
#define PGLOG_LEVEL_WARN 30
#define PGLOG_LEVEL_ERROR 40

static const char* __PLAIN_LEVELS[] = { "DEBUG", "INFO ", "WARN ", "ERROR" };
static const char* __COLOR_LEVELS[] = { "\033[32mDEBUG\033[0m",
                                        "\033[34mINFO \033[0m",
                                        "\033[33mWARN \033[0m",
                                        "\033[31mERROR\033[0m" };

#define __levelStr(level)                                                     \
  (level == PGLOG_LEVEL_DEBUG   ? __PLAIN_LEVELS[0]                           \
   : level == PGLOG_LEVEL_INFO  ? __PLAIN_LEVELS[1]                           \
   : level == PGLOG_LEVEL_WARN  ? __PLAIN_LEVELS[2]                           \
   : level == PGLOG_LEVEL_ERROR ? __PLAIN_LEVELS[3]                           \
                                : "UNKNOWN")

#define __levelColorStr(level)                                                \
  (level == PGLOG_LEVEL_DEBUG   ? __COLOR_LEVELS[0]                           \
   : level == PGLOG_LEVEL_INFO  ? __COLOR_LEVELS[1]                           \
   : level == PGLOG_LEVEL_WARN  ? __COLOR_LEVELS[2]                           \
   : level == PGLOG_LEVEL_ERROR ? __COLOR_LEVELS[3]                           \
                                : "UNKNOWN")

#if defined (__linux__) 
#define ISATTY(fp) isatty(fp->_fileno)
#else
#define ISAATTY(fp) 1
#endif


#define _STR(x) #x
#define STR(x) _STR(x)

#ifdef DEBUG
#define __header(file, line) "[" file ":" line "]"
#define __headerColor(file, line) "[\033[36m" file ":" line "\033[0m]"
#else
#define __header(file, line) ""
#define __headerColor(file, line) ""
#endif

extern int __pglog_default_level;
extern void log_set_level(int level);
extern int log_get_level();

static inline FILE*
__pglog_print_file(int level, FILE* stream, const char* header)
{
  if (level < __pglog_default_level)
    return NULL;
  char buff[30];
  time_t now = time(NULL);
  struct tm* t = localtime(&now);
  strftime(buff, sizeof(buff), "%Y-%m-%d %H:%M:%S", t);
  fprintf(stream, "[%s - %s]%s:", __levelStr(level), buff, header);
  return stream;
}

static inline FILE*
__pglog_print(int level, FILE* stream, const char* header)
{
  if (level < __pglog_default_level)
    return NULL;
  char buff[30];
  time_t now = time(NULL);
  struct tm* t = localtime(&now);
  strftime(buff, sizeof(buff), "%Y-%m-%d %H:%M:%S", t);
  fprintf(stream, "[%s - \033[36m%s\033[0m]%s:", __levelColorStr(level), buff,
          header);
  return stream;
}

#define __log2file(level, fp, ...)                                            \
  do {                                                                        \
    FILE* __fp =                                                              \
        __pglog_print_file(level, fp, __header(__FILE__, STR(__LINE__)));     \
    if (__fp) {                                                               \
      fprintf(__fp, __VA_ARGS__);                                             \
      fputc('\n', __fp);                                                      \
    }                                                                         \
  } while (0)

#define __log2terminal(level, ...)                                            \
  do {                                                                        \
    FILE* __tfp = level == PGLOG_LEVEL_ERROR ? stderr : stdout;               \
    if (!ISATTY(__tfp) || !ISATTY(stdout) || !ISATTY(stderr)) {               \
      __tfp = stdout;                                                         \
      FILE* __fp = __pglog_print_file(level, __tfp,                           \
                                      __header(__FILE__, STR(__LINE__)));     \
      if (__fp) {                                                             \
        fprintf(__fp, __VA_ARGS__);                                           \
        fputc('\n', __fp);                                                    \
      }                                                                       \
    } else {                                                                  \
      FILE* __fp = __pglog_print(level, __tfp,                                \
                                 __headerColor(__FILE__, STR(__LINE__)));     \
      if (__fp) {                                                             \
        fprintf(__fp, __VA_ARGS__);                                           \
        fputc('\n', __fp);                                                    \
      }                                                                       \
    }                                                                         \
  } while (0)

#define debug(...) __log2terminal(PGLOG_LEVEL_DEBUG, __VA_ARGS__)
#define info(...) __log2terminal(PGLOG_LEVEL_INFO, __VA_ARGS__)
#define warn(...) __log2terminal(PGLOG_LEVEL_WARN, __VA_ARGS__)
#define error(...) __log2terminal(PGLOG_LEVEL_ERROR, __VA_ARGS__)

#define fdebug(stream, ...)                                             \
  __log2file(PGLOG_LEVEL_DEBUG, stream, __VA_ARGS__)
#define finfo(stream, ...)                                              \
  __log2file(PGLOG_LEVEL_INFO, stream, __VA_ARGS__)
#define fwarn(stream, ...)                                              \
  __log2file(PGLOG_LEVEL_WARN, stream, __VA_ARGS__)
#define ferror(stream, ...)                                             \
  __log2file(PGLOG_LEVEL_ERROR, stream, __VA_ARGS__)
