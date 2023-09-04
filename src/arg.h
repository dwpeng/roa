#pragma once

#define arginit(name) void name(int argc, char** argv)

#define argstart()                                                            \
  int offset = 0;                                                             \
  while (1)

#define argend()                                                              \
  if (offset >= argc) {                                                       \
    break;                                                                    \
  } else {                                                                    \
    printf("Got unknown argument: %s\n", argv[offset]);                       \
    exit(1);                                                                  \
  }

#define argbreak()                                                            \
  if (offset >= argc) {                                                       \
    break;                                                                    \
  }

#define arglost_value(argv, offset, argname, count)                           \
  if (offset + count >= argc) {                                               \
    printf("%s requires %d arguments\n", argname, count);                     \
    exit(1);                                                                  \
  }                                                                           \
  if (offset + count >= argc) {                                               \
    printf("%s requires %d arguments\n", argname, count);                     \
    exit(1);                                                                  \
  }                                                                           \
  if (strcmp(argv[offset + count], "-") == 0) {                               \
    printf("%s requires %d arguments\n", argname, count);                     \
    exit(1);                                                                  \
  }

#define argbool(argname, argvalue)                                            \
  if (strcmp(argv[offset], argname) == 0) {                                   \
    arglost_value(argv, offset, argname, 1);                                  \
    argvalue = !!(atoi(argv[offset + 1]));                                    \
    offset += 2;                                                              \
    argbreak();                                                               \
    continue;                                                                 \
  }

#define argstring(argname, argvalue)                                          \
  if (strcmp(argv[offset], argname) == 0) {                                   \
    arglost_value(argv, offset, argname, 1);                                  \
    argvalue = argv[offset + 1];                                              \
    offset += 2;                                                              \
    argbreak();                                                               \
    continue;                                                                 \
  }

#define argint(argname, argvalue)                                             \
  if (strcmp(argv[offset], argname) == 0) {                                   \
    arglost_value(argv, offset, argname, 1);                                  \
    argvalue = atoi(argv[offset + 1]);                                        \
    offset += 2;                                                              \
    argbreak();                                                               \
    continue;                                                                 \
  }

#define argfloat(argname, argvalue)                                           \
  if (strcmp(argv[offset], argname) == 0) {                                   \
    argvalue = atof(argv[offset + 1]);                                        \
    arglost_value(argv, offset, argname, 1);                                  \
    offset += 2;                                                              \
    argbreak();                                                               \
    continue;                                                                 \
  }

#define argpass(argname)                                                      \
  if (strcmp(argv[offset], argname) == 0) {                                   \
    offset++;                                                                 \
    argbreak();                                                               \
    continue;                                                                 \
  }
