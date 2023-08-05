#pragma once

#include <errno.h>
#include <stdarg.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <zconf.h>
#include <zlib.h>

#include "alloc.h"
#include "file.h"

#define display_error strerror(errno)

#define MIN_BUFF_SIZE 128               // 128KB
#define FILE_BUFF_SIZE 1024 * 1024 * 64 // 64MB

#define LF '\n'
#define CR '\r'
#define CRLF "\r\n"

typedef struct {
  char* line;
  size_t size;
  size_t capacity;
} line_t;

typedef struct XFile XFile;
struct XFile {
  const char* path;
  const char* mode;
  int open;
  void* fp;
  // 文件读写指针
  long int offset;
  // 文件内容长度
  size_t length;
  // TODO
  char buff[FILE_BUFF_SIZE];
  // buff指针偏移量
  int buff_offset;
  // buff大小
  size_t buff_size;
  int no_more_data; // for zip file
  // store every line data
  line_t* line;
};

typedef struct XFileHandle {
  XFile* (*open)(const char* path, const char* mode);
  size_t (*read)(XFile* file, size_t size, char* buff);
  size_t (*write)(XFile* file, size_t size, char* buff);
  size_t (*write_format)(XFile* file, char* format, ...);
  int (*seek)(XFile* file, long int offset, int whence);
  int (*close)(XFile* file);
  size_t (*length)(XFile* file);
  line_t* (*readline)(XFile* file, line_t* line);
  size_t (*count)(XFile* file, char c);
  int (*reset)(XFile* file);
} XFileHandle;

extern XFileHandle plainFile;
extern XFileHandle zipFile;

static inline XFileHandle
choice_handle(const char* path)
{
  FILE* fp = fopen(path, "rb");
  if (fp == NULL) {
    printf("open file error: %s\n", display_error);
    exit(1);
  }
  unsigned char buff[2];
  size_t size = fread(buff, sizeof(char), 2, fp);
  if (size < 2) {
    return plainFile;
  }
  fclose(fp);
  if (buff[0] == 0x1f && buff[1] == 0x8b) {
    return zipFile;
  }
  return plainFile;
}

#define readline_iter(handle, file, line)                                     \
  while ((line = handle.readline(file, line)) != NULL)
