#include "file.h"

#define FP(file) ((file)->open ? (file)->fp : NULL)

XFile*
xfile(const char* path, const char* mode)
{
  XFile* file = (XFile*)dmalloc(sizeof(XFile));
  file->path = path;
  file->mode = mode;
  file->length = 0;
  file->offset = 0;
  file->buff_offset = 0;
  file->buff_size = FILE_BUFF_SIZE;
  file->no_more_data = 0;
  file->line = NULL;
  file->open = 0;
  memset(file->buff, 0, sizeof(char) * FILE_BUFF_SIZE);
  return file;
}

void
destory_xfile(XFile* file)
{
  file->open = 0;
  if (file) {
    dfree(file, sizeof(XFile));
  }
}

int
plain_seek(XFile* file, long int offset, int whence)
{
  void* fp = FP(file);
  if (fp == NULL) {
    return 0;
  }
  fseek(fp, offset, whence);
  return 1;
}

static inline int
plain_reset(XFile* file)
{
  file->buff_offset = 0;
  file->no_more_data = 0;
  file->buff_size = FILE_BUFF_SIZE;
  file->offset = 0;
  return plain_seek(file, 0, SEEK_SET);
}

size_t
plain_length(XFile* file)
{
  plain_seek(file, 0, SEEK_END);
  size_t length = ftell(file->fp);
  plain_seek(file, file->offset, SEEK_SET);
  return length;
}

XFile*
plain_open(const char* path, const char* mode)
{
  FILE* fp;
  if ((fp = fopen(path, mode)) == NULL) {
    fprintf(stderr, "Open file: %s failed. %s\n", path, display_error);
    return NULL;
  }
  XFile* file = xfile(path, mode);
  file->fp = fp;
  file->open = 1;
  file->length = plain_length(file);
  file->offset = 0;
  file->buff_offset = 0;
  file->buff_size = FILE_BUFF_SIZE;
  return file;
}

size_t
plain_read(XFile* file, size_t size, char* buff)
{
  void* fp = FP(file);
  if (fp == NULL) {
    return 0;
  }
  return fread(buff, sizeof(char), size, fp);
}

size_t
plain_write(XFile* file, size_t size, char* buff)
{
  void* fp = FP(file);
  if (fp == NULL) {
    return 0;
  }
  return fwrite(buff, sizeof(char), size, fp);
}

size_t
plain_write_format(XFile* file, char* format, ...)
{
  void* fp = FP(file);
  if (fp == NULL) {
    return 0;
  }

  // TODO
  char buff[1024] = { 0 };
  va_list valist;
  va_start(valist, format);
  vsnprintf(buff, sizeof(buff), format, valist);
  va_end(valist);
  return plain_write(fp, strlen(buff), buff);
}

int
plain_close(XFile* file)
{
  void* fp = FP(file);
  if (fp == NULL) {
    return 0;
  }
  fclose(fp);
  destory_xfile(file);
  return 1;
}

typedef size_t (*read_t)(XFile* file, size_t size, char* buff);
static __always_inline line_t*
readline(XFile* file, line_t* line, read_t read_func)
{
  void* fp = FP(file);

  if (fp == NULL) {
    return NULL;
  }

  if (line == NULL) {
    line = dmalloc(sizeof(line_t));
    line->line = NULL;
    line->size = 0;
    line->capacity = 0;
  }

  char c;
  char* s = NULL;

  if (file->no_more_data && file->buff_offset >= file->buff_size) {
    line->size = 0;
    goto clean;
  }

  size_t s_len = 0;
  size_t s_cap = 0;
  if (line->line == NULL) {
    s = line->line = dmalloc(MIN_BUFF_SIZE * sizeof(char) + 1);
    line->capacity = MIN_BUFF_SIZE;
    memset(s, 0, MIN_BUFF_SIZE * sizeof(char) + 1);
  } else {
    s = line->line;
    if (s) {
      memset(s, 0, strlen(s));
    }
  }
  s_cap = line->capacity;
  while (1) {

    if (file->no_more_data && file->buff_offset >= file->buff_size) {
      if (s_len) {
        goto success;
      }
      goto clean;
    }

    if (file->buff_offset >= file->buff_size) {
      file->buff_offset = 0;
    }

    int i = file->buff_offset;

    size_t read_size = 0;
    if (!i && !file->no_more_data) {
      read_size = read_func(file, FILE_BUFF_SIZE, file->buff);
      file->offset += read_size;
      if (read_size < FILE_BUFF_SIZE) {
        file->no_more_data = 1;
      }
      file->buff_size = read_size;
    }

    for (; i < file->buff_size; i++) {
      c = file->buff[i];
      if (c == LF) {
        // LF
        file->buff_offset++;
        goto success;
      } else if (c == CR) {
        if (file->buff[i + 1] == LF) {
          // CR LF
          file->buff_offset += 2;
        } else {
          // CR
          file->buff_offset++;
        }
        goto success;
      }
      s[s_len++] = c;
      file->buff_offset++;
      if (s_len == s_cap - 1) {
        size_t old_cap = s_cap;
        s_cap += s_cap / 5;
        line->line = s =
            drealloc(s, old_cap * sizeof(char) + 1, s_cap * sizeof(char) + 1);
        line->capacity = s_cap;
      }
    }
  }

clean:
  if (line) {
    if (line->line) {
      dfree(line->line, line->capacity + 1);
    }
    dfree(line, sizeof(line_t));
  }
  file->line = NULL;
  return NULL;

success:
  line->size = s_len;
  line->line[s_len] = '\0';
  file->line = line;
  return line;
}

line_t*
plain_readline(XFile* file, line_t* line)
{
  return readline(file, line, plain_read);
}

size_t
plain_count(XFile* file, char c)
{
  size_t count = 0;
  while (!file->no_more_data) {
    size_t read_size = plain_read(file, FILE_BUFF_SIZE, file->buff);
    if (read_size < FILE_BUFF_SIZE) {
      file->no_more_data = 1;
    }
    for (int i = 0; i < read_size; i++) {
      if (file->buff[i] == c) {
        count++;
      }
    }
  }
  plain_reset(file);
  return count;
}

int
zip_close(XFile* file)
{
  gzFile fp = FP(file);
  if (fp == NULL) {
    return 0;
  }
  gzclose(fp);
  destory_xfile(file);
  return 1;
}

static inline size_t
zip_length(XFile* file)
{
  return 0;
}

static inline int
isgzip(const char* path, const char* mode)
{
  FILE* tmpfp = fopen(path, mode);
  if (tmpfp == NULL) {
    // catch error and display err msg
    fprintf(stderr, "Open file: %s failed. %s\n", path, display_error);
    return 0;
  }
  unsigned char gzip_flag[2] = { 0 };
  if (fread(gzip_flag, sizeof(char), 2, tmpfp) < 2) {
    fprintf(stderr, "Open file: %s failed.\n", path);
    fclose(tmpfp);
    return 0;
  }
  if (!(gzip_flag[0] = 0x1f && gzip_flag[1] == 0x8b)) {
    fprintf(stderr,
            "Open file: %s failed. It looks like a plain format file.\n",
            path);
    fclose(tmpfp);
    return 0;
  }
  fclose(tmpfp);
  return 1;
}

XFile*
zip_open(const char* path, const char* mode)
{
  if (!isgzip(path, mode)) {
    return NULL;
  }
  gzFile fp;
  if ((fp = gzopen(path, mode)) == NULL) {
    fprintf(stderr, "Open file: %s failed. %s\n", path, display_error);
    return NULL;
  }
  XFile* file = xfile(path, mode);
  file->fp = fp;
  file->open = 1;
  file->length = zip_length(file);
  file->offset = 0;
  file->buff_offset = 0;
  file->buff_size = FILE_BUFF_SIZE;
  return file;
}

int
zip_seek(XFile* file, long int offset, int whence)
{
  void* fp = FP(file);
  if (fp == NULL) {
    return 0;
  }
  gzseek(fp, offset, whence);
  return 1;
}

static inline int
zip_reset(XFile* file)
{
  file->buff_offset = 0;
  file->no_more_data = 0;
  file->buff_size = FILE_BUFF_SIZE;
  file->offset = 0;
  return zip_seek(file, 0, SEEK_SET);
}

size_t
zip_read(XFile* file, size_t size, char* buff)
{
  void* fp = FP(file);
  if (fp == NULL) {
    return 0;
  }
  return gzread(fp, buff, sizeof(char) * size);
}

size_t
zip_write(XFile* file, size_t size, char* buff)
{
  void* fp = FP(file);
  if (fp == NULL) {
    return 0;
  }
  return gzwrite(fp, buff, sizeof(char) * size);
}

size_t
zip_write_format(XFile* file, char* format, ...)
{
  void* fp = FP(file);
  if (fp == NULL) {
    return 0;
  }

  char buff[1024] = { 0 };
  va_list valist;
  va_start(valist, format);
  vsnprintf(buff, sizeof(buff), format, valist);
  va_end(valist);
  return zip_write(fp, strlen(buff), buff);
}

line_t*
zip_readline(XFile* file, line_t* line)
{
  return readline(file, line, zip_read);
}

size_t
zip_count(XFile* file, char c)
{
  size_t count = 0;
  while (!file->no_more_data) {
    size_t read_size = zip_read(file, FILE_BUFF_SIZE, file->buff);
    if (read_size < FILE_BUFF_SIZE) {
      file->no_more_data = 1;
    }
    for (int i = 0; i < read_size; i++) {
      if (file->buff[i] == c) {
        count++;
      }
    }
  }
  zip_reset(file);
  return count;
}

XFileHandle plainFile = {
  .close = plain_close,
  .open = plain_open,
  .read = plain_read,
  .readline = plain_readline,
  .write = plain_write,
  .write_format = plain_write_format,
  .seek = plain_seek,
  .length = plain_length,
  .count = plain_count,
  .reset = plain_reset,
};

XFileHandle zipFile = {
  .close = zip_close,
  .open = zip_open,
  .read = zip_read,
  .seek = zip_seek,
  .write = zip_write,
  .write_format = zip_write_format,
  .readline = zip_readline,
  .length = zip_length,
  .count = zip_count,
  .reset = zip_reset,
};
