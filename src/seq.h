#pragma once

#include "alloc.h"
#include "file.h"
#include "log.h"

#define SeqInitSize 80

typedef enum {
  FASTA,
  FASTQ,
} SeqType;

typedef struct {
  char* name;
  char* seq;
  char* qual;
  size_t len;
  size_t cap;
} Seq;

static inline Seq*
init_seq(int size, SeqType type)
{
  Seq* seq = (Seq*)dmalloc(sizeof(Seq));
  seq->name = NULL;
  seq->len = 0;
  if (size) {
    if (type == FASTQ) {
      seq->qual = (char*)dmalloc(size + 1);
      seq->seq = (char*)dmalloc(size + 1);
      seq->cap = size;
    } else {
      seq->qual = NULL;
      seq->seq = (char*)dmalloc(size + 1);
      seq->cap = size;
    }
  } else {
    seq->qual = NULL;
    seq->seq = NULL;
    seq->cap = 0;
  }
  return seq;
}

static inline void
free_seq(Seq* seq)
{
  if (seq->name) {
    dfree(seq->name, sizeof(char) * (strlen(seq->name) + 1));
  }
  if (seq->seq) {
    dfree(seq->seq, sizeof(char) * (seq->cap + 1));
  }
  if (seq->qual) {
    dfree(seq->qual, sizeof(char) * (seq->cap + 1));
  }
  dfree(seq, sizeof(Seq));
}

static inline Seq*
seq_shrink(Seq* seq)
{
  if (seq->len == seq->cap) {
    return seq;
  }
  seq->seq = (char*)drealloc(seq->seq, seq->cap + 1, seq->len + 1);
  seq->cap = seq->len;
  seq->seq[seq->len] = '\0';
  return seq;
}

static inline Seq*
__read_fasta(XFile* file, XFileHandle* handle, Seq* s)
{
  Seq* seq = s;
  line_t* line = NULL;
  if (file->offset && file->line == NULL) {
    if (seq) {
      handle->close(file);
      free_seq(seq);
    }
    return NULL;
  }
  if (seq == NULL) {
    seq = init_seq(SeqInitSize, FASTA);
  } else {
    seq->len = 0;
    dfree(seq->name, sizeof(char) * (strlen(seq->name) + 1));
    seq->name = NULL;
  }
  if (file->line && file->line->size && file->line->line[0] == '>') {
    line = file->line;
    seq->name = (char*)dmalloc(line->size);
    memcpy(seq->name, line->line + 1, line->size - 1);
    seq->name[line->size - 1] = '\0';
  }
  while ((line = handle->readline(file, line)) != NULL) {
    if (line->size == 0) {
      continue;
    }
    if (line->size && line->line && line->line[0] == '>') {
      if (seq->name != NULL) {
        return seq;
      }
      seq->name = (char*)dmalloc(line->size);
      memcpy(seq->name, line->line + 1, line->size - 1);
      seq->name[line->size - 1] = '\0';
      continue;
    }
    if (seq->len + line->size > seq->cap) {
      size_t old_cap = seq->cap;
      seq->cap = roundup(seq->len + line->size);
      if (seq->seq) {
        seq->seq = (char*)drealloc(seq->seq, old_cap + 1, seq->cap + 1);
      } else {
        seq->seq = (char*)dmalloc(seq->cap + 1);
      }
    }
    memcpy(seq->seq + seq->len, line->line, line->size);
    seq->len += line->size;
    seq->seq[seq->len] = '\0';
  }
  if (seq && !seq->len) {
    debug("finish read fasta file: %s\n", file->path);
    free_seq(seq);
    handle->close(file);
    return NULL;
  }
  return seq;
}

static inline Seq*
__read_fastq(XFile* file, XFileHandle* handle, Seq* s)
{
  Seq* seq = s;
  line_t* line = NULL;
  if (file->offset && file->line == NULL) {
    if (seq) {
      free_seq(seq);
    }
    return NULL;
  }
  if (seq == NULL) {
    seq = init_seq(0, FASTQ);
  } else {
    seq->len = 0;
    dfree(seq->name, sizeof(char) * (strlen(seq->name) + 1));
    seq->name = NULL;
  }
  int count = 0;
  for (; count < 4; count++) {
    line = handle->readline(file, line);
    if (line == NULL || line->size == 0) {
      goto clean;
    }
    if (count == 0) {
      seq->name = (char*)dmalloc(line->size);
      memcpy(seq->name, line->line + 1, line->size - 1);
      seq->name[line->size - 1] = '\0';
      continue;
    }

    if (count == 1) {
      if (line->size > seq->cap) {
        if (seq->seq) {
          seq->seq = (char*)drealloc(seq->seq, seq->cap + 1, line->size + 1);
          seq->qual = (char*)drealloc(seq->qual, seq->cap + 1, line->size + 1);
        } else {
          seq->seq = (char*)dmalloc(line->size + 1);
          seq->qual = (char*)dmalloc(line->size + 1);
        }
        seq->seq[0] = '\0';
        seq->qual[0] = '\0';
        seq->cap = line->size;
      }
      memcpy(seq->seq, line->line, line->size);
      seq->len = line->size;
      seq->seq[line->size] = '\0';
      continue;
    }
    if (count == 2) {
      // + check here
      if (line->size != 1) {
        fprintf(stderr, "fastq file <%s> format error\n", file->path);
        goto clean;
      }
      continue;
    }
    if (count == 3) {
      memcpy(seq->qual, line->line, line->size);
      seq->qual[line->size] = '\0';
      continue;
    }
  }
  if (!seq->len) {
    goto clean;
  }
  return seq;

clean:
  if (seq) {
    free_seq(seq);
    handle->close(file);
  }
  return NULL;
}

#define iter_fasta(path, seq)                                                 \
  XFileHandle handle = choice_handle(path);                                   \
  XFile* file = handle.open(path, "rb");                                      \
  if (file == NULL) {                                                         \
    fprintf(stderr, "Open file: %s failed.\n", path);                         \
    exit(1);                                                                  \
  }                                                                           \
  while ((seq = __read_fasta(file, &handle, seq)) != NULL)

#define iter_fastq(path, seq)                                                 \
  XFileHandle handle = choice_handle(path);                                   \
  XFile* file = handle.open(path, "rb");                                      \
  if (file == NULL) {                                                         \
    fprintf(stderr, "Open file: %s failed.\n", path);                         \
    return 1;                                                                 \
  }                                                                           \
  while ((seq = __read_fastq(file, &handle, seq)) != NULL)
