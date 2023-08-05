#include "alloc.h"
#include "array.h"
#include "bitarray.h"
#include "file.h"
#include "log.h"
#include "seq.h"

#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <zlib.h>

// I follow the rule of primer design from Qiagen
// https://www.qiagen.com/zh-us/knowledge-and-support/knowledge-hub/bench-guide/pcr/introduction/pcr-primer-design

// clang-format off
static inline void create_base2int(unsigned char* map){
  map['A'] = 0;
  map['C'] = 1;
  map['G'] = 2;
  map['T'] = 3;
  map['a'] = 0;
  map['c'] = 1;
  map['g'] = 2;
  map['t'] = 3;
  map['N'] = 4;
  map['n'] = 4;
}

static unsigned char int2base[5] = {
  'A', 'C', 'G', 'T', 'N'
};

static inline void int2KmerString(uint64_t kmer, int k, char *buff) {
  for (int i = 0; i < k; i++) {
    buff[k - 1 - i] = int2base[kmer & 0x3];
    kmer >>= 2;
  }
  buff[k] = '\0';
}

// clang-format on

#define SEQ_MAX_LEN (1UL << 29) - 1
#define KMER_MASK (1UL << 32) - 1
#define KMER_LEN 16
#define KMER_LONG_LEN 20
#define KMER_PER_CIRCLE 4

static inline int
isFileExist(const char* path)
{
  FILE* fp = fopen(path, "rb");
  if (fp == NULL) {
    return 0;
  }
  fclose(fp);
  return 1;
}

typedef struct {
  // 16-bit kmer
  int pos : 30;
  int drop : 1;
  int strand : 1;
  uint32_t kmer;
} Kmer;

typedef struct {
  Array* kmers; // Array<Array<Kmer>>
  Array* seqs;  // Array<Seq*>
} Query;

typedef struct {
  char* name;      // name of the segment
  size_t start;    // start position of the segment
  size_t end;      // end position of the segment
  BitArray* bases; // bases of the segment
  int vaild;       // if the segment is vaild
  float Tm;        // pcr melting temperature
} Segment;

typedef struct {
  const char* path;
  BitArray* index;
} Index;

static inline Index*
createIndex(Index* index, const char* path)
{
  if (index == NULL) {
    index = dmalloc(sizeof(Index));
    index->path = path;
    index->index = bitarrayNew(KMER_MASK + 1, 1);
  }
  unsigned char basemap[128] = { 4 };
  create_base2int(basemap);
  Seq* seq = NULL;
  {
    iter_fasta(path, seq)
    {
      uint32_t kmer = 0;
      uint32_t reverse_kmer = 0;
      char c = 4;
      size_t count = 0;
      for (size_t i = 0; i < seq->len - KMER_LEN + 1; i++) {
        c = basemap[seq->seq[i]];
        if (c == 4) {
          kmer = 0;
          reverse_kmer = 0;
          count = 0;
          continue;
        }
        kmer = (kmer << 2) | (c & 0x3) & KMER_MASK;
        reverse_kmer = (reverse_kmer >> 2) | ((0x00000003 - c) << 30);
        count++;
        if (count < KMER_LEN) {
          continue;
        }
        bitarraySet(index->index, reverse_kmer, 1);
        bitarraySet(index->index, kmer, 1);
      }
    }
  }
  return index;
}

static inline void
dumpIndex(Index* index, const char* path)
{
  gzFile fp = gzopen(path, "wb");
  gzwrite(fp, &index->index->size, sizeof(size_t));
  gzwrite(fp, &index->index->mask, sizeof(char));
  gzwrite(fp, &index->index->nbit, sizeof(int));
  gzwrite(fp, &index->index->__realCols, sizeof(size_t));
  gzwrite(fp, index->index->data, sizeof(uint8_t) * index->index->__realCols);
  gzclose(fp);
}

static inline Index*
loadIndex(const char* path)
{
  Index* index = dmalloc(sizeof(Index));
  BitArray* b = bitarrayNew(1, 1);
  dfree(b->data, sizeof(uint8_t) * b->__realCols);
  index->index = b;
  gzFile fp = gzopen(path, "rb");
  gzread(fp, &index->index->size, sizeof(size_t));
  gzread(fp, &index->index->mask, sizeof(char));
  gzread(fp, &index->index->nbit, sizeof(int));
  gzread(fp, &index->index->__realCols, sizeof(size_t));
  index->index->data = dmalloc(sizeof(uint8_t) * index->index->__realCols);
  gzread(fp, index->index->data, sizeof(uint8_t) * index->index->__realCols);
  gzclose(fp);
  return index;
}

static inline void
freeIndex(Index* index)
{
  bitarrayFree(index->index);
  dfree(index, sizeof(Index));
}

#ifdef parallel
#include <omp.h>
#endif

static inline Query*
createQuery(const char* path)
{
  uint32_t kmer = 0x00000000;
  unsigned char basemap[128] = { 4 };
  create_base2int(basemap);
  char c = 4;
  Query* query = dmalloc(sizeof(Query));
  query->kmers = arrayNew(10);
  query->seqs = arrayNew(10);
  Seq* seq = NULL;
  {
    iter_fasta(path, seq)
    {
      // create 16 mer
      int count = 0;

      // backup seq
      Seq* copy_seq = dmalloc(sizeof(Seq));
      copy_seq->name = dmalloc(sizeof(char) * (strlen(seq->name) + 1));
      strcpy(copy_seq->name, seq->name);
      copy_seq->name[strlen(seq->name)] = '\0';
      copy_seq->len = seq->len;
      copy_seq->cap = seq->len;
      copy_seq->seq = dmalloc(sizeof(char) * (seq->len + 1));
      copy_seq->qual = NULL;
      strcpy(copy_seq->seq, seq->seq);
      copy_seq->seq[seq->len] = '\0';
      arrayPush(query->seqs, copy_seq);

      Array* kmers = arrayNew(10);
      if (seq->len < KMER_LEN) {
        continue;
      }
      for (size_t i = 0; i < seq->len - KMER_LEN + 1; i++) {
        c = basemap[seq->seq[i]];
        if (c == 4) {
          kmer = 0x00000000;
          count = 0;
          continue;
        }
        kmer = (kmer << 2) | (c & 0x3) & KMER_MASK;
        count++;
        if (count < KMER_LEN) {
          continue;
        }
        Kmer* kmer_t = dmalloc(sizeof(Kmer));
        kmer_t->kmer = kmer;
        kmer_t->pos = i - KMER_LEN + 1;
        kmer_t->drop = 0;
        kmer_t->strand = 0;
        arrayPush(kmers, kmer_t);
      }
      arrayPush(query->kmers, kmers);
    }
  }
  return query;
}

static inline void
freeQuery(Query* query)
{
  for (size_t i = 0; i < query->kmers->size; i++) {
    Array* kmers = arrayGet(query->kmers, i);
    for (size_t j = 0; j < kmers->size; j++) {
      dfree(kmers->data[j], sizeof(Kmer));
    }
    arrayFree(kmers);
  }
  arrayFree(query->kmers);
  for (size_t i = 0; i < query->seqs->size; i++) {
    free_seq(query->seqs->data[i]);
  }
  arrayFree(query->seqs);
  dfree(query, sizeof(Query));
}

static inline void
vaildKmers(Query* query, Index* index)
{
  for (size_t i = 0; i < query->kmers->size; i++) {
    Array* kmers = query->kmers->data[i];
#ifdef parallel
    size_t t = 0;
    size_t start = 0;
    size_t end = kmers->size;
#pragma omp parallel private(start, end, t)
    t = omp_get_thread_num();
    start = (kmers->size / omp_get_num_threads()) * t;
    end = (kmers->size / omp_get_num_threads()) * (t + 1);
    if (omp_get_thread_num() == omp_get_num_threads() - 1) {
      end = kmers->size;
    }
#else
    size_t start = 0;
    size_t end = kmers->size;
#endif
    for (size_t j = start; j < end; j++) {
      Kmer* kmer = kmers->data[j];
      if (kmer->drop == 1) {
        continue;
      }
      if (bitarrayGet(index->index, kmer->kmer)) {
        kmer->drop = 1;
        continue;
      }
    }
#ifdef parallel
#pragma omp barrier
#endif
    // set not continuous kmer to drop
    size_t window_start = 0;
    size_t window_end = 0;
    int prev_status = -1;
    for (size_t i = 0; i < kmers->size; i++) {
      Kmer* kmer = kmers->data[i];
      if (prev_status == -1) {
        if (kmer->drop == 1) {
          continue;
        }
        prev_status = kmer->drop;
        window_start = i;
        continue;
      }
      if (prev_status != kmer->drop) {
        window_end = i;
        if (window_end - window_start < 4) {
          for (size_t j = window_start; j < window_end; j++) {
            ((Kmer*)kmers->data[j])->drop = 1;
          }
        }
        window_start = i;
        prev_status = kmer->drop;
        continue;
      }
    }
    if (prev_status != -1) {
      window_end = kmers->size;
      if (window_end - window_start < 4) {
        for (size_t j = window_start; j < window_end; j++) {
          ((Kmer*)kmers->data[j])->drop = 1;
        }
      }
    }
  }
}

#define collectSegmentMacro(segments, start, end)                             \
  do {                                                                        \
    Segment* segment = dmalloc(sizeof(Segment));                              \
    segment->start = start;                                                   \
    segment->end = end;                                                       \
    segment->bases = bitarrayNew(end - start + KMER_LEN + 1, 2);              \
    segment->name = ((Seq*)query->seqs->data[i])->name;                       \
    segment->vaild = 1;                                                       \
    {                                                                         \
      Kmer* kmer = kmers->data[start];                                        \
      uint32_t kmerint = kmer->kmer;                                          \
      for (size_t k = 0; k < KMER_LEN; k++) {                                 \
        bitarraySet(segment->bases, KMER_LEN - k - 1, kmerint & 0x3);         \
        kmerint >>= 2;                                                        \
      }                                                                       \
      for (size_t k = start + 1; k <= end; k++) {                             \
        kmer = kmers->data[k];                                                \
        bitarraySet(segment->bases, KMER_LEN + k - 1, kmer->kmer & 0x3);      \
      }                                                                       \
    }                                                                         \
    arrayPush(segments, segment);                                             \
  } while (0)

static inline Array*
collectSegment(Query* query)
{
  Array* segments = arrayNew(10);
  for (size_t i = 0; i < query->kmers->size; i++) {
    Array* kmers = query->kmers->data[i];
    long int start = -1;
    long int end = -1;
    for (size_t j = 0; j < kmers->size; j++) {
      Kmer* kmer = kmers->data[j];
      // collect all continuous kmers that not drop
      // and represent a segment
      if (kmer->drop == 1) {
        if (start != -1) {
          if (end - start > 5) {
            collectSegmentMacro(segments, start, end);
          }
          start = -1;
          end = -1;
        }
        continue;
      }
      if (start == -1) {
        start = kmer->pos;
      }
      end = kmer->pos;
    }
    // collect the last segment
    if (start != -1 && end != -1 && ((Kmer*)kmers->data[start])->drop == 0) {
      if (end - start > 5) {
        collectSegmentMacro(segments, start, end);
      }
    }
  }
  debug("collect %zu segments", segments->size);
  return segments;
}

static float GC_RATE[21] = {
  0.0,         1.0 / 20.0,  2.0 / 20.0,  3.0 / 20.0,  4.0 / 20.0,  5.0 / 20.0,
  6.0 / 20.0,  7.0 / 20.0,  8.0 / 20.0,  9.0 / 20.0,  10.0 / 20.0, 11.0 / 20.0,
  12.0 / 20.0, 13.0 / 20.0, 14.0 / 20.0, 15.0 / 20.0, 16.0 / 20.0, 17.0 / 20.0,
  18.0 / 20.0, 19.0 / 20.0, 1.0
};

typedef struct {
  float minGC;
  float maxGC;
  int minTm;
  int maxTm;
  int deComplementarity;
  int homeopolymer;
  int avoidCGIn3;
  int avoidTIn3;
} FilterOpts;

static inline Array*
filterSegment(Array* segments, FilterOpts* opts)
{
  Array* result = arrayNew(10);
  unsigned char bases[KMER_LONG_LEN] = { 4 };
  for (size_t i = 0; i < segments->size; i++) {
    Segment* s = segments->data[i];
    debug("filter segment %zu", s->bases->size);
    size_t j = 0;
    size_t je = s->bases->size - KMER_LONG_LEN + 1;
    if (s->bases->size < KMER_LONG_LEN * 2) {
      continue;
    }
    for (; j < je; j++) {
      int gc = 0;
      for (size_t z = 0; z < KMER_LONG_LEN; z++) {
        bases[z] = bitarrayGet(s->bases, z + j);
        if (int2base[bases[z]] == 'G' || int2base[bases[z]] == 'C') {
          gc++;
        }
      }
      if (GC_RATE[gc] < opts->minGC || GC_RATE[gc] > opts->maxGC) {
        continue;
      }

      // compute Tm
      // Wallace formula: Tm = 64.9 +41*(yG+zC-16.4)/(wA+xT+yG+zC)
      // Wallace RB et al.(1979)Nucleic Acids Res 6 : 3543 - 3557,PMID 158748
      float tm = 64.9 + 41 * ((float)gc - 16.4) / (float)KMER_LONG_LEN;
      if (tm < opts->minTm || tm > opts->maxTm) {
        continue;
      }

      // avoid CG in 3' end more than 3 times
      if (opts->avoidCGIn3) {
        int cg = 0;
        for (int i = 0; i < 3; i++) {
          if (int2base[bases[i]] == 'G' || int2base[bases[i]] == 'C') {
            cg++;
          }
        }
        if (cg == 3) {
          continue;
        }
      }
      // avoid T in 3' end
      if (opts->avoidTIn3) {
        char left = int2base[bases[0]];
        char right = int2base[bases[KMER_LONG_LEN - 1]];
        if (left == 'T' || right == 'T' || left == 'A' || right == 'A') {
          continue;
        }
      }
      // homopolymer
      char prevc = bases[0];
      int maxHome = 1;
      int isHome = 0;
      for (int i = 1; i < KMER_LONG_LEN; i++) {
        char c = bases[i];
        if (c == prevc) {
          maxHome++;
        } else {
          maxHome = 1;
        }
        if (maxHome >= opts->homeopolymer) {
          isHome = 1;
          break;
        }
        prevc = c;
      }
      if (isHome) {
        continue;
      }

      Segment* segment = dmalloc(sizeof(Segment));
      segment->start = s->start + j;
      segment->end = s->start + j + KMER_LONG_LEN - 1;
      segment->bases = bitarrayNew(KMER_LONG_LEN, 2);
      segment->name = s->name;
      segment->vaild = 1;
      segment->Tm = tm;
      for (size_t k = 0; k < KMER_LONG_LEN; k++) {
        bitarraySet(segment->bases, k, bases[k]);
      }
      arrayPush(result, segment);
    }
  }
  return result;
}

static inline int
cmpSegmentVaild(const void* a, const void* b)
{
  Segment* sa = (Segment*)a;
  Segment* sb = (Segment*)b;
  return sb->vaild - sa->vaild;
}

static inline void
rankSegmentByCount(Array* segments, Query* query)
{
  Array* seqs = query->seqs;
  uint64_t longKmer = 0;
  uint64_t longKmerMask = (1ULL << (KMER_LONG_LEN * 2)) - 1ULL;
  unsigned char basemap[128] = { 4 };
  create_base2int(basemap);
  for (size_t i = 0; i < seqs->size; i++) {
    Seq* seq = seqs->data[i];
    for (size_t j = 0; j < seq->len - KMER_LONG_LEN + 1; j++) {
      longKmer <<= 2;
      longKmer &= longKmerMask;
      longKmer |= (basemap[seq->seq[j]] & 0x3);
      if (j < KMER_LONG_LEN) {
        continue;
      }
      // check if longKmer is in segment
      for (size_t k = 0; k < segments->size; k++) {
        Segment* s = segments->data[k];
        uint64_t skmer = 0;
        for (size_t z = 0; z < KMER_LONG_LEN; z++) {
          skmer <<= 2;
          skmer |= bitarrayGet(s->bases, z);
        }
        if (skmer == longKmer) {
          s->vaild++;
        }
      }
    }
  }
  for (size_t i = 0; i < segments->size; i++) {
    Segment* s = segments->data[i];
    s->vaild--;
  }
  // sort segments by vaild
  qsort(segments->data, segments->size, sizeof(Segment*), cmpSegmentVaild);
}

static inline Array*
createCircle(Array* segment, int count)
{
  size_t segmentSize = segment->size;
  Array* result = arrayNew(count);
  return NULL;
}

static inline void
freeSegments(Array* segments)
{
  // free segments
  for (size_t i = 0; i < segments->size; i++) {
    Segment* s = segments->data[i];
    bitarrayFree(s->bases);
    dfree(s, sizeof(Segment));
  }
  arrayFree(segments);
}

static inline void
writeSegments(Array* segments, const char* output)
{
  FILE* fp = fopen(output, "w");
  if (fp == NULL) {
    error("open file %s failed.", output);
    exit(1);
  }
  char buff1[100];
  char buff2[100];
  uint64_t kmer = 0;
  uint64_t reverseKmer = 0;
  uint64_t kmerMask = (1ULL << (KMER_LONG_LEN * 2)) - 1;
  uint64_t kmerShift = (KMER_LONG_LEN - 1) * 2;
  fprintf(fp, "id\tchr\tstart\tend\tTm\tkmer\treverse_kmer\n");
  for (size_t i = 0; i < segments->size; i++) {
    Segment* s = segments->data[i];
    for (int j = 0; j < KMER_LONG_LEN; j++) {
      kmer = (kmer << 2) & kmerMask;
      kmer |= bitarrayGet(s->bases, j);
      reverseKmer = (reverseKmer >> 2) & kmerMask;
      reverseKmer |= (3ULL - bitarrayGet(s->bases, j)) << kmerShift;
    }
    int2KmerString(kmer, KMER_LONG_LEN, buff1);
    int2KmerString(reverseKmer, KMER_LONG_LEN, buff2);
    fprintf(fp, "%ld\t%s\t%ld\t%ld\t%.2f\t%s\t%s\n", i, s->name, s->start + 1,
            s->end + 1, s->Tm, buff1, buff2);
  }
  fclose(fp);
}

int
main(int argc, const char* argv[])
{
  if (argc != 4) {
    fprintf(stderr, "Usage: %s <index> <query> <output>\n", argv[0]);
    exit(1);
  }
  const char* index_path = argv[1];
  const char* query_path = argv[2];
  const char* output_path = argv[3];
  log_set_level(PGLOG_LEVEL_DEBUG);
  Index* index = NULL;
  debug("start.");
  char buff[1024] = { 0 };
  sprintf(buff, "%s.index", index_path);
  int exist = isFileExist(buff);
  if (exist) {
    index = loadIndex(buff);
  } else {
    index = createIndex(index, index_path);
    dumpIndex(index, buff);
  }
  debug("index created.");
  Query* query = createQuery(query_path);
  debug("query created.");
  vaildKmers(query, index);
  debug("kmers vailded.");
  Array* segments = collectSegment(query);
  debug("segments collected.");
  FilterOpts filterOpts = {
    .avoidCGIn3 = 1,
    .avoidTIn3 = 1,
    .minGC = 0.45,
    .maxGC = 0.55,
    .minTm = 50,
    .maxTm = 60,
    .deComplementarity = 1,
    .homeopolymer = 4,
  };
  Array* filtered = filterSegment(segments, &filterOpts);
  rankSegmentByCount(filtered, query);
  debug("segments filtered.");
  writeSegments(filtered, output_path);
  freeSegments(segments);
  freeSegments(filtered);
  freeQuery(query);
  freeIndex(index);
  debug("memory usage: %lu", getUsedMemory());
}
