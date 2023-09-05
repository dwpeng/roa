#include "alloc.h"
#include "arg.h"
#include "array.h"
#include "bitarray.h"
#include "file.h"
#include "log.h"
#include "seq.h"

#include <stdarg.h>
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
#define KMER_LONG_MASK (1ULL << 40) - 1
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
  uint64_t kmer : 32;
  uint64_t reverse_kmer : 32;
} Kmer;

typedef struct {
  Array* kmers; // Array<Array<Kmer>>
  Array* seqs;  // Array<Seq*>
} Query;

typedef struct Segment Segment;
struct Segment {
  int id;          // id of the segment
  char* name;      // name of the segment
  size_t start;    // start position of the segment
  size_t end;      // end position of the segment
  BitArray* bases; // bases of the segment
  int vaild;       // if the segment is vaild
  float Tm;        // pcr melting temperature
};

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
createQuery(const char* path, Index* index)
{
  uint32_t kmer = 0x00000000;
  uint32_t reverse_kmer = 0x00000000;
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
          reverse_kmer = 0x00000000;
          count = 0;
          continue;
        }
        kmer = (kmer << 2) | (c & 0x3) & KMER_MASK;
        reverse_kmer =
            (reverse_kmer >> 2) | ((0x00000003 - c) << 30) & KMER_MASK;
        count++;
        if (count < KMER_LEN) {
          continue;
        }
        Kmer* kmer_t = dmalloc(sizeof(Kmer));
        kmer_t->kmer = kmer;
        kmer_t->reverse_kmer = reverse_kmer;
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
      if (bitarrayGet(index->index, kmer->reverse_kmer)) {
        kmer->drop = 1;
        continue;
      }
    }
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
  float minTm;
  float maxTm;
  int deComplementarity;
  int homeopolymer;
  int avoidCGIn3;
  int avoidTIn3;
  int ncircle;
} FilterOpts;

static inline Array*
filterSegment(Array* segments, FilterOpts* opts)
{
  Array* result = arrayNew(10);
  unsigned char bases[KMER_LONG_LEN] = { 4 };
#ifdef parallel
#pragma omp parallel for
#endif
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
#ifdef parallel
#pragma omp critical
#endif
      // set id
      segment->id = result->size;
      arrayPush(result, segment);
    }
  }
  return result;
}

int
cmpSegments(const void* a, const void* b)
{
  Segment* sa = *(Segment**)a;
  Segment* sb = *(Segment**)b;
  return sb->vaild - sa->vaild;
}

static inline Array*
pairJoinCheck(Array* segments, Index* index)
{
  Array* pair = arrayNew(segments->size);
  // connect every two segments and check if the connection is vaild
  // if vaild, then join the two segments, insert the pair into the pair array
  // connect
#ifdef parallel
#pragma omp parallel for
#endif
  for (size_t i = 0; i < segments->size; i++) {
    Segment* s1 = segments->data[i];
    arrayPush(pair, bitarrayNew(segments->size, 1));
    for (size_t j = 0; j < segments->size; j++) {
      if (i == j) {
        continue;
      }
      Segment* s2 = segments->data[j];
      uint64_t kmer = 0ULL;
      uint64_t reverseKmer = 0ULL;
      int c = 0;
      int succ = 1;
      for (int k = 1; k < s1->bases->size + s2->bases->size - 1; k++) {
        char base = 4;
        if (k >= s1->bases->size) {
          base = bitarrayGet(s2->bases, k - s1->bases->size);
        } else {
          base = bitarrayGet(s1->bases, k);
        }
        kmer = ((kmer << 2) | base) & KMER_MASK;
        reverseKmer =
            (reverseKmer >> 2) | ((3ULL - base) << ((KMER_LEN - 1) * 2));
        c++;
        if (c < KMER_LEN) {
          continue;
        }
        // query index
        if (bitarrayGet(index->index, kmer)
            || bitarrayGet(index->index, reverseKmer)) {
          succ = 0;
          break;
        }
      }
      if (succ) {
        s1->vaild++;
        bitarraySet(arrayGet(pair, i), j, 1);
      }
    }
  }
  // remove all segments that have circle deps
  for (int i = 0; i < segments->size; i++) {
    for (int j = 0; j < segments->size; j++) {
      if (i == j) {
        continue;
      }
      if (bitarrayGet(arrayGet(pair, i), j)
          && bitarrayGet(arrayGet(pair, j), i)) {
        bitarraySet(arrayGet(pair, j), i, 0);
        ((Segment*)segments->data[i])->vaild--;
      }
    }
  }
  qsort(segments->data, segments->size, sizeof(Segment**), cmpSegments);
  return pair;
}

static inline void
printPair(Array* pair, const char* path)
{
  FILE* fp = fopen(path, "w");
  if (fp == NULL) {
    error("open file %s failed.", path);
    exit(1);
  }
  for (size_t i = 0; i < pair->size; i++) {
    for (size_t j = 0; j < pair->size; j++) {
      fprintf(fp, "%d ", bitarrayGet(arrayGet(pair, i), j));
    }
    fprintf(fp, "\n");
  }
  fclose(fp);
}

static inline void
freePairArray(Array* pair)
{
  for (size_t i = 0; i < pair->size; i++) {
    bitarrayFree(pair->data[i]);
  }
  arrayFree(pair);
}

#define checkJoin(pair, s1, s2)                                               \
  (pair == NULL ? 1 : bitarrayGet(arrayGet(pair, (s1)->id), (s2)->id))

typedef struct Node Node;
struct Node {
  Segment* segment;
  int next_count;
  Node** next;
};

static inline Array*
createCircle(Array* segments, Array* pair, int count)
{
  size_t segmentSize = segments->size;
  Array* result = arrayNew(count * KMER_PER_CIRCLE);
  int offset = 0;
  if (pair) {
    Segment *a, *b;
    int span = 1000;
    Node** nodes = dmalloc(sizeof(Node*) * segmentSize);
    for (size_t i = 0; i < segmentSize; i++) {
      nodes[i] = dmalloc(sizeof(Node));
      nodes[i]->segment = segments->data[i];
      nodes[i]->next_count = 0;
      nodes[i]->next = NULL;
    }

    for (int i = 0; i < segmentSize; i++) {
      a = segments->data[i];
      nodes[i]->next = dmalloc(sizeof(Node*) * segmentSize);
      for (int j = 0; j < segmentSize; j++) {
        if (i == j) {
          continue;
        }
        b = segments->data[j];
        if (checkJoin(pair, a, b)) {
          if (b->start - a->start < span || a->start - b->start < span) {
            continue;
          }
          nodes[i]->next[nodes[i]->next_count] = nodes[j];
          nodes[i]->next_count++;
        }
      }
    }

    int circle_count = 0;
    // dfs
    Node** stack = dmalloc(sizeof(Node*) * segmentSize);
    int stack_size = 0;
    int* visited = dmalloc(sizeof(int) * segmentSize);
    for (int i = 0; i < segmentSize; i++) {
      visited[i] = 0;
    }

    for (int i = 0; i < segmentSize; i++) {
      if (visited[i]) {
        continue;
      }
      if (circle_count >= count) {
        break;
      }
      stack[stack_size] = nodes[i];
      stack_size++;
      while (stack_size > 0) {
        if (circle_count >= count) {
          break;
        }
        Node* node = stack[stack_size - 1];
        visited[node->segment->id] = 1;
        if (node->next_count == 0) {
          // find a circle
          for (int i = 0; i < KMER_PER_CIRCLE; i++) {
            arrayPush(result, node->segment);
            visited[node->segment->id] = 1;
          }
          circle_count++;
          stack_size = 0;
          continue;
        }
        if (stack_size == KMER_PER_CIRCLE) {
          // find a circle
          for (int i = 0; i < KMER_PER_CIRCLE; i++) {
            arrayPush(result, stack[i]->segment);
            visited[stack[i]->segment->id] = 1;
          }
          circle_count++;
          stack_size = 0;
          continue;
        }
        int find = 0;
        for (int i = 0; i < node->next_count; i++) {
          if (visited[node->next[i]->segment->id]) {
            continue;
          }
          stack[stack_size] = node->next[i];
          stack_size++;
          find = 1;
          break;
        }
        if (find) {
          continue;
        }
        stack_size--;
      }
    }

    // free
    for (int i = 0; i < segmentSize; i++) {
      dfree(nodes[i]->next, sizeof(Node*) * segmentSize);
      dfree(nodes[i], sizeof(Node));
    }
    dfree(nodes, sizeof(Node*) * segmentSize);
    dfree(stack, sizeof(Node*) * segmentSize);
    dfree(visited, sizeof(int) * segmentSize);

  } else {
    for (int i = 0; i < count; i++) {
      // generate a random circle
      for (int j = 0; j < KMER_PER_CIRCLE; j++) {
        // find_circle;
        arrayPush(result, segments->data[offset]);
        offset++;
      }
    }
  }

  return result;
}

#define segmentToKmer(segment, kstr, rstr)                                    \
  do {                                                                        \
    uint64_t kmer = 0;                                                        \
    uint64_t reverseKmer = 0;                                                 \
    for (int i = 0; i < KMER_LONG_LEN; i++) {                                 \
      kmer = (kmer << 2) | bitarrayGet(segment->bases, i);                    \
      reverseKmer = (reverseKmer >> 2)                                        \
                    | ((3ULL - bitarrayGet(segment->bases, i))                \
                       << ((KMER_LONG_LEN - 1) * 2));                         \
    }                                                                         \
    int2KmerString(kmer, KMER_LONG_LEN, (kstr));                              \
    int2KmerString(reverseKmer, KMER_LONG_LEN, (rstr));                       \
  } while (0)

static inline void
saveCircle(Array* circle, int count, const char* outpath)
{
  FILE* fp = fopen(outpath, "w");
  char circle_template[100] = { 0 };
  char kmer_str[100] = { 0 };
  char reverseKmer_str[100] = { 0 };
  int circle_id = 1;
  int circle_sub_id = 0;
  int offset = 0;
  int max_count = (circle->size + 3) / KMER_PER_CIRCLE;
  if (count > max_count) {
    info("count %d is larger than max count %d", count, max_count);
    info("set count to %d", max_count);
    count = max_count;
  }
  for (int i = 0; i < count * KMER_PER_CIRCLE; i++) {
    Segment* s = (Segment*)circle->data[i];
    segmentToKmer(s, kmer_str, reverseKmer_str);
    fprintf(fp, ">probe-%d/%d %s:%ld\n%s\n", circle_id, circle_sub_id + 1,
            s->name, s->start, reverseKmer_str);
    circle_sub_id++;
    for (int j = 0; j < KMER_LONG_LEN; j++) {
      // copy from reverseKmer_str
      circle_template[offset] = reverseKmer_str[j];
      offset++;
    }
    if (circle_sub_id == KMER_PER_CIRCLE) {
      fprintf(fp, ">circle-%d\n%s\n", circle_id, circle_template);
      info("save circle %d/%d", circle_id, count);
      circle_id++;
      circle_sub_id = 0;
      offset = 0;
    }
  }
  fclose(fp);
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
  fprintf(fp, "id\tchr\tstart\tend\tTm\tkmer\treverse_kmer\tcount\n");
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
    fprintf(fp, "%d\t%s\t%ld\t%ld\t%.2f\t%s\t%s\t%d\n", s->id, s->name,
            s->start + 1, s->end + 1, s->Tm, buff1, buff2, s->vaild);
  }
  fclose(fp);
}

#define p(...)                                                                \
  do {                                                                        \
    fprintf(stderr, __VA_ARGS__);                                             \
  } while (0)

void
design_usage()
{
  p("ROA Template Designer.\n");
  p("Usage:\n");
  p("  ./roa design <options>\n");
  p("Example:\n");
  p("  ./roa design -i index.index -q query.fa -o template.fa -pairCheck 1\n");
  p("Options:\n");
  p("  -i <index>    index file path\n");
  p("  -q <query>    query file path\n");
  p("  -o <output>   output file path [template.fa]\n");
  p("  -homopolymer  homopolymer length [3]\n");
  p("  -minGC        min GC rate [0.45]\n");
  p("  -maxGC        max GC rate [0.55]\n");
  p("  -minTm        min melting temperature [52.4]\n");
  p("  -maxTm        max melting temperature [55.4]\n");
  p("  -avoidCGIn3   avoid CG in 3' end [1]\n");
  p("  -avoidTIn3    avoid T in 3' end [1]\n");
  p("  -pairCheck    check pair [0] maybe cost a long time\n");
  p("  -ncircle      number of circles [5]\n");
  p("  -h            show this help message\n");
}

void
index_usage()
{
  p("ROA Template Designer.\n");
  p("Usage:\n");
  p("  ./roa index <index> <fa>...\n");
  p("Example:\n");
  p("  ./roa index index.index ref1.fa ref2.fa ...\n");
  p("Options:\n");
  p("  -h            show this help message\n");
}

int
usage(int argc, char* argv[])
{
  p("ROA Template Designer.\n");
  p("Usage:\n");
  p("  ./roa <command> <options>\n");
  p("Commands:\n");
  p("  index         create index file\n");
  p("  design        design ROA template\n");
  return 0;
}

int
invoke_help(int argc, char* argv[])
{
  for (int i = 0; i < argc; i++) {
    if (strcmp(argv[i], "-h") == 0) {
      return 1;
    }
  }
  return 0;
}

arginit(do_design)
{
  if (invoke_help(argc, argv)) {
    design_usage();
    exit(1);
  }
  if (argc < 2) {
    design_usage();
    exit(1);
  }
  const char* index_path = NULL;
  const char* query_path = NULL;
  const char* output_path = "template.fa";
  int homopolymer = 3;
  float minGC = 0.45;
  float maxGC = 0.55;
  float minTm = 52.4;
  float maxTm = 55.4;
  int avoidCGIn3 = 1;
  int avoidTIn3 = 1;
  int ncircle = 5;
  int pairCheck = 0;
  argstart()
  {
    argpass("-h");
    argstring("-i", index_path);
    argstring("-q", query_path);
    argstring("-o", output_path);
    argint("-homopolymer", homopolymer);
    argfloat("-minGC", minGC);
    argfloat("-maxGC", maxGC);
    argfloat("-minTm", minTm);
    argfloat("-maxTm", maxTm);
    argbool("-avoidCGIn3", avoidCGIn3);
    argbool("-avoidTIn3", avoidTIn3);
    argint("-ncircle", ncircle);
    argbool("-pairCheck", pairCheck);
    argend();
  }
  log_set_level(PGLOG_LEVEL_DEBUG);
  if (index_path == NULL || query_path == NULL || output_path == NULL) {
    design_usage();
    exit(1);
  }
  info("index_path: %s", index_path);
  info("query_path: %s", query_path);
  info("output_path: %s", output_path);
  info("homopolymer: %d", homopolymer);
  info("minGC: %.2f", minGC);
  info("maxGC: %.2f", maxGC);
  info("minTm: %.2f", minTm);
  info("maxTm: %.2f", maxTm);
  info("avoidCGIn3: %d", avoidCGIn3);
  info("avoidTIn3: %d", avoidTIn3);
  info("ncircle: %d", ncircle);
  info("pairCheck: %d", pairCheck);
  Index* index = loadIndex(index_path);
  Query* query = createQuery(query_path, index);
  vaildKmers(query, index);
  Array* segments = collectSegment(query);
  FilterOpts filterOpts = { .avoidCGIn3 = avoidCGIn3,
                            .avoidTIn3 = avoidTIn3,
                            .minGC = minGC,
                            .maxGC = maxGC,
                            .minTm = minTm,
                            .maxTm = maxTm,
                            .deComplementarity = 1,
                            .homeopolymer = homopolymer,
                            .ncircle = ncircle };
  Array* filtered = filterSegment(segments, &filterOpts);
  freeSegments(segments);
  debug("fileter %zu segments", filtered->size);
  if (filtered->size) {
    Array* pair = NULL;
    if (pairCheck) {
      pair = pairJoinCheck(filtered, index);
    }
    Array* circles = createCircle(filtered, pair, ncircle);

    saveCircle(circles, ncircle, output_path);
    if (pairCheck) {
      freePairArray(pair);
    }
    arrayFree(circles);
  } else {
    info("no specific kmer found.");
  }
  freeSegments(filtered);
  freeQuery(query);
  freeIndex(index);
  debug("useMemory: %zu", getUsedMemory());
}

void
do_index(int argc, char* argv[])
{
  if (invoke_help(argc, argv)) {
    index_usage();
    exit(1);
  }
  if (argc < 2) {
    index_usage();
    exit(1);
  }
  const char* index_path = argv[0];
  Index* index = NULL;
  for (int i = 1; i < argc; i++) {
    const char* ref_path = argv[i];
    if (!isFileExist(ref_path)) {
      error("file %s not exist.", ref_path);
      continue;
    }
    info("indexing %s", ref_path);
    char buff[1024] = { 0 };
    sprintf(buff, "%s.index", ref_path);
    int exist = isFileExist(buff);
    if (exist) {
      if (index == NULL) {
        index = loadIndex(buff);
      } else {
        Index* tmp = loadIndex(buff);
        bitarrayOr(index->index, tmp->index);
        freeIndex(tmp);
      }
    } else {
      Index* tmp = createIndex(NULL, ref_path);
      dumpIndex(tmp, buff);
      if (index == NULL) {
        index = tmp;
      } else {
        bitarrayOr(index->index, tmp->index);
        freeIndex(tmp);
      }
    }
  }
  info("Saving to %s", index_path);
  dumpIndex(index, index_path);
  freeIndex(index);
}

int
main(int argc, char* argv[])
{
  log_set_level(PGLOG_LEVEL_DEBUG);
  if (argc < 2) {
    usage(argc, argv);
    return 0;
  }
  if (strcmp(argv[1], "index") == 0) {
    do_index(argc - 2, argv + 2);
    return 0;
  }
  if (strcmp(argv[1], "design") == 0) {
    do_design(argc - 2, argv + 2);
    return 0;
  }
  usage(argc, argv);
  return 0;
}
