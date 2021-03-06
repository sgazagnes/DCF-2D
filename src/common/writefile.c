#include <stddef.h>
#include <stdio.h>
#include <errno.h>
#include <stdlib.h>
#include "writefile.h"
#include "types.h"
#include "checks.h"
#include "queue.h"
#include "logc.h"

void write_boundary_file_ascii(Boundary *b, const char *filename) {
  FILE *fp = fopen(filename, "w");
  check_not_null(fp, errno);

  fprintf(fp, "array: %p\n", (void*) b->array);
  fprintf(fp, "size: %zu\n", b->size);
  fprintf(fp, "offsets: %zu %zu %zu %zu %zu\n", b->offset[0], b->offset[1], b->offset[2], b->offset[3], b->offset[4]);
  fprintf(fp, "+ b-index   +   index   +   parent  +  bparent  +   area    + gval + filt + flags +\n");
  size_t size = b->size;
  for (size_t i = 0; i < size; i++) {
    fprintf(fp, "|%11ld", b->array[i].borderindex);
    fprintf(fp, "|%11ld", b->array[i].index);
    fprintf(fp, "|%11ld", b->array[i].parent);
    fprintf(fp, "|%11ld", b->borderparent[i].i);
    fprintf(fp, "|%11lu", b->array[i].area);
    fprintf(fp, "|%5d", b->array[i].gval);
    fprintf(fp, "|%6d", b->array[i].filter);
    // fprintf(fp, "|%7x", b->array[i].reached);

    fprintf(fp, "|\n");
  }
  fprintf(fp, "\n");

  int err = fclose(fp);
  check_file_close(err, filename);
} /* write_boundary_file_ascii */

void write_boundary_file_binary(Boundary *b, const char *filename) {
  FILE *fp = fopen(filename, "wb");
  check_not_null(fp, errno);

  fprintf(fp, "%zu\n", b->size);
  fprintf(fp, "%zu %zu %zu %zu %zu\n", b->offset[0], b->offset[1], b->offset[2], b->offset[3], b->offset[4]);
  size_t size = b->size;
  for (size_t i = 0; i < size; i++) {
    fprintf(fp, "%ld ", b->array[i].borderindex);
    fprintf(fp, "%ld ", b->array[i].index);
    fprintf(fp, "%ld ", b->array[i].parent);
    fprintf(fp, "%11ld ", b->borderparent[i].i);
    fprintf(fp, "%lu ", b->array[i].area);
    fprintf(fp, "%d ", b->array[i].gval);
    fprintf(fp, "%d ", b->array[i].filter);
    //  fprintf(fp, "%x ", b->array[i].reached);
  }
  fprintf(fp, "\n");

  int err = fclose(fp);
  check_file_close(err, filename);
} /* write_boundary_file_binary */

void write_maxtree_file_ascii(MaxNode *m, size_t size, const char *filename) {
  FILE *fp = fopen(filename, "w");
  check_not_null(fp, errno);

  fprintf(fp, "size: %zu\n", size);
  fprintf(fp, "+ b-index   +   index   +   parent  +   area    + gval + filt + reached +\n");

  for (size_t i = 0; i < size; i++) {
    fprintf(fp, "|%11ld", m[i].borderindex);
    fprintf(fp, "|%11ld", m[i].index);
    fprintf(fp, "|%11ld", m[i].parent);
    fprintf(fp, "|%11lu", m[i].area);
    fprintf(fp, "|%5d", m[i].gval);
    fprintf(fp, "|%6d", m[i].filter);
    //  fprintf(fp, "|%7x", m[i].reached);
    fprintf(fp, "|\n");
  }
  fprintf(fp, "\n");

  int err = fclose(fp);
  check_file_close(err, filename);
} /* write_boundary_file_ascii */

void queues_to_file(Queue *hq, value maxvalue, const char *filename) {
  FILE *fp = fopen(filename, "w");
  check_not_null(fp, errno);

  for (value i = 0; i < maxvalue; i++) {
    fprintf(fp, "hq[%d]: %lu pixels head: %lu tail: %lu \n", i, hq[i].pixels[hq[i].head], hq[i].head, hq[i].tail);
  }

  int err = fclose(fp);
  check_file_close(err, filename);
}
