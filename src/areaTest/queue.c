#include <stddef.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include "queue.h"
#include "checks.h"
#include "logc.h"

Queue *create_queue(size_t imgsize, unsigned long maxvalue) {
  Queue *hq;
  hq = calloc(maxvalue, sizeof(Queue));
  check_alloc(hq, 120);
  hq->pixels = calloc(imgsize, sizeof(unsigned long));
  check_alloc(hq->pixels, 119);
  return hq;
}

void set_queue_offsets(Queue *hq, unsigned long *numpixelsperlevel, unsigned long maxvalue) {
  ulong i;
  hq->head = hq->tail = 0;

  for (i = 1; i < maxvalue; i++) {
    /* set index of queue corresponding to level i to an offset (sum of number of pixels in all previous levels: n-2 levels plus level n-1) */
    /* trace("hq[%d].pixels %p\n", i, (void *)hq[i-1].pixels); */
    hq[i].pixels = hq[i - 1].pixels + numpixelsperlevel[i - 1];
    hq[i].head = hq[i].tail = 0;
  }
}

void free_queue(Queue *hq) {
  free(hq->pixels);
  free(hq);
}
