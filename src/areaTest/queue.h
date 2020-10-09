#ifndef QUEUE_H
#define QUEUE_H
#include "types.h"

#define queue_first(hq,h)    (hq[h].pixels[hq[h].head++]) // dequeue
#define queue_add(hq,h,p)    (hq[h].pixels[hq[h].tail++] = p)
#define queue_is_not_empty(hq,h) (hq[h].head != hq[h].tail)

Queue *create_queue(size_t imgsize, unsigned long maxvalue);
void set_queue_offsets(Queue *hq, unsigned long *numpixelsperlevel, unsigned long maxvalue); // set offsets
void free_queue(Queue *hq);

#endif
