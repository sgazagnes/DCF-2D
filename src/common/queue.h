#ifndef QUEUE_H
#define QUEUE_H
#include "types.h"

#define queue_first(hq,h)    (hq[h].pixels[hq[h].head++]) // dequeue
#define queue_add(hq,h,p)    (hq[h].pixels[hq[h].tail++] = p)
#define queue_is_not_empty(hq,h) (hq[h].head != hq[h].tail)
#define IsEmptyStack(stack)        ((stack->curpos)==stack->maxsize)
#define IsFullStack(stack)         ( stack->curpos == 0 )
#define IsEmpty(queue)             ((queue->size)==0)
#define IsFull(queue)              ( queue->size == queue->maxsize )
#define pStackTop(stack)       ( stack->array[stack->curpos + 1])
#define pStackPop(stack)       (stack->array[++stack->curpos])
#define pStackPush(stack,elem)      (stack->array[stack->curpos--]=elem)
#define pQueueFront(queue)       (queue->array[1])

Queue *create_queue(size_t imgsize, unsigned long maxvalue);
void set_queue_offsets(Queue *hq, unsigned long *numpixelsperlevel, unsigned long maxvalue); // set offsets
void free_queue(Queue *hq);
pStack *pStackCreate(ulong maxsize, ulong *stackArray);
void pStackDelete(pStack *oldstack);
pQueue *pQueueCreate(ulong maxsize);
void pQueueDelete(pQueue *oldqueue);
int pQueuePop(pQueue *queue, ushort *priority);
void pQueuePush(pQueue *queue, ushort *priority, ulong pixpos);
void pQueuePush(pQueue *queue, ushort *priority, ulong pixpos);

#endif
