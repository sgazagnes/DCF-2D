#include <stddef.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include "queue.h"
#include "checks.h"
#include "logc.h"
#include "types.h"

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


pStack *pStackCreate(ulong maxsize, ulong *stackArray){
  pStack *newStack = malloc(sizeof(pStack));
  check_alloc(newStack, 002);
  newStack->array = stackArray;
  newStack->maxsize= maxsize;
  newStack->curpos = maxsize;

  return newStack;
}

void pStackDelete(pStack *oldstack){
  free(oldstack);
}


pQueue *pQueueCreate(ulong maxsize){
  pQueue *newQueue = (pQueue *) malloc(sizeof(pQueue));
  newQueue->size = 0;
  newQueue->array = malloc((maxsize+1)*sizeof(long));
  newQueue->maxsize=maxsize;
  return newQueue;
}


void pQueueDelete(pQueue *oldqueue){
  free(oldqueue->array);
  free(oldqueue);

}

int pQueuePop(pQueue *queue, ushort *priority){
  ulong outval = queue->array[1];
  ulong current = 1, moved;
  ushort curval;

  moved = queue->array[queue->size];
  queue->size--;
  curval = priority[moved];

  while ( ((current*2<=queue->size) &&
	   (curval< priority[queue->array[current*2]]))
	  ||
	  ((current*2+1<=queue->size) &&
	   (curval< priority[queue->array[current*2+1]]))
	  ){
    if ((current*2+1<=queue->size) && 
	(priority[queue->array[current*2]]< 
	 priority[queue->array[current*2+1]])){
      queue->array[current]= queue->array[current*2+1];
      current+=current+1;
    } else {
      queue->array[current]= queue->array[current*2];
      current+=current;
    }
  }
  queue->array[current]=moved;
  return outval;
}

void pQueuePush(pQueue *queue, ushort *priority, ulong pixpos){
  long current;
  ushort curval = priority[pixpos];
  queue->size++;
  current=queue->size;
  
  while ((current/2 !=0) && (priority[queue->array[current/2]]<curval)){
    queue->array[current]= queue->array[current/2];
    current=current/2;
  }

  queue->array[current]=pixpos;
}
