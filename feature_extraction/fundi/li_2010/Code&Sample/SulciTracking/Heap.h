#ifndef HEAP_H
#define HEAP_H

#ifdef __cplusplus
extern "C" {
#endif
 
#include "pgutil.h"

typedef int* XhTrace;

typedef struct {
   double value;
   int id;
   int*  p;
} XheapElement;

typedef PGlist Xheap;

int   xhSize(Xheap H);   /* get size of heap */
Xheap xhInitEmpty();     /* an empty heap */
Xheap xhInit(XheapElement *array, int N); /* init from an array (0,N-1) */
void  xhDestroy(Xheap H); /* destroy the heap and free the memory */
int   xhUpHeap(int k, Xheap H);
int   xhDownHeap(int k, Xheap H);
int   xhInsert(double value, int id, int *p, Xheap H);
XheapElement xhRemove(Xheap H);
XheapElement xhReplace(double value, int id, int *p, Xheap H);
XheapElement xhDelete(int k, Xheap H);
XheapElement xhChange(int k, double value, int id, int *p, Xheap H);
XheapElement xhChangeValue(int k, double value, Xheap H);
XheapElement xhGet(int k, Xheap H); /* k must be 1, 2, ... N */
#define xhIsEmpty(H)  (xhSize(H) == 0)

#ifdef __cplusplus
}
#endif

#endif
