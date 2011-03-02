/*============================================================================
//                        PGutil Summary
//============================================================================
//
//  DATA STRUCTURE:
//  PGlist ---> {elementSize, size, data, capacity, capacityIncrement}
//  
//  NOTE: the type of PGlist is a pointer.
//
//  Assume the following variable declarations:
//      PGlist list, list1, list2;
//      ListElement element;
//      ListElement *data;
//      int capacity, capacityIncrement, size, 
//          elementSize, index, returnFlag;
//
//  LIST MANIPULATIONS:
//
//  pgList ........................ default constructor
//           list  = pgList(sizeof(ListElement));
//
//  pgList1 ....................... default constructor with capacity
//           list1 = pgList(sizeof(ListElement), capacity);
//
//  pgList2 ....................... default constructor with capacity 
//                                  and capacity increment
//           list2 = pgList(sizeof(ListElement), capacity, capacityIncrement);
//
//  pgListOfSize .................. default constructor with a specified size
//           list  = pgListOfSize(size, sizeof(ListElement));
//
//  pgListDelete .................. delete the list
//           pgListDelete(list);
//
//  pgListAddElement .............. add an element to this list
//           pgListAddElement(list, &element);
//
//  pgListInsertElementAt ......... insert an element in the list at the given 
//                                  index. Each list element's index greater 
//                                  or equal to the specified index is
//                                  shifted upward than its previous value.
//           returnFlag = pgListInsertElementAt(list, index, &element);
//
//  pgListElementAt ............... retrieve an element at index
//           returnFlag = pgListElementAt(list, index, &element);
//
//  pgListSetElementAt ............ set the element at the specified index of 
//                                  this list by copying the value of 
//                                  given element.
//           returnFlag = pgListSetElementAt(list, index, &element);
//
//  pgListRemoveElementAt ......... Delete the element at the specified index.
//                                  The index of each element after the 
//                                  specified index is decreased by 1.
//           returnFlag = pgListRemoveElementAt(list, index);
//
//  pgListRemoveAllElements ....... removes all elements from this list
//                                  and sets its size to zero 
//           pgListRemoveAllElements(list);
//
//  pgListTrim .................... trim this list to its current size 
//           pgListTrim(list);
//
//  pgListIsEmpty ................. PG_TRUE if this list has no elements
//           returnFlag = pgListIsEmpty(list);
//
//  pgListElementSize ............. the element size of each component
//           elementSize = pgListElementSize(list);
//
//  pgListSize .................... the current size of this list
//           size = pgListSize(list);
//
//  pgListData .................... the data of this list
//           data = pgListData(list);
//
//  DATA STRUCTURE:
//  PGstack ---> implemented on top of PGlist
//
//  NOTE: the type of PGstack is a pointer.
//
//  STACK MANIPULATIONS:
//  pgStack ....................... default constructor
//  pgStackPush ................... push an element into the stack
//  pgStackPop .................... pop out the top element from the stack
//  pgStackIsEmpty ................ PG_TRUE if this stack is empty
//  pgStackRemoveAllElements ...... removes all elements from the stack
//                                  and sets its size to zero
//  pgStackTrim ................... trim the stack to its current size
//  pgStackDelete ................. delete the stack
//  pgStackElementSize ............ the element size of each component
//
//  QUEUE DATA STRUCTURE:
//  PGqueue ---> implemented on top of PGlist
//
//  QUEUE  MANIPULATIONS:
//  pgQueue ....................... default constructor
//  pgQueue1 ...................... constructor 1
//  pgQueue2 ...................... constructor 2
//  pgQueueDelete ................. delete the queue and its memory
//  pgQueueRemoveAllElements ...... removes all elements from the queue
//                                  without releasing memory
//  pgQueuePush ................... push an element to the end of the queue
//  pgQueuePop .................... pop out the first element from the queue
//  pgQueueTrim ................... trim the queue to its current size
//  pgQueueToArray ................ extract an array from the queue
//  pgQueueSize ................... the current size of the queue
//  pgQueueIsEmpty ................ PG_TRUE if this stack is empty
//  pgQueueInfo ................... print the queue information
//  pgQueueElementSize ............ returns the size of an element
//==========================================================================*/

#ifndef PGUTIL_TOOLS
#define PGUTIL_TOOLS
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include "pgstd.h"

typedef struct {
   int   elementSize;
   int   size;
   void *data;
   int   capacity;
   int   capacityIncrement;
} PGlistStruct;

typedef PGlistStruct *PGlist;
#define PGstack PGlist

typedef struct {
   int start;
   int end;
   PGlist list;
} PGqueueStruct;

typedef PGqueueStruct *PGqueue;

#define PG_QUEUE_Q    2


#define pgListIsEmpty(list)            ((list)->size == 0 ? PG_TRUE: PG_FALSE)
#define pgListSize(list)               ((list)->size)
#define pgListData(list)               ((list)->data)
#define pgListElementSize(list)        ((list)->elementSize)

PGlist  pgList(int elementSize);
PGlist  pgList1(int elementSize, int capacity);
PGlist  pgList2(int elementSize, int capacity, int capacityIncrement);
PGlist  pgListOfSize(int size, int elementSize);
void    pgListDelete(PGlist list);
void    pgListAddElement(PGlist list, void *element);
void    pgListAddInt(PGlist list, int element);
void    pgListAddArray(PGlist list, void *array, int num);
PGerror pgListInsertElementAt(PGlist list, int index, void *element);
PGerror pgListSetElementAt(PGlist list, int index, void *element);
PGerror pgListElementAt(PGlist list, int index, 
                        /* stores the result at */ void *element);
PGerror pgListRemoveElementAt(PGlist list, int index);
void    pgListRemoveAllElements(PGlist list);
void    pgListTrim(PGlist list);
void    pgListInfo(PGlist list);

#define pgStack(elementSize)       pgList(elementSize)
#define pgStack1(elementSize,capacity)  pgList(elementSize,capacity)
#define pgStack2(elementSize,capacity,capacityIncrement) pgList(elementSize,capacity,capacityIncrement)
#define pgStackPush(stack, element) pgListAddElement(stack, element)
void    pgStackPop(PGstack stack, /* stores the result at */ void *element);
#define pgStackIsEmpty(stack)       pgListIsEmpty(stack)
#define pgStackRemoveAllElements(stack) pgListRemoveAllElements(stack)
#define pgStackTrim(stack)          pgListTrim(stack)
#define pgStackDelete(stack)        pgListDelete(stack)
#define pgStackElementSize(stack)  (stack->elementSize)

PGqueue pgQueue2(int elementSize, int capacity, int capacityIncrement);
PGqueue pgQueue1(int elementSize, int capacity);
PGqueue pgQueue(int elementSize);
void    pgQueueDelete(PGqueue queue);
void    pgQueueRemoveAllElements(PGqueue queue);
void    pgQueuePush(PGqueue queue, void* element);
void    pgQueuePushArray(PGqueue queue, void* array, int num);
PGerror pgQueuePop(PGqueue queue, void* element);
void    pgQueueTrim(PGqueue queue);
void*   pgQueueToArray(PGqueue queue);
void    pgQueueInfo(PGqueue queue);
#define pgQueueSize(queue)         (queue->end - queue->start + 1)
#define pgQueueIsEmpty(queue)      (pgQueueSize(queue) == 0) 
#define pgQueueElementSize(queue)  (pgListElementSize(queue->list))

/*============================================================================
//                        PGutil private functions
//                        USER PLEASE DO NOT USE
// NOTE: these numbers should not be changed once the list is created
//============================================================================
//  pgListCapacity ................ the current capacity of this list
//  pgListCapacityIncrement ....... the current capacityIncrement of this list
//  pgListSetElementSize .......... set the element size of this list
//  pgListSetSize ................. set the size of this list
//  pgListSetData ................. set the data of this list
//  pgListSetCapacity ............. set the capacity of this list
//  pgListSetCapacityIncrement .... set the capacityIncrement of this list
//  pgListInfo .................... print the information about the list
//==========================================================================*/

/*--------------------------------------------------------------------------
// Private functions
//------------------------------------------------------------------------*/
#define pgListCapacity(list)                  ((list)->capacity)
#define pgListCapacityIncrement(list)         ((list)->capacityIncrement)
#define pgListSetElementSize(list, val)       (list)->elementSize = val
#define pgListSetSize(list, val)              (list)->size = val
#define pgListSetData(list, val)              (list)->data = (void *)val
#define pgListSetCapacity(list, val)          (list)->capacity = val
#define pgListSetCapacityIncrement(list, val) (list)->capacityIncrement = val

#endif
