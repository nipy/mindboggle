/*----------------------------------------------------------------------------
//
//      File: pgutil.c                          
//      A PG utility library                    
//                                              
//--------------------------------------------------------------------------*/

#include "pgutil.h"
#include <string.h>
#include <stdio.h>
#include <stdlib.h>


/*---------------------------------------------------------------------------
// Construct an empty PGlist with specified capacity and capacity increment 
//-------------------------------------------------------------------------*/
PGlist pgList2(int elementSize, int capacity, int capacityIncrement)
{
   PGlist list;
   void * data;

   if (elementSize<1) {
      fprintf(stderr, "pgList(): elementSize must be a postive integer!\n");
      exit(1);
   }
   if (capacity<0) {
      fprintf(stderr, "pgList(): capacity must not be negative!\n");
      exit(1);
   }
   
   //list = (PGlist)pgMalloc(sizeof(PGlistStruct));
   list = (PGlist)malloc(sizeof(PGlistStruct));
   
   pgListSetElementSize(list, elementSize);
   pgListSetSize(list, 0);
   pgListSetCapacity(list, capacity);
   pgListSetCapacityIncrement(list, capacityIncrement);
   if (capacity == 0) 
      pgListSetData(list, NULL);
   else {
      //data = (void *)pgMalloc(elementSize*capacity);
	  data = (void *)malloc(elementSize*capacity);
      pgListSetData(list, data);
   }

   return(list);
}

/*---------------------------------------------------------------------------
// Construct an empty PGlist with specified capacity and default 
// capacity increment as 100 
//-------------------------------------------------------------------------*/
PGlist pgList1(int elementSize, int capacity)
{
   PGlist list;

   list = pgList2(elementSize, capacity, 100);

   return(list);
}

/*---------------------------------------------------------------------------
// Construct an empty PGlist with default capacity as 0 and
// capacity increment as 100 
//-------------------------------------------------------------------------*/
PGlist pgList(int elementSize)
{
   PGlist list;

   list = pgList2(elementSize, 0, 100);

   return(list);
}


/*---------------------------------------------------------------------------
// Construct an empty PGlist with specified size, all the elements are set to 
// zero
//-------------------------------------------------------------------------*/
PGlist pgListOfSize(int size, int elementSize)
{
   PGlist list;
   char *data;
   int i;
   int capacity, capacityIncrement;

   if (size<0) {
      fprintf(stderr, "pgListOfSize(): size must not be negative!\n");
      exit(1);
   }

   capacity = size;
   capacityIncrement = 100;
   list = pgList2(elementSize, capacity, capacityIncrement);
   pgListSetSize(list, size);
   data = (char *)pgListData(list);
   for (i=0; i<elementSize*size; i++) 
      data[i] = 0;

   return(list);
}

/*---------------------------------------------------------------------------
// Delete this list
//-------------------------------------------------------------------------*/
void pgListDelete(PGlist list)
{
   void *data;
   
   data = pgListData(list);
   free(data);
   free(list);
}

/*---------------------------------------------------------------------------
// Add an element to this list
//-------------------------------------------------------------------------*/
void pgListAddElement(PGlist list, void *element)
{
   int size, capacity, elementSize, capacityIncrement;
   void *data;

   size        = pgListSize(list);
   capacity    = pgListCapacity(list);
   elementSize = pgListElementSize(list);
   data        = pgListData(list);
   if (size >= capacity) {
      capacityIncrement = pgListCapacityIncrement(list);
      capacity += capacityIncrement;
      pgListSetCapacity(list, capacity);
      if (data == NULL) {
	 /* initial list */
	 //data = (void *)pgMalloc(elementSize * capacity);
	 data = (void *)malloc(elementSize * capacity);
      } else {
	 /* allocate a larger list */
	 //data = (void *)pgRealloc(data, elementSize*capacity);
	 data = (void *)realloc(data, elementSize*capacity);
      }      
      pgListSetData(list, data);
   }

   memcpy((char *)data+size*elementSize, (char *)element, elementSize); 
   pgListSetSize(list, size+1);
}

/*---------------------------------------------------------------------------
// Add an integer to this list (must be a list consists of only integers)
//-------------------------------------------------------------------------*/
void pgListAddInt(PGlist list, int element)
{
   int size, capacity, elementSize, capacityIncrement;
   int *data;

   size        = pgListSize(list);
   capacity    = pgListCapacity(list);
   elementSize = pgListElementSize(list);
   data        = (int *)pgListData(list);
   
   if (size >= capacity) {
      capacityIncrement = pgListCapacityIncrement(list);
      capacity += capacityIncrement;
      pgListSetCapacity(list, capacity);
      if (data == NULL) {
	 /* initial list */
	 //data = (int *)pgMalloc(elementSize * capacity);
	 data = (int *)malloc(elementSize * capacity);
      } else {
	 /* allocate a larger list */
	 //data = (int *)pgRealloc(data, elementSize*capacity);
	 data = (int *)realloc(data, elementSize*capacity);
      }      
      pgListSetData(list, data);
   }

   data[size] = element;
   pgListSetSize(list, size+1);
}

/*---------------------------------------------------------------------------
// Add an array to this list
//-------------------------------------------------------------------------*/
void    pgListAddArray(PGlist list, void *array, int num)
{
   int size, capacity, elementSize, capacityIncrement, actualIncrement;
   void *data;

   size        = pgListSize(list);
   capacity    = pgListCapacity(list);
   elementSize = pgListElementSize(list);
   data        = pgListData(list);
   if (size + num > capacity) {   
      capacityIncrement = pgListCapacityIncrement(list);
      actualIncrement = (capacityIncrement > num)? capacityIncrement: num;
      capacity += actualIncrement;
      pgListSetCapacity(list, capacity);

      if (data == NULL) {
	 /* initial list */
	 //data = (void *)pgMalloc(elementSize * capacity);
	 data = (void *)malloc(elementSize * capacity);
      } else {
	 /* allocate a larger list */
	 //data = (void *)pgRealloc(data, elementSize*capacity);
	 data = (void *)realloc(data, elementSize*capacity);
      }
      pgListSetData(list, data);
   }

   memcpy((char *)data+size*elementSize, (char *)array, num*elementSize);
   pgListSetSize(list, size+num);      
}

/*---------------------------------------------------------------------------
// Insert an element into the list at the specified index
//-------------------------------------------------------------------------*/
PGerror pgListInsertElementAt(PGlist list, int index, void *element)
{
   int size, elementSize;
   void *data;
   void *tempPtr;
   char *currentPtr, *nextPtr;
   int i;

   size        = pgListSize(list);
   elementSize = pgListElementSize(list);
   
   if (index<0 || index>size-1) {
      return(PG_ERROR); /* out of bound error */
   }

   //tempPtr = (void *)pgMalloc(elementSize);
   tempPtr = (void *)malloc(elementSize);
   pgListAddElement(list, tempPtr);
   
   data        = pgListData(list);

   for (i=size-1; i>=index; i--) {
      currentPtr = (char *)data+i*elementSize;
      nextPtr    = (char *)currentPtr + elementSize;
      memcpy(nextPtr, currentPtr, elementSize);          
   }

   memcpy((char *)data+index*elementSize, (char *)element, elementSize);    

   return(PG_OK);
}



/*---------------------------------------------------------------------------
// Retrieve an element from this list at a given index
//-------------------------------------------------------------------------*/
PGerror pgListElementAt(PGlist list, int index, void *element)
{
   int size, elementSize;
   void *data;

   size        = pgListSize(list);
   elementSize = pgListElementSize(list);
   data        = pgListData(list);
   
   if (index<0 || index>size-1) {
      return(PG_ERROR); /* out of bound error */
   }
   memcpy((char *)element, (char *)data+index*elementSize, elementSize);    

   return(PG_OK);
}

/*---------------------------------------------------------------------------
// Sets a list element at a given index
//-------------------------------------------------------------------------*/
PGerror pgListSetElementAt(PGlist list, int index, void *element)
{
   int size, elementSize;
   void *data;

   size        = pgListSize(list);
   elementSize = pgListElementSize(list);
   data        = pgListData(list);
   
   if (index<0 || index>size-1) {
      return(PG_ERROR); /* out of bound error */
   }

   memcpy((char *)data+index*elementSize, (char *)element, elementSize);    

   return(PG_OK);
}

/*---------------------------------------------------------------------------
// Removes all elements from this list and sets its size to zero
//-------------------------------------------------------------------------*/
PGerror pgListRemoveElementAt(PGlist list, int index)
{
   int size, elementSize;
   void *data;
   char *currentPtr, *nextPtr;
   int i;

   size        = pgListSize(list);
   elementSize = pgListElementSize(list);
   data        = pgListData(list);
   
   if (index<0 || index>size-1) {
      return(PG_ERROR); /* out of bound error */
   }

   for (i=index; i<size-1; i++) {
      currentPtr = (char *)data+i*elementSize;
      nextPtr    = (char *)currentPtr + elementSize;
      memcpy(currentPtr, nextPtr, elementSize);          
   }

   pgListSetSize(list, size-1);

   return(PG_OK);   
}

/*---------------------------------------------------------------------------
// Removes all elements from this list and sets its size to zero
//-------------------------------------------------------------------------*/
void pgListRemoveAllElements(PGlist list)
{
   pgListSetSize(list, 0);
}

/*---------------------------------------------------------------------------
// Trim this list to current size
//-------------------------------------------------------------------------*/
void pgListTrim(PGlist list)
{
   void *data;
   int size, elementSize;
   
   size        = pgListSize(list);
   elementSize = pgListElementSize(list);
   data = pgListData(list);
   
   //data = (void *)pgRealloc(data, elementSize*size);
   data = (void *)realloc(data, elementSize*size);
   pgListSetData(list, data);
   pgListSetCapacity(list, size);
}

void pgListInfo(PGlist list)
{
   int elementSize, size, capacity, capacityIncrement;

   elementSize = pgListElementSize(list);
   size        = pgListSize(list);
   capacity    = pgListCapacity(list);
   capacityIncrement = pgListCapacityIncrement(list);

   printf("         elementSize = %d\n", elementSize);
   printf("                size = %d\n", size);
   printf("            capacity = %d\n", capacity);
   printf("   capacityIncrement = %d\n", capacityIncrement);
   printf("\n");
}


/*---------------------------------------------------------------------------
// pop out the top element from the stack
//-------------------------------------------------------------------------*/
void pgStackPop(PGstack stack, void *element)
{
   int size, elementSize;
   void *data;
   
   size        = pgListSize(stack);
   elementSize = pgListElementSize(stack);
   data        = pgListData(stack);
   
   memcpy((char *)element, (char *)data+(size-1)*elementSize, elementSize);    
   pgListSetSize(stack, size-1);
}

/*---------------------------------------------------------------------------
// Construct an empty PGqueue with specified capacity and capacity increment 
//-------------------------------------------------------------------------*/
PGqueue pgQueue2(int elementSize, int capacity, int capacityIncrement)
{
   PGqueue queue;
   
   //queue = (PGqueue)pgMalloc(sizeof(PGqueueStruct));
   queue = (PGqueue)malloc(sizeof(PGqueueStruct));

   queue->start = 0;
   queue->end   = -1;
   queue->list  = pgList2(elementSize, capacity, capacityIncrement);

   return(queue);
}

/*---------------------------------------------------------------------------
// Construct an empty PGqueue with specified capacity and default 
// capacity increment as 100 
//-------------------------------------------------------------------------*/
PGqueue pgQueue1(int elementSize, int capacity)
{
   PGqueue queue;
   
   queue = pgQueue2(elementSize, capacity, 100);

   return(queue);
}

/*---------------------------------------------------------------------------
// default constructor
//-------------------------------------------------------------------------*/
PGqueue pgQueue(int elementSize)
{
   PGqueue queue;
   
   queue = pgQueue2(elementSize, 0, 100);

   return(queue);
}

/* private functinos */

void    pgQueueEnsureSize(PGqueue queue)
{
   pgListSetSize(queue->list, pgQueueSize(queue));
}

/*---------------------------------------------------------------------------
// delete the queue and its memory
//-------------------------------------------------------------------------*/
void    pgQueueDelete(PGqueue queue)
{
   pgListDelete(queue->list);
   free(queue);
}

/*---------------------------------------------------------------------------
// removes all elements from the queue without releasing memory
//-------------------------------------------------------------------------*/
void    pgQueueRemoveAllElements(PGqueue queue) 
{
   queue->start = 0;
   queue->end   = -1;
   pgQueueEnsureSize(queue);
}

/* a private function for clean queue, ie. move things to the front */
void    pgQueueMoveToFront(PGqueue queue)
{
   void *data;
   void *s1, *s2;
   int elementSize, start, end;

   elementSize = pgListElementSize(queue->list);
   data  = pgListData(queue->list);
   start = queue->start;
   end   = queue->end;
   s2 = (char*)data + start*elementSize;
   s1 = data;
   memmove(s1, s2, (end - start + 1)*elementSize);
   queue->end   = end - start;
   queue->start = 0;
   pgQueueEnsureSize(queue);
}


/*---------------------------------------------------------------------------
// push an element to the end of the queue
//-------------------------------------------------------------------------*/
void    pgQueuePush(PGqueue queue, void* element)
{  
   int size, capacity;
   int q = PG_QUEUE_Q;   
                /* this is the factor determines the frequency of cleaning */
                /* a factor of n denotes that (n-1)*memory(queue) will be 
		   wasted */

   size     = pgQueueSize(queue);
   capacity = pgListCapacity(queue->list);

   if (queue->end >= (capacity - 1) && q*size < capacity ) 
      /* move the block to front, and release more memory */
      pgQueueMoveToFront(queue);

   /* just keep adding the element */
   queue->end = queue->end + 1;
   pgListAddElement(queue->list, element);
}

/*---------------------------------------------------------------------------
// push an array to the end of the queue
//-------------------------------------------------------------------------*/
void    pgQueuePushArray(PGqueue queue, void* array, int num)
{
   int size, capacity;
   int q = PG_QUEUE_Q;   
                /* this is the factor determines the frequency of cleaning */
                /* a factor of n denotes that (n-1)*memory(queue) will be 
		   wasted */

   size     = pgQueueSize(queue);
   capacity = pgListCapacity(queue->list);

   if (queue->end >= (capacity - 1) && q*size < capacity ) 
      /* move the block to front, and release more memory */
      pgQueueMoveToFront(queue);

   /* just keep adding the array */
   queue->end = queue->end + num;
   pgListAddArray(queue->list, array, num);   
}

/*---------------------------------------------------------------------------
// pop out the first element from the queue
//-------------------------------------------------------------------------*/
PGerror pgQueuePop(PGqueue queue, void* element)
{
   int start;

   start = queue->start;
   if (pgQueueIsEmpty(queue) == PG_FALSE) {
      pgListElementAt(queue->list, start, element); /* get the first one */   
      queue->start = start + 1;
      return(PG_TRUE); /* correctly got the result */
   } else {
      return(PG_FALSE); /* element is undefined */
   }
}

/*---------------------------------------------------------------------------
// trim the queue to its current size
//-------------------------------------------------------------------------*/
void    pgQueueTrim(PGqueue queue)
{
   pgQueueMoveToFront(queue);   

   pgListTrim(queue->list);
}

/*---------------------------------------------------------------------------
// extract an array from the queue
//-------------------------------------------------------------------------*/
void*   pgQueueToArray(PGqueue queue)
{
   void *arr, *data;
   int size, elementSize;

   size = pgQueueSize(queue);
   elementSize = pgListElementSize(queue->list);

   if (size == 0) {
      arr = NULL;
   } 
   else {
      data = pgListData(queue->list);
      //arr = (void *)pgMalloc(elementSize*size);
	  arr = (void *)malloc(elementSize*size);
      memcpy(arr, (char*)data + queue->start*elementSize, size*elementSize);
   }

   return(arr);
}


void    pgQueueInfo(PGqueue queue)
{
   printf("               start = %d\n", queue->start);
   printf("               end   = %d\n", queue->end);
   pgListInfo(queue->list);
}
