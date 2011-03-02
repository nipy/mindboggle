#include "pgutil.h"
#include "Heap.h"

#define INCREMENT 1000

Xheap xhInitEmpty()
{
   Xheap H;
   XheapElement he;

   H = (Xheap)pgList2(sizeof(XheapElement),0,INCREMENT);
   
   /* assume that a[0] = smallest floating point such as -1e33 
      serve as sentinel for stopping condition. 
      heap starts from 1 */      
   he.value = 1e-34;  /* not used in UpHeap anymore */
   he.id = -1;

   pgListAddElement(H, &he);

   return(H);
}

Xheap xhInit(XheapElement *array, int N)
{
   Xheap H;
   XheapElement *data, he;
   int i, *p;
   
   H = (Xheap)pgList2(sizeof(XheapElement),0,INCREMENT);

   he.value = 1e-34; 
   he.id = -1;

   pgListAddElement(H, &he);
   pgListSetSize(H,N+1);

   data = (XheapElement *)pgListData(H);
   for (i = 1; i <= N; i++) {
      data[i] = array[i-1];
      p = data[i].p;
      *p = i;
   }
   /* down build */
   for (i = N/2; i >= 1; i--)
      xhDownHeap(i, H); 

   return(H);
}

/* destroy the heap and free the memory */
void  xhDestroy(Xheap H)
{
   pgListDelete(H);
}


int xhUpHeap(int k, Xheap H)
{
   XheapElement *a, v;
   int k_father;
   int *p;

   a = (XheapElement *)pgListData(H);
   
   v = a[k];
   k_father = k/2;  /* integer divsion to retrieve its parent */
   while (k_father > 0 && a[k_father].value > v.value) {
      a[k] = a[k_father];
      p = a[k].p;
      *p = k;
      k = k_father;
      k_father = k/2;
   }
   a[k] = v;
   p = a[k].p;
   *p = k;

   return(k);
}

int xhDownHeap(int k, Xheap H)
{
   XheapElement *a, v;
   int N, k_minson;
   int *p;

   a = (XheapElement *)pgListData(H);
   N = xhSize(H);

   v = a[k];
   while ( k <= N/2 ) {
      k_minson = k+k;
      if ( k_minson < N ) {
	 if (a[k_minson].value > a[k_minson+1].value)
	    k_minson = k_minson + 1;  /* always locate the smallest son */
      }      
      if ( v.value <= a[k_minson].value ) 
	break; /* break out the loop */
      a[k] = a[k_minson]; 
      p = a[k].p;
      *p = k;
      k = k_minson;      /* go down one level */
   }
   a[k] = v;
   p = a[k].p;
   *p = k;

   return(k);
}


int xhInsert(double value, int id, int *p, Xheap H)
{
   XheapElement *a, v;
   int N, k;

   a = (XheapElement *)pgListData(H);

   v.value = value;
   v.id = id;
   v.p = p;

   pgListAddElement(H, &v);
   N = xhSize(H);
   k = xhUpHeap(N, H);

   return(k);
}

/* remove the smallest element */
XheapElement xhRemove(Xheap H)
{
   XheapElement v, *a;
   int N;
   
   N = xhSize(H);
   a = (XheapElement *)pgListData(H);

   v = a[1];
   a[1] = a[N];
   pgListSetSize(H, N);  
   /* the size of list is always 1 more than the size of heap 
      since the heap starts at 1 */

   xhDownHeap(1,H);

   return(v);
}

/* replace the smallest value with a new value if the new value is smaller
   otherwise the new value is returned and the heap is unchanged */
XheapElement xhReplace(double value, int id, int* p, Xheap H)
{
   XheapElement *a, v;

   a = (XheapElement *)pgListData(H);

   if ( value < a[1].value ) {
      v = a[1];
      a[1].value = value;
      a[1].id = id;
      a[1].p = p;
      xhDownHeap(1,H);
   } else {
      v.value = value;
      v.id = id;
      v.p = p;
   }
   return(v);
}

/* delete an item in the heap and its value is returned */
XheapElement xhDelete(int k, Xheap H)
{
   XheapElement *a, v;
   int N;

   N = xhSize(H);
   a = (XheapElement *)pgListData(H);

   v = a[k];

   a[k] = a[N];
   pgListSetSize(H, N);  
   /* the size of list is always 1 more than the size of heap 
      since the heap starts at 1 */

   xhDownHeap(k, H);

   return(v);
}

/* change the value of an item and its original value is returned */
XheapElement xhChange(int k, double value, int id, int *p, Xheap H)
{ 
   XheapElement *a, v;

   a = (XheapElement *)pgListData(H);
   v = a[k];
   
   if (value != a[k].value) {
      a[k].value = value;
      a[k].id = id;
      a[k].p = p;
      if ( value < v.value ) 
	 xhUpHeap(k, H);
      else 
	 xhDownHeap(k, H);
   }

   return(v);
}

/* change the value of an item and its original value is returned */
XheapElement xhChangeValue(int k, double value, Xheap H)
{ 
   XheapElement *a, v;

   a = (XheapElement *)pgListData(H);
   v = a[k];
   
   if (value != a[k].value) {
      a[k].value = value;
      if ( value < v.value ) 
	 xhUpHeap(k, H);
      else 
	 xhDownHeap(k, H);
   }

   return(v);
}

XheapElement xhGet(int k, Xheap H)
{
   XheapElement v;

   pgListElementAt(H, k, &v);

   return(v);
}

int xhSize(Xheap H)
{
   return( pgListSize(H) - 1 );
}

