/*
 *    File:  quicksort.c
 * Purpose:  order a linear linked list of type TOTFUNC in
 *           decreasing order according to real values "xord"
 *           stored in each TOTVALUE element of the list
 *   Usage:  void quicksort(PTRTOTVALUE left, PTRTOTVALUE right)
 *
 * History:  2001 Aug 23  PB  adapted for linked lists from
 *                            quicksort.c in "A Book on C" (Kelley/Pohl)
 *
 */

#include <stdlib.h>
#include "bgtotal.h"

#define swap(x,y)   {float t; t = x; x = y; y = t;}
#define order(x,y)  if (x < y) swap (x,y)
#define o2(x,y)     order(x,y)
#define o3(x,y,z)   o2(x,y); o2(x,z); o2(y,z)

typedef enum {yes, no} yes_no;

static yes_no find_pivot(PTRTOTVALUE left, PTRTOTVALUE right,
                         float *pivot_ptr, int *n_ptr);
static PTRTOTVALUE partition(PTRTOTVALUE left, PTRTOTVALUE right,
                             float pivot, int n);


void quicksort(PTRTOTVALUE left, PTRTOTVALUE right)
{
   PTRTOTVALUE p;
   float pivot;
   int n;

   if (find_pivot(left, right, &pivot, &n) == yes) {
      p = partition(left, right, pivot, n);
      quicksort(left, p -> prev);
      quicksort(p, right);
   }
}

static yes_no find_pivot(PTRTOTVALUE left, PTRTOTVALUE right,
                         float *pivot_ptr, int *n_ptr)
{
   float a, b, c;
   PTRTOTVALUE p;
   int i, m;

   if (left == right) {
      *n_ptr = 1;
      return no;
   }

   a = left -> xord;

   m = 1;
   for (p = left; p != right; p = p -> next) {
      m++;
   }
   p = left;
   for (i = 0; i < m/2; i++) {
      p = p -> next;
   }
   *n_ptr = m;

   b = p -> xord;

   c = right -> xord;

   o3(a, b, c);

   if (a > b) {
      *pivot_ptr = b;
      return yes;
   }
   if (b > c) {
      *pivot_ptr = c;
      return yes;
   }
   for (p = left -> next; p != right; p = p -> next) {
      if (p -> xord != left -> xord) {
         *pivot_ptr = (p -> xord > left -> xord) ? left -> xord : p -> xord;
         return yes;
      }
   }
   return no;
}
      
static PTRTOTVALUE partition(PTRTOTVALUE left, PTRTOTVALUE right,
                             float pivot, int n)
{
   int i;
   float dnt,dat,SNt,xordt;
   PTRDKVALUE ptrpvt,ptrkp2t,ptrtdt,ptrkm2tt,ptrkm2bt;
   PTRBMVALUE ptrbm1t,ptrbm2t;

   i = 0;
   while (i < n) {
      while (left -> xord > pivot) {
         left = left -> next;
         i++;
      }
      while (right -> xord <= pivot) {
         right = right -> prev;
         i++;
      }
      if (i < n) {
         dnt = left -> dn;
         dat = left -> da;
         SNt = left -> SN;
         xordt = left -> xord;
         ptrpvt = left -> ptrpv;
         ptrkp2t = left -> ptrkp2;
         ptrtdt = left -> ptrtd;
         ptrkm2tt = left -> ptrkm2t;
         ptrkm2bt = left -> ptrkm2b;
         ptrbm1t = left -> ptrbm1;
         ptrbm2t = left -> ptrbm2;

         left -> dn = right -> dn;
         left -> da = right -> da;
         left -> SN = right -> SN;
         left -> xord = right -> xord;
         left -> ptrpv = right -> ptrpv;
         left -> ptrkp2 = right -> ptrkp2;
         left -> ptrtd = right -> ptrtd;
         left -> ptrkm2t = right -> ptrkm2t;
         left -> ptrkm2b = right -> ptrkm2b;
         left -> ptrbm1 = right -> ptrbm1;
         left -> ptrbm2 = right -> ptrbm2;

         right -> dn = dnt;
         right -> da = dat;
         right -> SN = SNt;
         right -> xord = xordt;
         right -> ptrpv = ptrpvt;
         right -> ptrkp2 = ptrkp2t;
         right -> ptrtd = ptrtdt;
         right -> ptrkm2t = ptrkm2tt;
         right -> ptrkm2b = ptrkm2bt;
         right -> ptrbm1 = ptrbm1t;
         right -> ptrbm2 = ptrbm2t;

         left = left -> next;
         i++;
         right = right -> prev;
         i++;
      }
   }
   return left;
}

