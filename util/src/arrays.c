/*----------------------------------------------------------------------
  File    : arrays.c
  Contents: some basic array operations, especially for pointer arrays
  Author  : Christian Borgelt
  History : 1996.09.16 file created as arrays.c
            1999.02.04 long int changed to int (later to size_t)
            2001.06.03 function ptr_shuffle() added
            2002.01.02 functions for basic data types added
            2002.03.03 functions ptr_reverse() etc. added
            2003.08.21 function ptr_heapsort() added
            2007.01.16 shuffle functions for basic data types added
            2007.12.02 bug in reverse functions fixed
            2008.08.01 renamed to arrays.c, some functions added
            2008.08.11 main function added (sortargs)
            2008.08.12 functions ptr_unique() etc. added
            2008.08.17 binary search functions improved
            2008.10.05 functions to clear arrays added
            2010.07.31 index array sorting functions added
            2010.12.07 added several explicit type casts
            2011.09.28 function ptr_mrgsort() added (merge sort)
            2011.09.30 merge sort combined with insertion sort
            2012.06.03 functions for data type long int added
            2013.03.07 direction parameter added to sorting functions
            2013.03.10 binary search and bisection functions separated
            2013.03.20 adapted return values and arguments to ptrdiff_t
            2013.03.27 index sorting for types int, long and ptrdiff_t
            2013.07.24 bug in move functions fixed (forward move)
----------------------------------------------------------------------*/
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include "arrays.h"

/*----------------------------------------------------------------------
  Preprocessor Definitions
----------------------------------------------------------------------*/
#define BUFSIZE     1024        /* size of fixed buffer for moving */
#define TH_INSERT   16          /* threshold for insertion sort */

/*----------------------------------------------------------------------
  Functions for Arrays of Basic Data Types
----------------------------------------------------------------------*/

#define MOVE(name,type) \
void name##_move (type *array, size_t off, size_t n, size_t pos)       \
{                               /* --- move a number array section */  \
  size_t i, end;                /* loop variable, end index */         \
  type   *src, *dst;            /* to traverse the array */            \
  type   fxd[BUFSIZE], *buf;    /* buffer for copying */               \
                                                                       \
  assert(array);                /* check the function arguments */     \
  if ((pos >= off) && (pos < off +n))                                  \
    return;                     /* check whether moving is necessary */\
  if (pos < off) { end = off +n; off = pos; pos = end -n; }            \
  else           { end = pos +1;            pos = off +n; }            \
  buf = fxd;                    /* normalize the indices */            \
  if (pos +pos > end +off) {    /* if first section is smaller */      \
/*if (pos -off < end -pos) { */                                        \
    while (pos > off) {         /* while there are elements to shift */\
      n = pos -off;             /* get the number of elements */       \
      if (n > BUFSIZE) {        /* if the fixed buffer is too small */ \
        buf = (type*)malloc(n *sizeof(type));                          \
        if (!buf) { buf = fxd; n = BUFSIZE; }                          \
      }                         /* try to allocate a fitting buffer */ \
      src = array + pos -n;     /* get number of elements and */       \
      dst = buf;                /* copy source to the buffer */        \
      for (i = n;        i > 0; i--) *dst++ = *src++;                  \
      dst = array +pos -n;      /* shift down/left second section */   \
      for (i = end -pos; i > 0; i--) *dst++ = *src++;                  \
      src = buf;                /* copy buffer to destination */       \
      for (i = n;        i > 0; i--) *dst++ = *src++;                  \
      pos -= n; end -= n;       /* second section has been shifted */  \
    } }                         /* down/left cnt elements */           \
  else {                        /* if second section is smaller */     \
    while (end > pos) {         /* while there are elements to shift */\
      n = end -pos;             /* get the number of elements */       \
      if (n > BUFSIZE) {        /* if the fixed buffer is too small */ \
        buf = (type*)malloc(n *sizeof(type));                          \
        if (!buf) { buf = fxd; n = BUFSIZE; }                          \
      }                         /* try to allocate a fitting buffer */ \
      src = array +pos +n;      /* get number of elements and */       \
      dst = buf +n;             /* copy source to the buffer */        \
      for (i = n;        i > 0; i--) *--dst = *--src;                  \
      dst = array +pos +n;      /* shift up/right first section */     \
      for (i = pos -off; i > 0; i--) *--dst = *--src;                  \
      src = buf +n;             /* copy buffer to destination */       \
      for (i = n;        i > 0; i--) *--dst = *--src;                  \
      pos += n; off += n;       /* first section has been shifted */   \
    }                           /* up/right cnt elements */            \
  }                                                                    \
  if (buf != fxd) free(buf);    /* delete an allocated buffer */       \
}  /* move() */

/*--------------------------------------------------------------------*/

MOVE(sht, short)
MOVE(int, int)
MOVE(lng, long)
MOVE(dif, diff_t)
MOVE(siz, size_t)
MOVE(flt, float)
MOVE(dbl, double)

/*--------------------------------------------------------------------*/

#define SELECT(name,type) \
void name##_select (type *array, size_t n, size_t k, RANDFN *rand)     \
{                               /* --- shuffle array entries */        \
  size_t i;                     /* array index */                      \
  type   t;                     /* exchange buffer */                  \
                                                                       \
  assert(array && (n >= k));    /* check the function arguments */     \
  k = (k < n) ? k+1 : n;        /* adapt the number of selections */   \
  while (--k > 0) {             /* shuffle loop (k selections) */      \
    i = (size_t)(rand() *(double)n);     /* compute a random index */  \
    if (i > --n) i = n;         /* and exchange the array elements */  \
    t = array[i]; array[i] = *array; *array++ = t;                     \
  }                                                                    \
}  /* select() */

/*--------------------------------------------------------------------*/

SELECT(sht, short)
SELECT(int, int)
SELECT(lng, long)
SELECT(dif, diff_t)
SELECT(siz, size_t)
SELECT(flt, float)
SELECT(dbl, double)

/*--------------------------------------------------------------------*/

#define SHUFFLE(name,type) \
void name##_shuffle (type *array, size_t n, RANDFN *rand) \
{ name##_select(array, n, n-1, rand); }

/*--------------------------------------------------------------------*/

SHUFFLE(sht, short)
SHUFFLE(int, int)
SHUFFLE(lng, long)
SHUFFLE(dif, diff_t)
SHUFFLE(siz, size_t)
SHUFFLE(flt, float)
SHUFFLE(dbl, double)

/*--------------------------------------------------------------------*/

#define REVERSE(name, type) \
void name##_reverse (type *array, size_t n)                            \
{                               /* --- reverse a number array */       \
  type *end = array +n;         /* end of the array to reverse */      \
  type t;                       /* exchange buffer */                  \
                                                                       \
  while (--end > array) {       /* reverse the order of the elems. */  \
    t = *end; *end = *array; *array++ = t; }                           \
}  /* reverse */

/*--------------------------------------------------------------------*/

REVERSE(sht, short)
REVERSE(int, int)
REVERSE(lng, long)
REVERSE(dif, diff_t)
REVERSE(siz, size_t)
REVERSE(flt, float)
REVERSE(dbl, double)

/*--------------------------------------------------------------------*/

#define QSORT(name,type) \
static void name##_qrec (type *array, size_t n)                        \
{                               /* --- recursive part of sort */       \
  type   *l, *r;                /* pointers to exchange positions */   \
  type   x, t;                  /* pivot element and exchange buffer */\
  size_t m;                     /* number of elements in sections */   \
                                                                       \
  do {                          /* sections sort loop */               \
    l = array; r = l +n -1;     /* start at left and right boundary */ \
    if (*l > *r) { t = *l; *l = *r; *r = t; }                          \
    x = array[n/2];             /* get the middle element as pivot */  \
    if      (x < *l) x = *l;    /* compute median of three */          \
    else if (x > *r) x = *r;    /* to find a better pivot */           \
    while (1) {                 /* split and exchange loop */          \
      while (*++l < x);         /* skip smaller elems. on the left */  \
      while (*--r > x);         /* skip greater elems. on the right */ \
      if (l >= r) {             /* if less than two elements left, */  \
        if (l <= r) { l++; r--; } break; }       /* abort the loop */  \
      t = *l; *l = *r; *r = t;  /* otherwise exchange elements */      \
    }                                                                  \
    m = (size_t)(array +n -l);  /* compute the number of elements */   \
    n = (size_t)(r -array +1);  /* right and left of the split */      \
    if (n > m) {                /* if right section is smaller, */     \
      if (m >= TH_INSERT)       /* but larger than the threshold, */   \
        name##_qrec(l, m); }    /* sort it by an recursive call */     \
    else {                      /* if the left section is smaller, */  \
      if (n >= TH_INSERT)       /* but larger than the threshold, */   \
        name##_qrec(array, n);  /* sort it by an recursive call, */    \
      array = l; n = m;         /* then switch to the right section */ \
    }                           /* keeping its size m in variable n */ \
  } while (n >= TH_INSERT);     /* while greater than threshold */     \
}  /* qrec() */                                                        \
                                                                       \
/*------------------------------------------------------------------*/ \
                                                                       \
void name##_qsort (type *array, size_t n, int dir)                     \
{                               /* --- sort a number array */          \
  size_t i, k;                  /* loop variable, first section */     \
  type   *l, *r;                /* to traverse the array */            \
  type   t;                     /* exchange buffer */                  \
                                                                       \
  assert(array);                /* check the function arguments */     \
  if (n < 2) return;            /* do not sort less than two elems. */ \
  if (n < TH_INSERT)            /* if less elements than threshold */  \
    k = n;                      /* for insertion sort, note the */     \
  else {                        /* number of elements, otherwise */    \
    name##_qrec(array, n);      /* call the recursive sort function */ \
    k = TH_INSERT -1;           /* and get the number of elements */   \
  }                             /* in the first array section */       \
  for (l = r = array; --k > 0;) /* find position of smallest element */\
    if (*++r < *l) l = r;       /* within the first k elements */      \
  r = array;                    /* swap the smallest element */        \
  t = *l; *l = *r; *r = t;      /* to the front as a sentinel */       \
  for (i = n; --i > 0; ) {      /* standard insertion sort */          \
    t = *++r;                   /* note the number to insert */        \
    for (l = r; *--l > t; )     /* shift right all numbers that are */ \
      l[1] = *l;                /* greater than the one to insert */   \
    l[1] = t;                   /* and store the number to insert */   \
  }                             /* in the place thus found */          \
  if (dir < 0)                  /* if descending order requested, */   \
    name##_reverse(array, n);   /* reverse the element order */        \
}  /* qsort() */

/*--------------------------------------------------------------------*/

QSORT(sht, short)
QSORT(int, int)
QSORT(lng, long)
QSORT(dif, diff_t)
QSORT(siz, size_t)
QSORT(flt, float)
QSORT(dbl, double)

/*--------------------------------------------------------------------*/

#define HEAPSORT(name,type) \
static void name##_sift (type *array, size_t l, size_t r)              \
{                               /* --- let element sift down in heap */\
  size_t i;                     /* index of first successor in heap */ \
  type   t;                     /* buffer for an array element */      \
                                                                       \
  t = array[l];                 /* note the sift element */            \
  i = l +l +1;                  /* compute index of first successor */ \
  do {                          /* sift loop */                        \
    if ((i < r) && (array[i] < array[i+1]))                            \
      i++;                      /* if second successor is greater */   \
    if (t >= array[i])          /* if the successor is greater */      \
      break;                    /* than the sift element, */           \
    array[l] = array[i];        /* let the successor ascend in heap */ \
    l = i; i += i +1;           /* compute index of first successor */ \
  } while (i <= r);             /* while still within heap */          \
  array[l] = t;                 /* store the sift element */           \
}  /* sift() */                                                        \
                                                                       \
/*------------------------------------------------------------------*/ \
                                                                       \
void name##_heapsort (type *array, size_t n, int dir)                  \
{                               /* --- heap sort for number arrays */  \
  size_t l, r;                  /* boundaries of heap section */       \
  type   t;                     /* exchange buffer */                  \
                                                                       \
  assert(array);                /* check the function arguments */     \
  if (n < 2) return;            /* do not sort less than two elems. */ \
  l = n /2;                     /* at start, only the second half */   \
  r = n -1;                     /* of the array has heap structure */  \
  while (l > 0)                 /* while the heap is not complete, */  \
    name##_sift(array, --l, r); /* extend it by one element */         \
  while (1) {                   /* heap reduction loop */              \
    t = array[0];               /* swap the greatest element */        \
    array[0] = array[r];        /* to the end of the array */          \
    array[r] = t;                                                      \
    if (--r <= 0) break;        /* if the heap is empty, abort */      \
    name##_sift(array, 0, r);   /* let swapped element sift down */    \
  }                                                                    \
  if (dir < 0)                  /* if descending order requested, */   \
    name##_reverse(array, n);   /* reverse the element order */        \
}  /* heapsort() */

/*--------------------------------------------------------------------*/

HEAPSORT(sht, short)
HEAPSORT(int, int)
HEAPSORT(lng, long)
HEAPSORT(dif, diff_t)
HEAPSORT(siz, size_t)
HEAPSORT(flt, float)
HEAPSORT(dbl, double)

/*--------------------------------------------------------------------*/

#define UNIQUE(name,type) \
size_t name##_unique (type *array, size_t n)                           \
{                               /* --- remove duplicate elements */    \
  type *s, *d;                  /* to traverse the array */            \
                                                                       \
  assert(array);                /* check the function arguments */     \
  if (n <= 1) return n;         /* check for 0 or 1 element */         \
  for (d = s = array; --n > 0;) /* traverse the (sorted) array and */  \
    if (*++s != *d) *++d = *s;  /* collect the unique elements */      \
  return (size_t)(++d -array);  /* return new number of elements */    \
}  /* unique() */

/*--------------------------------------------------------------------*/

UNIQUE(sht, short)
UNIQUE(int, int)
UNIQUE(lng, long)
UNIQUE(dif, diff_t)
UNIQUE(siz, size_t)
UNIQUE(flt, float)
UNIQUE(dbl, double)

/*--------------------------------------------------------------------*/

#define BSEARCH(name,type) \
diff_t name##_bsearch (type key, const type *array, size_t n)          \
{                               /* --- do a binary search */           \
  size_t l, r, m;               /* array indices */                    \
  type   t;                     /* array element */                    \
                                                                       \
  assert(array);                /* check the function arguments */     \
  for (l = 0, r = n; l < r; ) { /* while search range is not empty */  \
    t = array[m = (l+r)/2];     /* compare the given key */            \
    if      (key > t) l = m+1;  /* to the middle element and */        \
    else if (key < t) r = m;    /* adapt the search range */           \
    else return (diff_t)m;      /* according to the result */          \
  }                             /* if match found, return index */     \
  return (diff_t)-1;            /* return 'not found' */               \
}  /* bsearch() */

/*--------------------------------------------------------------------*/

BSEARCH(sht, short)
BSEARCH(int, int)
BSEARCH(lng, long)
BSEARCH(dif, diff_t)
BSEARCH(siz, size_t)
BSEARCH(flt, float)
BSEARCH(dbl, double)

/*--------------------------------------------------------------------*/

#define BISECT(name,type) \
size_t name##_bisect (type key, const type *array, size_t n)           \
{                               /* --- do a bisection search */        \
  size_t l, r, m;               /* array indices */                    \
  type   t;                     /* array element */                    \
                                                                       \
  assert(array);                /* check the function arguments */     \
  for (l = 0, r = n; l < r; ) { /* while search range is not empty */  \
    t = array[m = (l+r) /2];    /* compare the given key */            \
    if      (key > t) l = m+1;  /* to the middle element and */        \
    else if (key < t) r = m;    /* adapt the search range */           \
    else return m;              /* according to the result */          \
  }                             /* if match found, return index */     \
  return l;                     /* return the insertion position */    \
}  /* bisect() */

/*--------------------------------------------------------------------*/

BISECT(sht, short)
BISECT(int, int)
BISECT(lng, long)
BISECT(dif, diff_t)
BISECT(siz, size_t)
BISECT(flt, float)
BISECT(dbl, double)

/*----------------------------------------------------------------------
  Functions for Pointer Arrays
----------------------------------------------------------------------*/

void ptr_move (void *array, size_t off, size_t n, size_t pos)
{                               /* --- move a pointer array section */
  size_t i, end;                /* loop variable, end index */
  void   **src, **dst;          /* to traverse the array */
  void   *fxd[BUFSIZE], **buf;  /* buffer for copying */

  assert(array);                /* check the function arguments */
  if ((pos >= off) && (pos < off +n))
    return;                     /* check whether moving is necessary */
  if (pos < off) { end = off +n; off = pos; pos = end -n; }
  else           { end = pos +1;            pos = off +n; }
  buf = fxd;                    /* normalize the indices */
  if (pos +pos > end +off) {    /* if first section is smaller */
/*if (pos -off < end -pos) { */
    while (pos > off) {         /* while there are elements to shift */
      n = pos -off;             /* get the number of elements */
      if (n > BUFSIZE) {        /* if the fixed buffer is too small */
        buf = (void**)malloc(n *sizeof(void*));
        if (!buf) { buf = fxd; n = BUFSIZE; }
      }                         /* try to allocate a fitting buffer */
      src = (void**)array+pos-n;/* get number of elements and */
      dst = buf;                /* copy source to the buffer */
      for (i = n;        i > 0; i--) *dst++ = *src++;
      dst = (void**)array+pos-n;/* shift down/left second section */
      for (i = end -pos; i > 0; i--) *dst++ = *src++;
      src = buf;                /* copy buffer to destination */
      for (i = n;        i > 0; i--) *dst++ = *src++;
      pos -= n; end -= n;       /* second section has been shifted */
    } }                         /* down/left cnt elements */
  else {                        /* if second section is smaller */
    while (end > pos) {         /* while there are elements to shift */
      n = end -pos;             /* get the number of elements */
      if (n > BUFSIZE) {        /* if the fixed buffer is too small */
        buf = (void**)malloc(n *sizeof(void*));
        if (!buf) { buf = fxd; n = BUFSIZE; }
      }                         /* try to allocate a fitting buffer */
      src = (void**)array+pos+n;/* get number of elements and */
      dst = buf +n;             /* copy source to the buffer */
      for (i = n;        i > 0; i--) *--dst = *--src;
      dst = (void**)array+pos+n;/* shift up/right first section */
      for (i = pos -off; i > 0; i--) *--dst = *--src;
      src = buf +n;             /* copy buffer to destination */
      for (i = n;        i > 0; i--) *--dst = *--src;
      pos += n; off += n;       /* first section has been shifted */
    }                           /* up/right cnt elements */
  }
  if (buf != fxd) free(buf);    /* delete an allocated buffer */
}  /* ptr_move() */

/*--------------------------------------------------------------------*/

void ptr_select (void *array, size_t n, size_t k, RANDFN *rand)
{                               /* --- select random array entries */
  size_t i;                     /* array index */
  void   **a = (void**)array;   /* array to sort */
  void   *t;                    /* exchange buffer */

  assert(array && rand && (n >= k));  /* check the function arguments */
  k = (k < n) ? k+1 : n;        /* adapt the number of selections */
  while (--k > 0) {             /* shuffle loop (k selections) */
    i = (size_t)(rand() *(double)n); /* compute a random index */
    if (i > --n) i = n;         /* and clamp it to a valid range */
    t = a[i]; a[i] = *a; *a++ = t;
  }                             /* exchange the array elements */
}  /* ptr_select() */

/*--------------------------------------------------------------------*/

void ptr_shuffle (void *array, size_t n, RANDFN *rand)
{ ptr_select(array, n, n-1, rand); }

/*--------------------------------------------------------------------*/

void ptr_reverse (void *array, size_t n)
{                               /* --- reverse a pointer array */
  void **a = (void**)array;     /* array to reverse */
  void **e = a +n;              /* end of array to reverse */
  void *t;                      /* exchange buffer */

  assert(array);                /* check the function arguments */
  while (--e > a) {             /* reverse the order of the elements */
    t = *e; *e = *a; *a++ = t; }
}  /* ptr_reverse() */

/*--------------------------------------------------------------------*/

static void qrec (void **array, size_t n, CMPFN *cmp, void *data)
{                               /* --- recursive part of quicksort */
  void   **l, **r;              /* pointers to exchange positions */
  void   *x,  *t;               /* pivot element and exchange buffer */
  size_t m;                     /* number of elements in 2nd section */

  do {                          /* sections sort loop */
    l = array; r = l +n -1;     /* start at left and right boundary */
    if (cmp(*l, *r, data) > 0){ /* bring the first and last */
      t = *l; *l = *r; *r = t;} /* element into proper order */
    x = array[n /2];            /* get the middle element as pivot */
    if      (cmp(x, *l, data) < 0) x = *l;  /* try to find a */
    else if (cmp(x, *r, data) > 0) x = *r;  /* better pivot */
    while (1) {                 /* split and exchange loop */
      while (cmp(*++l, x, data) < 0)      /* skip left  elements that */
        ;                       /* are smaller than the pivot element */
      while (cmp(*--r, x, data) > 0)      /* skip right elements that */
        ;                       /* are greater than the pivot element */
      if (l >= r) {             /* if less than two elements left, */
        if (l <= r) { l++; r--; } break; }       /* abort the loop */
      t = *l; *l = *r; *r = t;  /* otherwise exchange elements */
    }
    m = (size_t)(array +n -l);  /* compute the number of elements */
    n = (size_t)(r -array +1);  /* right and left of the split */
    if (n > m) {                /* if right section is smaller, */
      if (m >= TH_INSERT)       /* but larger than the threshold, */
        qrec(l, m, cmp, data); }/* sort it by a recursive call, */
    else {                      /* if the left section is smaller, */
      if (n >= TH_INSERT)       /* but larger than the threshold, */
        qrec(array,n,cmp,data); /* sort it by a recursive call, */
      array = l; n = m;         /* then switch to the right section */
    }                           /* keeping its size m in variable n */
  } while (n >= TH_INSERT);     /* while greater than threshold */
}  /* qrec() */

/*--------------------------------------------------------------------*/

void ptr_qsort (void *array, size_t n, int dir, CMPFN *cmp, void *data)
{                               /* --- quicksort for pointer arrays */
  size_t i, k;                  /* loop variable, first section */
  void   **l, **r;              /* to traverse the array */
  void   *t;                    /* exchange buffer */

  assert(array && cmp);         /* check the function arguments */
  if (n <= 1) return;           /* do not sort less than two elements */
  if (n < TH_INSERT)            /* if fewer elements than threshold */
    k = n;                      /* for insertion sort, note the */
  else {                        /* number of elements, otherwise */
    qrec((void**)array, n, cmp, data);
    k = TH_INSERT -1;           /* call the recursive function and */
  }                             /* get size of first array section */
  for (l = r = (void**)array; --k > 0; )
    if (cmp(*++r, *l, data) < 0)
      l = r;                    /* find smallest of first k elements */
  r = (void**)array;            /* swap the smallest element */
  t = *l; *l = *r; *r = t;      /* to the front as a sentinel */
  for (i = n; --i > 0; ) {      /* standard insertion sort */
    t = *++r;                   /* note the element to insert */
    for (l = r; cmp(*--l, t, data) > 0; ) /* shift right elements */
      l[1] = *l;                /* that are greater than the one */
    l[1] = t;                   /* to insert and store this element */
  }                             /* in the place thus found */
  if (dir < 0)                  /* if descending order requested, */
    ptr_reverse(array, n);      /* reverse the element order */
}  /* ptr_qsort() */

/*--------------------------------------------------------------------*/

static void sift (void **array, size_t l, size_t r,
                  CMPFN *cmp, void *data)
{                               /* --- let element sift down in heap */
  size_t i;                     /* index of first successor in heap */
  void   *t;                    /* buffer for an array element */

  t = array[l];                 /* note the sift element */
  i = l +l +1;                  /* compute index of first successor */
  do {                          /* sift loop */
    if ((i < r)                 /* if second successor is greater */
    &&  (cmp(array[i], array[i+1], data) < 0))
      i++;                      /* go to the second successor */
    if (cmp(t, array[i], data) >= 0) /* if the successor is greater */
      break;                         /* than the sift element, */
    array[l] = array[i];        /* let the successor ascend in heap */
    l = i; i += i +1;           /* compute index of first successor */
  } while (i <= r);             /* while still within heap */
  array[l] = t;                 /* store the sift element */
}  /* sift() */

/*--------------------------------------------------------------------*/

void ptr_heapsort (void *array, size_t n, int dir,
                   CMPFN *cmp, void *data)
{                               /* --- heap sort for pointer arrays */
  size_t l, r;                  /* boundaries of heap section */
  void   *t, **a;               /* exchange buffer, array */

  assert(array && cmp);         /* check the function arguments */
  if (n <= 1) return;           /* do not sort less than two elements */
  l = n /2;                     /* at start, only the second half */
  r = n -1;                     /* of the array has heap structure */
  a = (void**)array;            /* type the array pointer */
  while (l > 0)                 /* while the heap is not complete, */
    sift(a, --l, r, cmp, data); /* extend it by one element */
  while (1) {                   /* heap reduction loop */
    t = a[0]; a[0] = a[r];      /* swap the greatest element */
    a[r] = t;                   /* to the end of the array */
    if (--r <= 0) break;        /* if the heap is empty, abort */
    sift(a, 0, r, cmp, data);   /* let the element that has been */
  }                             /* swapped to front sift down */
  if (dir < 0)                  /* if descending order requested, */
    ptr_reverse(array, n);      /* reverse the element order */
}  /* ptr_heapsort() */

/*--------------------------------------------------------------------*/

static void mrgsort (void **array, void **buf, size_t n,
                     CMPFN *cmp, void *data)
{                               /* --- merge sort for pointer arrays */
  size_t k, a, b;               /* numbers of objects in sections */
  void   **sa, **sb, **ea, **eb;/* starts and ends of sorted sections */
  void   **d, *t;               /* merge destination, exchange buffer */

  assert(array && buf && cmp);  /* check the function arguments */
  if (n <= 8) {                 /* if only few elements to sort */
    for (sa = array; --n > 0;){ /* insertion sort loop */
      t = *(d = ++sa);          /* note the element to insert */
      while ((--d >= array)     /* while not at the array start, */
      &&     (cmp(*d, t, data) > 0))     /* shift right elements */
        d[1] = *d;              /* that are greater than the one */
      d[1] = t;                 /* to insert and store the element */
    } return;                   /* to insert in the place thus found */
  }                             /* aftwards sorting is done, so abort */
  /* Using insertion sort for less than eight elements is not only */
  /* slightly faster, but also ensures that all subsections sorted */
  /* recursively in the code below contain at least two elements.  */

  k = n/2; d = buf;             /* sort two subsections recursively */
  mrgsort(sa = array,   d,   a = k/2, cmp, data);
  mrgsort(sb = sa+a,    d+a, b = k-a, cmp, data);
  for (ea = sb, eb = sb+b; 1;){ /* traverse the sorted sections */
    if (cmp(*sa, *sb, data) <= 0)
         { *d++ = *sa++; if (sa >= ea) break; }
    else { *d++ = *sb++; if (sb >= eb) break; }
  }                             /* copy smaller element to dest. */
  while (sa < ea) *d++ = *sa++; /* copy remaining elements */
  while (sb < eb) *d++ = *sb++; /* from source to destination */

  n -= k; d = buf+k;            /* sort two subsections recursively */
  mrgsort(sa = array+k, d,   a = n/2, cmp, data);
  mrgsort(sb = sa+a,    d+a, b = n-a, cmp, data);
  for (ea = sb, eb = sb+b; 1;){ /* traverse the sorted sections */
    if (cmp(*sa, *sb, data) <= 0)
         { *d++ = *sa++; if (sa >= ea) break; }
    else { *d++ = *sb++; if (sb >= eb) break; }
  }                             /* copy smaller element to dest. */
  while (sa < ea) *d++ = *sa++; /* copy remaining elements */
  while (sb < eb) *d++ = *sb++; /* from source to destination */

  sa = buf; sb = sa+k; d = array;
  for (ea = sb, eb = sb+n; 1;){ /* traverse the sorted sections */
    if (cmp(*sa, *sb, data) <= 0)
         { *d++ = *sa++; if (sa >= ea) break; }
    else { *d++ = *sb++; if (sb >= eb) break; }
  }                             /* copy smaller element to dest. */
  while (sa < ea) *d++ = *sa++; /* copy remaining elements */
  while (sb < eb) *d++ = *sb++; /* from source to destination */
}  /* mrgsort() */

/*--------------------------------------------------------------------*/

int ptr_mrgsort (void *array, size_t n, int dir,
                 CMPFN *cmp, void *data, void *buf)
{                               /* --- merge sort for pointer arrays */
  void **b;                     /* (allocated) buffer */

  assert(array && cmp);         /* check the function arguments */
  if (n < 2) return 0;          /* do not sort less than two objects */
  if (!(b = (void**)buf) && !(b = (void**)malloc(n *sizeof(void*))))
    return -1;                  /* allocate a buffer if not given */
  mrgsort(array, buf, n, cmp, data);
  if (!buf) free(b);            /* sort the array recursively */
  if (dir < 0)                  /* if descending order requested, */
    ptr_reverse(array, n);      /* reverse the element order */
  return 0;                     /* return 'ok' */
}  /* ptr_mrgsort() */

/* This implementation of merge sort is stable, that is, it does not  */
/* change the relative order of elements that are considered equal by */
/* the comparison function. Thus it maintains the order of a previous */
/* sort with another comparison function as long as the order imposed */
/* by the current comparison function does not override this order.   */

/*--------------------------------------------------------------------*/

size_t ptr_unique (void *array, size_t n,
                   CMPFN *cmp, void *data, OBJFN *del)
{                               /* --- remove duplicate elements */
  void **s, **d;                /* to traverse the pointer array */

  assert(array && cmp);         /* check the function arguments */
  if (n <= 1) return n;         /* check for 0 or 1 element */
  for (d = s = (void**)array; --n > 0; ) {
    if (cmp(*++s, *d, data) != 0) *++d = *s;
    else if (del) del(*s);      /* traverse the (sorted) array */
  }                             /* and collect unique elements */
  return (size_t)(++d -(void**)array);
}  /* ptr_unique() */           /* return the new number of elements */

/*--------------------------------------------------------------------*/

diff_t ptr_bsearch (const void *key, const void *array, size_t n,
                    CMPFN *cmp, void *data)
{                               /* --- do a binary search */
  size_t l, r, m;               /* array indices */
  int    c;                     /* comparison result */
  void   **a;                   /* typed array */

  assert(key && array && cmp);  /* check the function arguments */
  a = (void**)array;            /* type the array pointer */
  for (l = 0, r = n; l < r; ) { /* while search range is not empty */
    m = (l+r)/2;                /* compare the given key */
    c = cmp(key, a[m], data);   /* to the middle element */
    if      (c > 0) l = m+1;    /* adapt the search range */
    else if (c < 0) r = m;      /* according to the result */
    else return (diff_t)m;      /* if match found, return index */
  }
  return (diff_t)-1;            /* return 'not found' */
}  /* ptr_bsearch() */

/*--------------------------------------------------------------------*/

size_t ptr_bisect (const void *key, const void *array, size_t n,
                   CMPFN *cmp, void *data)
{                               /* --- do a binary search */
  size_t l, r, m;               /* array indices */
  int    c;                     /* comparison result */
  void   **a;                   /* typed array */

  assert(key && array && cmp);  /* check the function arguments */
  a = (void**)array;            /* type the array pointer */
  for (l = 0, r = n; l < r; ) { /* while search range is not empty */
    m = (l+r)/2;                /* compare the given key */
    c = cmp(key, a[m], data);   /* to the middle element */
    if      (c > 0) l = m+1;    /* adapt the search range */
    else if (c < 0) r = m;      /* according to the result */
    else return m;              /* if match found, return index */
  }
  return l;                     /* return the insertion position */
}  /* ptr_bisect() */

/*----------------------------------------------------------------------
  Functions for Index Arrays
----------------------------------------------------------------------*/

#define IDX_QSORT(name,tin,tidx,type)                                  \
static void name##_qrec (tidx *index, size_t n, const type *array)     \
{                               /* --- recursive part of sort */       \
  tidx   *l, *r;                /* pointers to exchange positions */   \
  tidx   x, t;                  /* pivot element and exchange buffer */\
  size_t m;                     /* number of elements in sections */   \
  type   p, a, z;               /* buffers for array elements */       \
                                                                       \
  do {                          /* sections sort loop */               \
    l = index; r = l +n -1;     /* start at left and right boundary */ \
    a = array[*l];              /* get the first and last elements */  \
    z = array[*r];              /* and bring them into right order */  \
    if (a > z) { t = *l; *l = *r; *r = t; }                            \
    x = index[n /2];            /* get the middle element as pivot */  \
    p = array[x];               /* and array element referred to */    \
    if      (p < a) { p = a; x = *l; }  /* compute median of three */  \
    else if (p > z) { p = z; x = *r; }  /* to find a better pivot */   \
    while (1) {                 /* split and exchange loop */          \
      while (array[*++l] < p);  /* skip smaller elems. on the left */  \
      while (array[*--r] > p);  /* skip greater elems. on the right */ \
      if (l >= r) {             /* if less than two elements left, */  \
        if (l <= r) { l++; r--; } break; }       /* abort the loop */  \
      t = *l; *l = *r; *r = t;  /* otherwise exchange elements */      \
    }                                                                  \
    m = (size_t)(index +n -l);  /* compute the number of elements */   \
    n = (size_t)(r -index +1);  /* right and left of the split */      \
    if (n > m) {                /* if right section is smaller, */     \
      if (m >= TH_INSERT)       /* but larger than the threshold, */   \
        name##_qrec(l,     m, array); }   /* sort it recursively, */   \
    else {                      /* if the left section is smaller, */  \
      if (n >= TH_INSERT)       /* but larger than the threshold, */   \
        name##_qrec(index, n, array);     /* sort it recursively, */   \
      index = l; n = m;         /* then switch to the right section */ \
    }                           /* keeping its size m in variable n */ \
  } while (n >= TH_INSERT);     /* while greater than threshold */     \
}  /* qrec() */                                                        \
                                                                       \
/*------------------------------------------------------------------*/ \
                                                                       \
void name##_qsort (tidx *index, size_t n, int dir, const type *array)  \
{                               /* --- sort an index array */          \
  size_t i, k;                  /* loop variable, first section */     \
  tidx   *l, *r;                /* to traverse the array */            \
  tidx   x;                     /* exchange buffer */                  \
  type   t;                     /* buffer for element referred to */   \
                                                                       \
  assert(index && array);       /* check function arguments */         \
  if (n <= 1) return;           /* do not sort less than two elems. */ \
  if (n < TH_INSERT)            /* if less elements than threshold */  \
    k = n;                      /* for insertion sort, note the */     \
  else {                        /* number of elements, otherwise */    \
    name##_qrec(index, n, array);     /* call recursive function */    \
    k = TH_INSERT -1;           /* and get the number of elements */   \
  }                             /* in the first array section */       \
  for (l = r = index; --k > 0;) /* find the position */                \
    if (array[*++r] < array[*l])/* of the smallest element */          \
      l = r;                    /* within the first k elements */      \
  r = index;                    /* swap the smallest element */        \
  x = *l; *l = *r; *r = x;      /* to front as a sentinel */           \
  for (i = n; --i > 0; ) {      /* standard insertion sort */          \
    t = array[x = *++r];        /* note the number to insert */        \
    for (l = r; array[*--l] > t; ) /* shift right all that are */      \
      l[1] = *l;                /* greater than the one to insert */   \
    l[1] = x;                   /* and store the number to insert */   \
  }                             /* in the place thus found */          \
  if (dir < 0)                  /* if descending order requested, */   \
    tin##_reverse(index, n);    /* reverse the element order */        \
}  /* qsort() */

/*--------------------------------------------------------------------*/

IDX_QSORT(i2i, int, int,    int)
IDX_QSORT(i2l, int, int,    long)
IDX_QSORT(i2x, int, int,    diff_t)
IDX_QSORT(i2z, int, int,    size_t)
IDX_QSORT(i2f, int, int,    float)
IDX_QSORT(i2d, int, int,    double)

IDX_QSORT(l2i, lng, long,   int)
IDX_QSORT(l2l, lng, long,   long)
IDX_QSORT(l2x, lng, long,   diff_t)
IDX_QSORT(l2z, lng, long,   size_t)
IDX_QSORT(l2f, lng, long,   float)
IDX_QSORT(l2d, lng, long,   double)

IDX_QSORT(x2i, dif, diff_t, int)
IDX_QSORT(x2l, dif, diff_t, long)
IDX_QSORT(x2x, dif, diff_t, diff_t)
IDX_QSORT(x2z, dif, diff_t, size_t)
IDX_QSORT(x2f, dif, diff_t, float)
IDX_QSORT(x2d, dif, diff_t, double)

/*--------------------------------------------------------------------*/

#define IDX_HEAPSORT(name,tin,tidx,type)                               \
static void name##_sift (tidx *index, size_t l, size_t r,              \
                         const type *array)                            \
{                               /* --- let element sift down in heap */\
  size_t i;                     /* index of first successor in heap */ \
  tidx   x;                     /* buffer for an index element */      \
  type   t;                     /* buffer for an array element */      \
                                                                       \
  t = array[x = index[l]];      /* note the sift element */            \
  i = l +l +1;                  /* compute index of first successor */ \
  do {                          /* sift loop */                        \
    if ((i < r) && (array[index[i]] < array[index[i+1]]))              \
      i++;                      /* if second successor is greater */   \
    if (t >= array[index[i]])   /* if the successor is greater */      \
      break;                    /* than the sift element, */           \
    index[l] = index[i];        /* let the successor ascend in heap */ \
    l = i; i += i +1;           /* compute index of first successor */ \
  } while (i <= r);             /* while still within heap */          \
  index[l] = x;                 /* store the sift element */           \
}  /* sift() */                                                        \
                                                                       \
/*------------------------------------------------------------------*/ \
                                                                       \
void name##_heapsort (tidx *index, size_t n, int dir,                  \
                      const type *array)                               \
{                               /* --- heap sort for index arrays */   \
  size_t l, r;                  /* boundaries of heap section */       \
  tidx   t;                     /* exchange buffer */                  \
                                                                       \
  assert(index && array);       /* check the function arguments */     \
  if (n <= 1) return;           /* do not sort less than 2 elements */ \
  l = n /2;                     /* at start, only the second half */   \
  r = n -1;                     /* of the array has heap structure */  \
  while (l > 0)                 /* while the heap is not complete, */  \
    name##_sift(index, --l, r, array); /* extend it by one element */  \
  while (1) {                   /* heap reduction loop */              \
    t        = index[0];        /* swap the greatest element */        \
    index[0] = index[r];        /* to the end of the array */          \
    index[r] = t;                                                      \
    if (--r <= 0) break;        /* if the heap is empty, abort */      \
    name##_sift(index, 0, r, array);                                   \
  }                             /* let the swap element sift down */   \
  if (dir < 0)                  /* if descending order requested, */   \
    tin##_reverse(index, n);    /* reverse the element order */        \
}  /* heapsort() */

/*--------------------------------------------------------------------*/

IDX_HEAPSORT(i2i, int, int,    int)
IDX_HEAPSORT(i2l, int, int,    long)
IDX_HEAPSORT(i2x, int, int,    diff_t)
IDX_HEAPSORT(i2z, int, int,    size_t)
IDX_HEAPSORT(i2f, int, int,    float)
IDX_HEAPSORT(i2d, int, int,    double)

IDX_HEAPSORT(l2i, lng, long,   int)
IDX_HEAPSORT(l2l, lng, long,   long)
IDX_HEAPSORT(l2x, lng, long,   diff_t)
IDX_HEAPSORT(l2z, lng, long,   size_t)
IDX_HEAPSORT(l2f, lng, long,   float)
IDX_HEAPSORT(l2d, lng, long,   double)

IDX_HEAPSORT(x2i, dif, diff_t, int)
IDX_HEAPSORT(x2l, dif, diff_t, long)
IDX_HEAPSORT(x2x, dif, diff_t, diff_t)
IDX_HEAPSORT(x2z, dif, diff_t, size_t)
IDX_HEAPSORT(x2f, dif, diff_t, float)
IDX_HEAPSORT(x2d, dif, diff_t, double)

/*--------------------------------------------------------------------*/

#define I2P_QSORT(name,tin,tidx)                                       \
static void name##_qrec (tidx *index, size_t n,                        \
                         const void **array, CMPFN *cmp, void *data)   \
{                               /* --- recursive part of quicksort */  \
  tidx   *l, *r;                /* pointers to exchange positions */   \
  tidx   x,  t;                 /* pivot element and exchange buffer */\
  size_t m;                     /* number of elements in 2nd section */\
  const void *p, *a, *z;        /* buffers for array elements */       \
                                                                       \
  do {                          /* sections sort loop */               \
    l = index; r = l +n -1;     /* start at left and right boundary */ \
    a = array[*l];              /* get the first and last elements */  \
    z = array[*r];              /* and bring them into right order */  \
    if (cmp(a, z, data) > 0) { t = *l; *l = *r; *r = t; }              \
    x = index[n /2];            /* get the middle element as pivot */  \
    p = array[x];               /* compute median of 3 to improve */   \
    if      (cmp(p, a, data) < 0) { p = a; x = *l; }                   \
    else if (cmp(p, z, data) > 0) { p = z; x = *r; }                   \
    while (1) {                 /* split and exchange loop */          \
      while (cmp(array[*++l], p, data) < 0)                            \
        ;                       /* skip elements smaller than pivot */ \
      while (cmp(array[*--r], p, data) > 0)                            \
        ;                       /* skip elements greater than pivot */ \
      if (l >= r) {             /* if less than two elements left, */  \
        if (l <= r) { l++; r--; } break; }       /* abort the loop */  \
      t = *l; *l = *r; *r = t;  /* otherwise exchange elements */      \
    }                                                                  \
    m = (size_t)(index +n -l);  /* compute the number of elements */   \
    n = (size_t)(r -index +1);  /* right and left of the split */      \
    if (n > m) {                /* if right section is smaller, */     \
      if (m >= TH_INSERT)       /* but larger than the threshold, */   \
        name##_qrec(l,     m, array, cmp, data); }    /* sort it, */   \
    else {                      /* if the left section is smaller, */  \
      if (n >= TH_INSERT)       /* but larger than the threshold, */   \
        name##_qrec(index, n, array, cmp, data);      /* sort it, */   \
      index = l; n = m;         /* then switch to the right section */ \
    }                           /* keeping its size m in variable n */ \
  } while (n >= TH_INSERT);     /* while greater than threshold */     \
}  /* qrec() */                                                        \
                                                                       \
/*------------------------------------------------------------------*/ \
                                                                       \
void name##_qsort (tidx *index, size_t n, int dir,                     \
                   const void **array, CMPFN *cmp, void *data)         \
{                               /* --- quick sort for index arrays */  \
  size_t i, k;                  /* loop variable, first section */     \
  tidx   *l, *r;                /* to traverse the array */            \
  tidx   x;                     /* exchange buffer */                  \
  const void *t;                /* buffer for pointer referred to */   \
                                                                       \
  assert(index && array && cmp);/* check the function arguments */     \
  if (n <= 1) return;           /* do not sort less than two elems. */ \
  if (n < TH_INSERT)            /* if fewer elements than threshold */ \
    k = n;                      /* for insertion sort, note the */     \
  else {                        /* number of elements, otherwise */    \
    name##_qrec(index, n, array, cmp, data); /* sort recursively */    \
    k = TH_INSERT -1;           /* and get the number of elements */   \
  }                             /* in the first array section */       \
  for (l = r = index; --k > 0;) /* find the smallest element */        \
    if (cmp(array[*++r], array[*l], data) < 0)                         \
      l = r;                    /* among the first k elements */       \
  r = index;                    /* swap the smallest element */        \
  x = *l; *l = *r; *r = x;      /* to the front as a sentinel */       \
  for (i = n; --i > 0; ) {      /* standard insertion sort */          \
    t = array[x = *++r];        /* note the element to insert */       \
    for (l = r; cmp(array[*--l], t, data) > 0; )                       \
      l[1] = *l;                /* shift right index elements */       \
    l[1] = x;                   /* that are greater than the one */    \
  }                             /* to insert and store this element */ \
  if (dir < 0)                  /* if descending order requested, */   \
    tin##_reverse(index, n);    /* reverse the element order */        \
}  /* qsort() */

/*--------------------------------------------------------------------*/

I2P_QSORT(i2p, int, int)
I2P_QSORT(l2p, lng, long)
I2P_QSORT(x2p, dif, diff_t)

/*--------------------------------------------------------------------*/

#define I2P_HEAPSORT(name,tin,tidx)                                    \
static void name##_sift (tidx *index, size_t l, size_t r,              \
                         const void **array, CMPFN *cmp, void *data)   \
{                               /* --- let element sift down in heap */\
  size_t i;                     /* index of first successor in heap */ \
  tidx   x;                     /* buffer for an index element */      \
  const void *t;                /* buffer for an array element */      \
                                                                       \
  t = array[x = index[l]];      /* note the sift element */            \
  i = l +l +1;                  /* compute index of first successor */ \
  do {                          /* sift loop */                        \
    if ((i < r)                 /* if second successor is greater */   \
    &&  (cmp(array[index[i]], array[index[i+1]], data) < 0))           \
      i++;                      /* go to the second successor */       \
    if (cmp(t, array[i], data) >= 0) /* if the successor is greater */ \
      break;                         /* than the sift element, */      \
    index[l] = index[i];        /* let the successor ascend in heap */ \
    l = i; i += i +1;           /* compute index of first successor */ \
  } while (i <= r);             /* while still within heap */          \
  index[l] = x;                 /* store the sift element */           \
}  /* sift() */                                                        \
                                                                       \
/*------------------------------------------------------------------*/ \
                                                                       \
void name##_heapsort (tidx *index, size_t n, int dir,                  \
                      const void **array, CMPFN *cmp, void *data)      \
{                               /* --- heap sort for index arrays */   \
  size_t l, r;                  /* boundaries of heap section */       \
  tidx   t;                     /* exchange buffer */                  \
                                                                       \
  assert(index && array && cmp);/* check the function arguments */     \
  if (n <= 1) return;           /* do not sort less than two elems. */ \
  l = n /2;                     /* at start, only the second half */   \
  r = n -1;                     /* of the array has heap structure */  \
  while (l > 0)                 /* while the heap is not complete, */  \
    name##_sift(index, --l, r, array, cmp, data);     /* extend it */  \
  while (1) {                   /* heap reduction loop */              \
    t        = index[0];        /* swap the greatest element */        \
    index[0] = index[r];        /* to the end of the array */          \
    index[r] = t;                                                      \
    if (--r <= 0) break;        /* if the heap is empty, abort */      \
    name##_sift(index, 0, r, array, cmp, data);                        \
  }                             /* let the swapped element sift down */\
  if (dir < 0)                  /* if descending order requested, */   \
    tin##_reverse(index, n);    /* reverse the element order */        \
}  /* heapsort() */

/*--------------------------------------------------------------------*/

I2P_HEAPSORT(i2p, int, int)
I2P_HEAPSORT(l2p, lng, long)
I2P_HEAPSORT(x2p, dif, diff_t)

/*--------------------------------------------------------------------*/

#define I2C_QSORT(name,tin,tidx)                                       \
static void name##_qrec (tidx *index, size_t n,                        \
                         tin##_CMPFN *cmp, void *data)                 \
{                               /* --- recursive part of quicksort */  \
  tidx   *l, *r;                /* pointers to exchange positions */   \
  tidx   x,  t;                 /* pivot element and exchange buffer */\
  size_t m;                     /* number of elements in 2nd section */\
                                                                       \
  do {                          /* sections sort loop */               \
    l = index; r = l +n -1;     /* start at left and right boundary */ \
    if (cmp(*l, *r, data) > 0){ /* bring the first and last */         \
      t = *l; *l = *r; *r = t;} /* element into proper order */        \
    x = index[n /2];            /* get the middle element as pivot */  \
    if      (cmp(x, *l, data) < 0) x = *l;  /* try to find a */        \
    else if (cmp(x, *r, data) > 0) x = *r;  /* better pivot */         \
    while (1) {                 /* split and exchange loop */          \
      while (cmp(*++l, x, data) < 0)  /* skip left  elements that */   \
        ;                       /* are smaller than pivot element */   \
      while (cmp(*--r, x, data) > 0)  /* skip right elements that */   \
        ;                       /* are greater than pivot element */   \
      if (l >= r) {             /* if less than two elements left, */  \
        if (l <= r) { l++; r--; } break; }       /* abort the loop */  \
      t = *l; *l = *r; *r = t;  /* otherwise exchange elements */      \
    }                                                                  \
    m = (size_t)(index +n -l);  /* compute the number of elements */   \
    n = (size_t)(r -index +1);  /* right and left of the split */      \
    if (n > m) {                /* if right section is smaller, */     \
      if (m >= TH_INSERT)       /* but larger than the threshold, */   \
        name##_qrec(l,     m, cmp, data);} /* sort it recursively, */  \
    else {                      /* if the left section is smaller, */  \
      if (n >= TH_INSERT)       /* but larger than the threshold, */   \
        name##_qrec(index, n, cmp, data);  /* sort it recursively, */  \
      index = l; n = m;         /* then switch to the right section */ \
    }                           /* keeping its size m in variable n */ \
  } while (n >= TH_INSERT);     /* while greater than threshold */     \
}  /* qrec() */                                                        \
                                                                       \
/*------------------------------------------------------------------*/ \
                                                                       \
void name##_qsort (tidx *index, size_t n, int dir,                     \
                   tin##_CMPFN *cmp, void *data)                       \
{                               /* --- quicksort for index arrays */   \
  size_t i, k;                  /* loop variable, first section */     \
  tidx   *l, *r;                /* to traverse the array */            \
  tidx   t;                     /* exchange buffer */                  \
                                                                       \
  assert(index && cmp);         /* check the function arguments */     \
  if (n <= 1) return;           /* do not sort less than two elems. */ \
  if (n < TH_INSERT)            /* if fewer elements than threshold */ \
    k = n;                      /* for insertion sort, note the */     \
  else {                        /* number of elements, otherwise */    \
    name##_qrec(index, n, cmp, data);  /* call recursice function */   \
    k = TH_INSERT -1;           /* and get the number of elements */   \
  }                             /* in the first array section */       \
  for (l = r = index; --k > 0;) /* find the smallest element within */ \
    if (cmp(*++r, *l, data) < 0) l = r;     /* the first k elements */ \
  r = index;                    /* swap the smallest element */        \
  t = *l; *l = *r; *r = t;      /* to the front as a sentinel */       \
  for (i = n; --i > 0; ) {      /* standard insertion sort */          \
    t = *++r;                   /* note the element to insert */       \
    for (l = r; cmp(*--l, t, data) > 0; )   /* shift right elements */ \
      l[1] = *l;                /* that are greater than the one to */ \
    l[1] = t;                   /* insert and store the element to */  \
  }                             /* insert in the place thus found */   \
  if (dir < 0)                  /* if descending order requested, */   \
    tin##_reverse(index, n);    /* reverse the element order */        \
}  /* qsort() */

/*--------------------------------------------------------------------*/

I2C_QSORT(i2c, int, int)
I2C_QSORT(l2c, lng, long)
I2C_QSORT(x2c, dif, diff_t)

/*--------------------------------------------------------------------*/

#define I2C_HEAPSORT(name,tin,tidx)                                    \
static void name##_sift (tidx *index, size_t l, size_t r,              \
                         tin##_CMPFN *cmp, void *data)                 \
{                               /* --- let element sift down in heap */\
  size_t i;                     /* index of first successor in heap */ \
  tidx   t;                     /* buffer for an array element */      \
                                                                       \
  t = index[l];                 /* note the sift element */            \
  i = l +l +1;                  /* compute index of first successor */ \
  do {                          /* sift loop */                        \
    if ((i < r)                 /* if second successor is greater */   \
    &&  (cmp(index[i], index[i+1], data) < 0))                         \
      i++;                      /* go to the second successor */       \
    if (cmp(t, index[i], data) >= 0) /* if the successor is greater */ \
      break;                         /* than the sift element, */      \
    index[l] = index[i];        /* let the successor ascend in heap */ \
    l = i; i += i +1;           /* compute index of first successor */ \
  } while (i <= r);             /* while still within heap */          \
  index[l] = t;                 /* store the sift element */           \
}  /* sift() */                                                        \
                                                                       \
/*------------------------------------------------------------------*/ \
                                                                       \
void name##_heapsort (tidx *index, size_t n, int dir,                  \
                      tin##_CMPFN *cmp, void *data)                    \
{                               /* --- heap sort for index arrays */   \
  size_t l, r;                  /* boundaries of heap section */       \
  tidx   t;                     /* exchange buffer */                  \
                                                                       \
  assert(index && cmp);         /* check the function arguments */     \
  if (n <= 1) return;           /* do not sort less than two elems. */ \
  l = n /2;                     /* at start, only the second half */   \
  r = n -1;                     /* of the index has heap structure */  \
  while (l > 0)                 /* while the heap is not complete, */  \
    name##_sift(index,--l, r, cmp, data); /* extend it by one elem. */ \
  while (1) {                   /* heap reduction loop */              \
    t        = index[0];        /* swap the greatest element */        \
    index[0] = index[r];        /* to the end of the index */          \
    index[r] = t;                                                      \
    if (--r <= 0) break;        /* if the heap is empty, abort */      \
    name##_sift(index, 0, r, cmp, data);                               \
  }                             /* let the swapped element sift down */\
  if (dir < 0)                  /* if descending order requested, */   \
    tin##_reverse(index, n);    /* reverse the element order */        \
}  /* heapsort() */

/*--------------------------------------------------------------------*/

I2C_HEAPSORT(i2c, int, int)
I2C_HEAPSORT(l2c, lng, long)
I2C_HEAPSORT(x2c, dif, diff_t)

/*----------------------------------------------------------------------
  Main Function for Testing
----------------------------------------------------------------------*/
#ifdef ARRAYS_MAIN

static int lexcmp (const void *p1, const void *p2, void *data)
{                               /* --- compare lexicographically */
  return strcmp(p1, p2);        /* use standard string comparison */
}  /* lexcmp() */

/*--------------------------------------------------------------------*/

static int numcmp (const void *p1, const void *p2, void *data)
{                               /* --- compare strings numerically */
  double d1 = strtod((const char*)p1, NULL);
  double d2 = strtod((const char*)p2, NULL);
  if (d1 < d2) return -1;       /* convert to numbers and */
  if (d1 > d2) return +1;       /* compare numerically */
  return strcmp(p1, p2);        /* if the numbers are equal, */
}  /* numcmp() */               /* compare strings lexicographically */

/*--------------------------------------------------------------------*/

int main (int argc, char *argv[])
{                               /* --- sort program arguments */
  int   i, n;                   /* loop variables */
  char  *s;                     /* to traverse the arguments */
  CMPFN *cmp = lexcmp;          /* comparison function */

  if (argc < 2) {               /* if no arguments are given */
    printf("usage: %s [arg [arg ...]]\n", argv[0]);
    printf("sort the list of program arguments\n");
    return 0;                   /* print a usage message */
  }                             /* and abort the program */
  for (i = n = 0; ++i < argc; ) {
    s = argv[i];                /* traverse the arguments */
    if (*s != '-') { argv[n++] = s; continue; }
    s++;                        /* store the arguments to sort */
    while (*s) {                /* traverse the options */
      switch (*s++) {           /* evaluate the options */
        case 'n': cmp = numcmp; break;
        default : printf("unknown option -%c\n", *--s); return -1;
      }                         /* set the option variables */
    }                           /* and check for known options */
  }
  ptr_qsort(argv, (size_t)n, +1, cmp, NULL);
  for (i = 0; i < n; i++) {     /* sort and print program arguments */
    fputs(argv[i], stdout); fputc('\n', stdout); }
  return 0;                     /* return 'ok' */
}  /* main() */

#endif
