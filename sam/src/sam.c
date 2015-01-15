/*----------------------------------------------------------------------
  File    : sam.c
  Contents: split and merge algorithm for finding frequent item sets
  Author  : Christian Borgelt
  History : 2008.10.15 file created from file relim.c
            2008.10.25 first version completed
            2008.11.10 framework for item insertions added
            2008.11.11 functions for unlimited item insertion added
            2008.11.13 functions for limited   item insertion added
            2008.11.18 handling of vanishing penalties improved
            2008.12.05 perfect extension pruning added (optional)
            2008.12.20 optional binary search based merging added
            2009.05.28 adapted to modified function tbg_filter()
            2009.10.15 adapted to item set counter in reporter
            2009.10.16 closed and maximal item set mining added
            2010.03.18 recording of actual number of occurrences added
            2010.04.07 threshold for both item set support and weight
            2010.07.14 output file made optional (for benchmarking)
            2010.08.19 item selection file added as optional input
            2010.08.22 adapted to modified modules tabread and tract
            2010.10.15 adapted to modified interface of module report
            2010.11.03 bug in function rec_lim() fixed (variant sam_lim)
            2010.11.05 clearer interpretation of minimum support
            2010.11.24 adapted to modified error reporting (tract)
            2010.12.11 adapted to a generic error reporting function
            2011.03.16 closed/maximal item sets with item insertions
            2011.03.20 optional integer transaction weights added
            2011.05.30 item weight combination with t-norms added
            2011.07.08 adapted to modified function tbg_recode()
            2011.08.11 bug in call/parameter list of sam_lim() fixed
            2011.08.28 output of item set counters per size added
            2011.08.29 16 items machine added (without item insertions)
            2011.12.02 version using a transaction prefix tree added
            2013.04.01 adapted to type changes in module tract
            2013.10.15 checks of return code of isr_report() added
            2013.10.18 optional pattern spectrum collection added
            2013.11.12 item insertion penalties changed to option -R#
            2013.11.20 function sam() added, bug in rec_tree() fixed
            2014.05.12 option -F# added (support border for filtering)
            2014.08.02 option -c renamed to -i, option -a renamed to -A
            2014.08.22 adapted to modified item set reporter interface
            2014.08.28 functions sam_data() and sam_repo() added
            2014.10.24 changed from LGPL license to MIT license
------------------------------------------------------------------------
  Reference for the SaM algorithm:
    C. Borgelt and X. Wang.
    SaM: A Split and Merge Algorithm for Fuzzy Frequent Item Set Mining.
    Proc. 13th Int. Fuzzy Systems Association World Congress and
    6th Conf. of the European Society for Fuzzy Logic and Technology
    (IFSA/EUSFLAT'09, Lisbon, Portugal), 968--973.
    IFSA/EUSFLAT Organization Committee, Lisbon, Portugal 2009
----------------------------------------------------------------------*/
#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <limits.h>
#include <float.h>
#include <math.h>
#include <time.h>
#include <assert.h>
#ifndef ISR_PATSPEC
#define ISR_PATSPEC
#endif
#ifdef SAM_MAIN
#ifndef PSP_REPORT
#define PSP_REPORT
#endif
#ifndef TA_READ
#define TA_READ
#endif
#endif
#include "sam.h"
#include "fim16.h"
#ifdef SAM_MAIN
#include "error.h"
#endif
#ifdef STORAGE
#include "storage.h"
#endif

#ifndef INFINITY
#define INFINITY    (DBL_MAX+DBL_MAX)
#endif                          /* MSC still does not support C99 */

/*----------------------------------------------------------------------
  Preprocessor Definitions
----------------------------------------------------------------------*/
#define PRGNAME     "sam"
#define DESCRIPTION "find frequent item sets " \
                    "with a split and merge algorithm"
#define VERSION     "version 3.11 (2014.10.24)        " \
                    "(c) 2008-2014   Christian Borgelt"

/* --- error codes --- */
/* error codes   0 to  -4 defined in tract.h */
#define E_STDIN      (-5)       /* double assignment of stdin */
#define E_OPTION     (-6)       /* unknown option */
#define E_OPTARG     (-7)       /* missing option argument */
#define E_ARGCNT     (-8)       /* too few/many arguments */
#define E_TARGET     (-9)       /* invalid target type */
#define E_SIZE      (-10)       /* invalid item set size */
#define E_SUPPORT   (-11)       /* invalid minimum item set support */
#define E_VARIANT   (-12)       /* invalid algorithm variant */
#define E_WEIGHT    (-13)       /* invalid minimum transaction weight */
#define E_MEASURE   (-14)       /* invalid evaluation measure */
#define E_TNORM     (-16)       /* invalid triangular norm */
/* error codes -15 to -26 defined in tract.h */

#define COPYERR       ((TTNODE*)-1)

#ifndef QUIET                   /* if not quiet version, */
#define MSG         fprintf     /* print messages */
#define XMSG        if (mode & SAM_VERBOSE) fprintf
#else                           /* if quiet version, */
#define MSG(...)                /* suppress messages */
#define XMSG(...)
#endif

#define SEC_SINCE(t)  ((double)(clock()-(t)) /(double)CLOCKS_PER_SEC)

/*----------------------------------------------------------------------
  Type Definitions
----------------------------------------------------------------------*/
typedef double TNORM (double a, double b);

typedef struct {                /* --- transaction array element --- */
  const ITEM *items;            /* items in the transaction */
  SUPP       occ;               /* number of occurrences */
} TAAE;                         /* (transaction array element) */

typedef struct {                /* --- transaction array element --- */
  const ITEM *items;            /* items in the transaction */
  SUPP       occ;                /* actual number of occurrences */
  double     wgt;                /* weight of transactions */
} TXAE;                         /* (transaction array element) */

typedef struct {                /* --- transaction array element --- */
  const ITEM *items;            /* items in the transaction */
  SUPP       occ;               /* actual number of occurrences */
  SUPP       cnt;               /* number of transactions */
  double     wgt;               /* weight per transaction */
} TZAE;                         /* (transaction array element) */

typedef struct ttnode {         /* --- transaction tree node --- */
  ITEM          item;           /* associated item (last item in set) */
  SUPP          supp;           /* support of represented item set */
  struct ttnode *children;      /* list of child nodes */
  struct ttnode *sibling;       /* successor node in sibling list */
} TTNODE;                       /* (transaction tree node) */

typedef struct {                /* --- recursion data --- */
  int       mode;               /* operation mode (e.g. pruning) */
  SUPP      smin;               /* minimum support of an item set */
  double    sins;               /* minimum support with insertions */
  double    twgt;               /* minimum transaction weight */
  MEMSYS    *mem;               /* memory system for tree version */
  TNORM     *tnorm;             /* t-norm for comb. item penalties */
  FIM16     *fim16;             /* 16 items machine */
  TID       merge;              /* threshold for source merging */
  void      *buf;               /* buffer for projection */
  ITEMBASE  *base;              /* underlying item base */
  ISREPORT  *report;            /* item set reporter */
} RECDATA;                      /* (recursion data) */

/*----------------------------------------------------------------------
  Constants
----------------------------------------------------------------------*/
#if !defined QUIET && defined SAM_MAIN
/* --- error messages --- */
static const char *errmsgs[] = {
  /* E_NONE      0 */  "no error",
  /* E_NOMEM    -1 */  "not enough memory",
  /* E_FOPEN    -2 */  "cannot open file %s",
  /* E_FREAD    -3 */  "read error on file %s",
  /* E_FWRITE   -4 */  "write error on file %s",
  /* E_STDIN    -5 */  "double assignment of standard input",
  /* E_OPTION   -6 */  "unknown option -%c",
  /* E_OPTARG   -7 */  "missing option argument",
  /* E_ARGCNT   -8 */  "wrong number of arguments",
  /* E_TARGET   -9 */  "invalid target type '%c'",
  /* E_SIZE    -10 */  "invalid item set size %"ITEM_FMT,
  /* E_SUPPORT -11 */  "invalid minimum support %g",
  /* E_VARIANT -12 */  "invalid sam variant '%c'",
  /* E_WEIGHT  -13 */  "invalid minimum transaction weight %g",
  /* E_MEASURE -14 */  "invalid evaluation measure '%c'",
  /* E_NOITEMS -15 */  "no (frequent) items found",
  /* E_TNORM   -16 */  "invalid triangular norm '%c'",
  /*           -18 */  "unknown error"
};
#endif

/*----------------------------------------------------------------------
  Global Variables
----------------------------------------------------------------------*/
#ifdef SAM_MAIN
#ifndef QUIET
static CCHAR    *prgname;       /* program name for error messages */
#endif
static TABREAD  *tread  = NULL; /* table/transaction reader */
static ITEMBASE *ibase  = NULL; /* item base */
static TABAG    *tabag  = NULL; /* transaction bag/multiset */
static ISREPORT *report = NULL; /* item set reporter */
static TABWRITE *twrite = NULL; /* table writer for pattern spectrum */
static double   *border = NULL; /* support border for filtering */
#endif

/*----------------------------------------------------------------------
  Auxiliary Functions (for debugging)
----------------------------------------------------------------------*/
#if !defined NDEBUG && defined SAM_MAIN

static void indent (int k)
{ while (--k >= 0) printf("   "); }

/*--------------------------------------------------------------------*/

#define SHOW(show,type,weight,fmt,term) \
static void show (const char *text, type *a, int ind)                  \
{                               /* --- show a transaction array */     \
  ITEM       n;                 /* item counter */                     \
  const ITEM *p;                /* to traverse the item array */       \
  weight     w;                 /* total weight */                     \
                                                                       \
  indent(ind);                  /* indent the output line */           \
  printf("%s\n", text);         /* print the given text */             \
  for (w = 0, n = 0; a->items; a++) { /* traverse the transactions */  \
    indent(ind);                /* indent the output line */           \
    OCCUR;                      /* print the number of occurrences */  \
    printf(fmt" :", term);      /* print the transaction weight */     \
    for (p = a->items; *p >= 0; p++)                                   \
      printf(" %s", ib_name(ibase, *p));                               \
    printf("\n");               /* print items in the transaction */   \
    w += term; n++;             /* sum the transaction weights */      \
  }                             /* and count the transactions */       \
  indent(ind);                  /* indent the output line */           \
  printf("total: %"ITEM_FMT"/"fmt"\n\n", n, w);                        \
}  /* show() */                 /* print total transaction weight */

/*--------------------------------------------------------------------*/

#define OCCUR
SHOW(show,     TAAE, SUPP, "%"SUPP_FMT, a->occ)
#undef  OCCUR
#define OCCUR  printf("%2"SUPP_FMT"/", a->occ)
SHOW(show_ins, TXAE, double, "%5.2f", a->wgt)
#undef  OCCUR
#define OCCUR  printf("%4"SUPP_FMT"/%4"TID_FMT"/%7.2f/", \
                      a->occ, a->cnt, a->wgt)
SHOW(show_lim, TZAE, double, "%7.2f", a->wgt *(double)a->cnt)

/*--------------------------------------------------------------------*/

static void show_tree (TTNODE *node, int ind)
{                               /* --- show a freq. pattern tree */
  ITEM i;                       /* mapped item identifier */

  assert(ind >= 0);             /* check the function arguments */
  while (node) {                /* traverse the node list */
    indent(ind);                /* indent the output line */
    i = node->item;             /* print the item name/bit pattern */
    if (i < 0) printf("%04x",  (int)(i & 0xffff));
    else       printf("%s/%"ITEM_FMT, ib_name(ibase, i), i);
    printf(":%"SUPP_FMT"\n",node->supp);/* print the node information */
    show_tree(node->children, ind+1);
    node = node->sibling;       /* recursively show the child nodes, */
  }                             /* then go to the next node */
}  /* show_tree() */

#endif
/*----------------------------------------------------------------------
  Triangular Norms (t-norms)
----------------------------------------------------------------------*/

static double t_min  (double a, double b)
{ return (a < b) ? a : b; }

static double t_nilp (double a, double b)
{ return (a+b <= 1) ? 0 : (a < b) ? a : b; }

static double t_prod (double a, double b)
{ return a*b; }

static double t_luka (double a, double b)
{ double x = a+b-1; return (x > 0) ? x : 0; }

static double t_hama (double a, double b)
{ double x = a+b-a*b; return (x > 0) ? (a*b)/x : 0; }

/*--------------------------------------------------------------------*/

static TNORM *tnorms[] = {      /* t-norms (triangular norms) */
  /* T_MIN    0 */  t_min,      /* minimum */
  /* T_NILP   1 */  t_nilp,     /* nil-potent minimum */
  /* T_PROD   2 */  t_prod,     /* product */
  /* T_LUKA   3 */  t_luka,     /* Lukasiewicz */
  /* T_HAMA   4 */  t_hama,     /* Hamacher product */
};

/*----------------------------------------------------------------------
  Comparison Function
----------------------------------------------------------------------*/

static int cmp (const ITEM *a, const ITEM *b)
{                               /* --- compare two transactions */
  assert(a && b);               /* check the function arguments */
  for ( ; 1; a++, b++) {        /* lexicographic comparison loop */
    if (*a < *b) return -1;     /* compare corresponding items */
    if (*a > *b) return +1;     /* and if one is greater, abort */
    if (*a <= TA_END) return 0; /* otherwise check for the sentinel */
  }                             /* and abort if it is reached */
}  /* cmp() */

/*----------------------------------------------------------------------
  Split and Merge: Basic Version
----------------------------------------------------------------------*/

static int recurse (TAAE *a, TID n, RECDATA *rd)
{                               /* --- split and merge recursion */
  int  r = 0;                   /* error status */
  ITEM i;                       /* current item */
  TAAE *proj;                   /* projected transaction database */
  TAAE *s, *t, *d;              /* to traverse the transactions */
  SUPP supp;                    /* support of (current) split item */
  SUPP pex;                     /* minimum support for perfect exts. */

  assert(a && (n > 0) && rd);   /* check the function arguments */
  pex  = (rd->mode & SAM_PERFECT) ? isr_supp(rd->report) : SUPP_MAX;
  proj = (TAAE*)malloc((size_t)(n+1) *sizeof(TAAE));
  if (!proj) return -1;         /* allocate the projection array */
  while (a->items) {            /* split and merge loop */
    i = a->items[0];            /* get the next split item */
    if (i < 0) {                /* if only packed items left */
      do { m16_add(rd->fim16, (BITTA)(a->items[0] & ~TA_END), a->occ);
      } while ((++a)->items);   /* add trans. to 16 items machine */
      r = m16_mine(rd->fim16);  /* mine with 16 items machine */
      break;                    /* and abort the search loop */
    }
    d = proj; s = a; supp = 0;  /* -- split the transaction array */
    while (s->items && (s->items[0] == i)) {
      d->items = ++s->items;    /* copy trans. with current item */
      supp += (d++)->occ = (s++)->occ;
    }                           /* sum the item occurrences */
    if (supp >= pex) {          /* if item is a perfect extension, */
      isr_addpex(rd->report,i); /* add it to the item set reporter */
      if ((s-1)->items[0] <= TA_END) (--s)->items = NULL;
      continue;                 /* remove an empty transaction, */
    }                           /* store a sentinel, and skip item */
    if ((d-1)->items[0] <= TA_END)
      --d;                      /* remove an empty transaction, */
    d->items = NULL;            /* store a sentinel at the end and */
    n = (TID)(d -proj);         /* compute the number of elements */
    d = a; t = proj;            /* -- merge the transaction arrays */
    while (s->items && t->items) {         /* compare transactions */
      int c = cmp(s->items, t->items);     /* from the two sources */
      if (c) *d++ = (c > 0) ? *s++ : *t++; /* and copy the smaller */
      else { *d = *s++; (d++)->occ += (t++)->occ; }
    }                           /* combine equal transactions */
    while (t->items) *d++ = *t++;  /* copy remaining transactions */
    while (s->items) *d++ = *s++;  /* from the non-empty source */
    d->items = NULL;            /* store a sentinel at the end */
    if (supp < rd->smin)        /* if the support is too low, */
      continue;                 /* skip the recursive processing */
    r = isr_add(rd->report, i, supp);
    if (r <  0) break;          /* add current item to the reporter */
    if (r <= 0) continue;       /* and check if it needs processing */
    if ((n > 0)                 /* if the projection is not empty */
    && isr_xable(rd->report,1)){/* and another item can be added, */
      r = recurse(proj, n, rd); /* search projection recursively */
      if (r < 0) break;         /* abort on a recursion error */
    }
    r = isr_report(rd->report); /* report the current item set */
    if (r < 0) break;           /* and check for an error */
    isr_remove(rd->report, 1);  /* remove the current item */
  }                             /* from the item set reporter */
  free(proj);                   /* deallocate the projection array */
  return r;                     /* return the error status */
}  /* recurse() */

/*--------------------------------------------------------------------*/

int sam_base (TABAG *tabag, SUPP smin, int mode, ISREPORT *report)
{                               /* --- split and merge algorithm */
  int     r;                    /* result of recursion */
  ITEM    k;                    /* number of items */
  TID     i, n;                 /* loop variable, number of trans. */
  TRACT   *t;                   /* to traverse the transactions */
  TAAE    *a;                   /* initial transaction array */
  RECDATA rd;                   /* recursion data */

  assert(tabag && report);      /* check the function arguments */
  rd.mode = mode;               /* make search mode consistent */
  rd.smin = (smin > 0) ? smin : 1; /* check and adapt the support */
  if (tbg_wgt(tabag) < rd.smin) /* check the total transaction weight */
    return 0;                   /* against the minimum support */
  k = tbg_itemcnt(tabag);       /* get the number of items */
  if (k <= 0) return isr_report(report);
  n = tbg_cnt(tabag);           /* get the number of transactions */
  a = (TAAE*)malloc((size_t)(n+1) *sizeof(TAAE));
  if (!a) return -1;            /* create initial transaction array */
  for (i = n; --i >= 0; ) {     /* traverse the transactions */
    t = tbg_tract(tabag, i);    /* initialize a new array element */
    a[i].items = ta_items(t);   /* for each transaction and */
    a[i].occ   = ta_wgt(t);     /* set the item array (pointer) */
  }                             /* and the transaction weight */
  if (a[n-1].items[0] <= TA_END)
    n--;                        /* remove an empty transaction */
  a[n].items = NULL;            /* store a sentinel at the end */
  rd.fim16   = NULL;            /* default: no 16 items machine */
  if (mode & SAM_FIM16) {       /* if to use a 16 items machine */
    rd.fim16 = m16_create(-1, rd.smin, report);
    if (!rd.fim16) { free(a); return -1; }
  }                             /* create a 16 items machine */
  rd.report = report;           /* note the item set reporter */
  r = recurse(a, n, &rd);       /* execute split and merge recursion */
  if (rd.fim16)                 /* if a 16 items machine was used, */
    m16_delete(rd.fim16);       /* delete the 16 items machine */
  free(a);                      /* deallocate the transaction array */
  if (r >= 0)                   /* if no error occurred, */
    r = isr_report(report);     /* report the empty item set */
  return r;                     /* return the error status */
}  /* sam_base() */

/*----------------------------------------------------------------------
  Split and Merge: Optimized Merge (Binary Search)
----------------------------------------------------------------------*/

static int rec_opt (TAAE *a, TID n, RECDATA *rd)
{                               /* --- split and merge recursion */
  int  r = 0;                   /* error status */
  ITEM i;                       /* current item */
  TAAE *proj;                   /* projected transaction database */
  TAAE *s, *t, *d;              /* to traverse the transactions */
  SUPP supp;                    /* support of (current) split item */
  SUPP pex;                     /* minimum support for perfect exts. */
  TID  k;                       /* number of remaining transactions */

  assert(a && (n > 0) && rd);   /* check the function arguments */
  pex  = (rd->mode & SAM_PERFECT) ? isr_supp(rd->report) : SUPP_MAX;
  proj = (TAAE*)malloc((size_t)(n+1) *sizeof(TAAE));
  if (!proj) return -1;         /* allocate the projection array */
  for (k = n; a->items; ) {     /* split and merge loop */
    i = a->items[0];            /* get the next split item */
    if (i < 0) {                /* if only packed items left */
      do { m16_add(rd->fim16, (BITTA)(a->items[0] & ~TA_END), a->occ);
      } while ((++a)->items);   /* add trans. to 16 items machine */
      r = m16_mine(rd->fim16);  /* mine with 16 items machine */
      break;                    /* and abort the search loop */
    }
    d = proj; s = a; supp = 0;  /* -- split the transaction array */
    while (s->items && (s->items[0] == i)) {
      d->items = ++s->items;    /* copy trans. with current item */
      supp += (d++)->occ = (s++)->occ;
    }                           /* sum the item occurrences */
    if (supp >= pex) {          /* if item is a perfect extension, */
      isr_addpex(rd->report,i); /* add it to the item set reporter */
      if ((s-1)->items[0] <= TA_END) { (--s)->items = NULL; --k; }
      continue;                 /* remove an empty transaction, */
    }                           /* store a sentinel, and skip item */
    k -= (TID)(d -proj);        /* remove transferred transactions */
    if ((d-1)->items[0] <= TA_END)
      --d;                      /* remove an empty transaction, */
    d->items = NULL;            /* store a sentinel at the end and */
    n = (TID)(d -proj);         /* compute the number of elements */
    d = a; t = proj;            /* -- merge the transaction arrays */
    if (((n << 4) > k)          /* if the size difference is limited */
    ||  !(rd->mode & SAM_BSEARCH)) {  /* or no binary search merging */
      while (s->items && t->items) {         /* compare transactions */
        int c = cmp(s->items, t->items);     /* from the two sources */
        if (c) *d++ = (c > 0) ? *s++ : *t++; /* and copy the smaller */
        else { *d = *s++; (d++)->occ += (t++)->occ; }
      } }                       /* combine equal transactions */
    else {                      /* if sizes differ considerably */
      while (t->items && k) {   /* while transactions to insert */
        TID l, m, r, x = 1;     /* binary search variables */
        for (l = 0, r = k; l < r; ) {
          int c = cmp(t->items, s[m = (l+r) >> 1].items);
          if (c < 0)  l = m+1;  /* find the insertion position */
          else x = c, r = m;    /* in the remaining transactions */
        }                       /* of the other source */
        for (k -= l; --l >= 0; )/* copy transactions before */
          *d++ = *s++;          /* the insertion position */
        *d++ = *t++;            /* copy transaction to insert */
        if (!x) { (d-1)->occ += (s++)->occ; k--; }
      }                         /* if equal transaction in source, */
    }                           /* combine it with ins. transaction */
    while (t->items) *d++ = *t++;  /* copy remaining transactions */
    while (s->items) *d++ = *s++;  /* from the non-empty source */
    d->items = NULL;            /* store a sentinel at the end and */
    k = (TID)(d -a);            /* compute the number of elements */
    if (supp < rd->smin)        /* if the support is too low, */
      continue;                 /* skip the recursive processing */
    r = isr_add(rd->report, i, supp);
    if (r <  0) break;          /* add current item to the reporter */
    if (r <= 0) continue;       /* and check if it needs processing */
    if ((n > 0)                 /* if the projection is not empty */
    && isr_xable(rd->report,1)){/* and another item can be added */
      r = rec_opt(proj, n, rd); /* search projection recursively */
      if (r < 0) break;         /* abort on a recursion error */
    }
    r = isr_report(rd->report); /* report the current item set */
    if (r < 0) break;           /* and check for an error */
    isr_remove(rd->report, 1);  /* remove the current item */
  }                             /* from the item set reporter */
  free(proj);                   /* deallocate the projection array */
  return r;                     /* return the error status */
}  /* rec_opt() */

/*--------------------------------------------------------------------*/

int sam_opt (TABAG *tabag, SUPP smin, int mode, ISREPORT *report)
{                               /* --- split and merge algorithm */
  int     r;                    /* result of recursion */
  ITEM    k;                    /* number of items */
  TID     i, n;                 /* loop variable, number of trans. */
  TRACT   *t;                   /* to traverse the transactions */
  TAAE    *a;                   /* initial transaction array */
  RECDATA rd;                   /* recursion data */

  assert(tabag && report);      /* check the function arguments */
  rd.mode = mode;               /* make search mode consistent */
  rd.smin = (smin > 0) ? smin : 1; /* check and adapt the support */
  if (tbg_wgt(tabag) < rd.smin) /* check the total transaction weight */
    return 0;                   /* against the minimum support */
  k = tbg_itemcnt(tabag);       /* get and check the number of items */
  if (k <= 0) return isr_report(report);
  n = tbg_cnt(tabag);           /* get the number of transactions */
  a = (TAAE*)malloc((size_t)(n+1) *sizeof(TAAE));
  if (!a) return -1;            /* create initial transaction array */
  for (i = n; --i >= 0; ) {     /* traverse the transactions */
    t = tbg_tract(tabag, i);    /* initialize a new array element */
    a[i].items = ta_items(t);   /* for each transaction and */
    a[i].occ   = ta_wgt(t);     /* set the item array (pointer) */
  }                             /* and the transaction weight */
  if (a[n-1].items[0] <= TA_END)
    n--;                        /* remove an empty transaction */
  a[n].items = NULL;            /* store a sentinel at the end */
  rd.fim16   = NULL;            /* default: no 16 items machine */
  if (mode & SAM_FIM16) {       /* if to use a 16 items machine */
    rd.fim16 = m16_create(-1, rd.smin, report);
    if (!rd.fim16) { free(a); return -1; }
  }                             /* create a 16 items machine */
  rd.report = report;           /* note the item set reporter */
  r = rec_opt(a, n, &rd);       /* execute split and merge recursion */
  if (rd.fim16)                 /* if a 16 items machine was used, */
    m16_delete(rd.fim16);       /* delete the 16 items machine */
  free(a);                      /* deallocate the transaction array */
  if (r >= 0)                   /* if no error occurred, */
    r = isr_report(report);     /* report the empty item set */
  return r;                     /* return the error status */
}  /* sam_opt() */

/*----------------------------------------------------------------------
  Split and Merge: Double Source Buffering
----------------------------------------------------------------------*/

static int rec_dsb (TAAE *a, TID n, RECDATA *rd)
{                               /* --- split and merge recursion */
  int  r = 0;                   /* error status */
  ITEM i;                       /* current item */
  TID  k, m;                    /* number of remaining transactions */
  TAAE *proj, *b;               /* projected transaction database */
  TAAE *s, *t, *d;              /* to traverse the transactions */
  TAAE *e, *f;                  /* ends of the transaction sources */
  SUPP supp;                    /* support of (current) split item */
  SUPP pex;                     /* minimum support for perfect exts. */

  assert(a && (n > 0) && rd);   /* check the function arguments */
  pex  = (rd->mode & SAM_PERFECT) ? isr_supp(rd->report) : SUPP_MAX;
  proj = (TAAE*)malloc((size_t)(n+n+2) *sizeof(TAAE));
  if (!proj) return -1;         /* allocate the projection array */
  e = a+n; f = b = proj+n+n+1;  /* get the ends of the input arrays */
  b->items = NULL;              /* init. the second source (clear it) */
  while (a->items || b->items){ /* split and merge loop */
    d = proj; s = a; t = b;     /* -- split the transaction array */
    supp = 0;                   /* init. the projection weight */
    if      (!t->items) i = s->items[0];
    else if (!s->items) i = t->items[0];
    else { i = ((s->items[0] > t->items[0]) ? s : t)->items[0]; }
    if (i < 0) {                /* if only packed items left */
      for ( ; s->items; s++)    /* add trans. to 16 items machine */
        m16_add(rd->fim16, (BITTA)(s->items[0] & ~TA_END), s->occ);
      for ( ; t->items; t++)    /* add trans. to 16 items machine */
        m16_add(rd->fim16, (BITTA)(t->items[0] & ~TA_END), t->occ);
      r = m16_mine(rd->fim16);  /* mine with 16 items machine */
      break;                    /* and abort the search loop */
    }
    while (s->items          &&  t->items
    &&    (s->items[0] == i) && (t->items[0] == i)) {
      int c = cmp(s->items+1, t->items+1);
      if      (c > 0) {         /* if transaction of first  source */
        d->items = ++s->items;  /* is larger, copy it to destination */
        supp += (d++)->occ = (s++)->occ; }
      else if (c < 0) {         /* if transaction of second source */
        d->items = ++t->items;  /* is larger, copy it to destination */
        supp += (d++)->occ = (t++)->occ; }
      else {                    /* if transactions in both sources */
        d->items = ++s->items;  /* are equal, combine them and */
        ++t->items;             /* skip the leading item in both */
        supp += (d++)->occ = (s++)->occ +(t++)->occ;
      }                         /* copy trasanction to destination */
    }                           /* and sum the transaction weights */
    while (s->items && (s->items[0] == i)) {
      d->items = ++s->items;    /* copy trans. with current item */
      supp += (d++)->occ = (s++)->occ;
    }                           /* sum the item occurrences */
    while (t->items && (t->items[0] == i)) {
      d->items = ++t->items;    /* copy trans. with current item */
      supp += (d++)->occ = (t++)->occ;
    }                           /* sum the item occurrences */
    if (supp >= pex) {          /* if item is a perfect extension, */
      isr_addpex(rd->report,i); /* add it to the item set reporter */
      if ((s > a) && ((s-1)->items[0] <= TA_END)) --s;
      s->items = NULL; e = s;   /* remove an empty transaction */
      if ((t > b) && ((t-1)->items[0] <= TA_END)) --t;
      t->items = NULL; f = t;   /* remove an empty transaction */
      continue;                 /* store a sentinel at the end, */
    }                           /* and skip the current item */
    if ((d-1)->items[0] <= TA_END)
      --d;                      /* remove an empty transaction, */
    d->items = NULL;            /* store a sentinel at the end and */
    n = (TID)(d-proj);          /* compute the number of elements */
    k = (TID)(e-s);             /* in the two source buffers */
    m = (TID)(f-t);             /* -- merge the transaction arrays */
    if (k < m) {  d = a = s-n; b = t; e = NULL; }
    else { k = m; d = b = t-n; a = s; s = t;    }
    t = proj;                   /* get the destination array */
    if (((n << 4) > k)          /* if size difference is limited */
    ||  !(rd->mode & SAM_BSEARCH)) {     /* compare transactions */
      while (s->items && t->items) {     /* from the two sources */
        int c = cmp(s->items, t->items); /* and copy the smaller */
        if (c) *d++ = (c > 0) ? *s++ : *t++;
        else { *d = *s++; (d++)->occ += (t++)->occ; }
      } }                       /* combine equal transactions */
    else {                      /* if sizes differ considerably */
      while (t->items && k) {   /* while transactions to insert */
        TID l, m, r, x = 1;     /* binary search variables */
        for (l = 0, r = k; l < r; ) {
          int c = cmp(t->items, s[m = (l+r) >> 1].items);
          if (c < 0)  l = m+1;  /* find the insertion position */
          else x = c, r = m;    /* in the remaining transactions */
        }                       /* of the other source */
        for (k -= l; --l >= 0; )/* copy transactions before */
          *d++ = *s++;          /* the insertion position */
        *d++ = *t++;            /* copy transaction to insert */
        if (!x) { (d-1)->occ += (s++)->occ; k--; }
      }                         /* if equal transaction in source, */
    }                           /* combine it with ins. transaction */
    while (t->items) *d++ = *t++;  /* copy remaining transactions */
    while (s->items) *d++ = *s++;  /* from the non-empty source */
    d->items = NULL;            /* store a sentinel at the end */
    if (e) f = d; else e = d;   /* set new end for modified source */
    if ((f-b > rd->merge)       /* if both sources are large */
    &&  (e-a > rd->merge)) {    /* (larger than given threshold) */
      s = a; d = a -= f-b;      /* -- merge the transaction sources */
      while (s->items && b->items) {        /* compare transactions */
        int c = cmp(s->items, b->items);    /* from the two sources */
        if (c) *d++ = (c > 0) ? *s++ : *b++;/* and copy the smaller */
        else { *d = *s++; (d++)->occ += (b++)->occ; }
      }                         /* combine equal transactions */
      while (b->items) *d++ = *b++; /* copy remaining transactions */
      while (s->items) *d++ = *s++; /* from the non-empty source */
      d->items = NULL; e = d;   /* store a sentinel at the end */
    }                           /* and note the new array end */
    if (supp < rd->smin)        /* if the support is too low, */
      continue;                 /* skip the recursive processing */
    r = isr_add(rd->report, i, supp);
    if (r <  0) break;          /* add current item to the reporter */
    if (r <= 0) continue;       /* check if item needs processing */
    if ((n > 0)                 /* if the projection is not empty */
    && isr_xable(rd->report,1)){/* and another item can be added */
      r = rec_dsb(proj, n, rd); /* search projection recursively */
      if (r < 0) break;         /* abort on a recursion error */
    }
    r = isr_report(rd->report); /* report the current item set */
    if (r < 0) break;           /* and check for an error */
    isr_remove(rd->report, 1);  /* remove the current item */
  }                             /* from the item set reporter */
  free(proj);                   /* deallocate the projection array */
  return r;                     /* return the error status */
}  /* rec_dsb() */

/*--------------------------------------------------------------------*/

int sam_dsb (TABAG *tabag, SUPP smin, TID merge, int mode,
             ISREPORT *report)
{                               /* --- split and merge algorithm */
  int     r = 0;                /* result of recursion */
  ITEM    k;                    /* number of items */
  TID     i, n;                 /* loop variable, number of trans. */
  TRACT   *t;                   /* to traverse the transactions */
  TAAE    *a;                   /* initial transaction array */
  RECDATA rd;                   /* recursion data */

  assert(tabag && report);      /* check the function arguments */
  rd.mode   = mode;             /* make search mode consistent */
  rd.smin   = (smin > 0) ? smin : 1; /* check and adapt the support */
  rd.merge  = merge;            /* initialize the recursion data */
  if (tbg_wgt(tabag) < rd.smin) /* check the total transaction weight */
    return 0;                   /* against the minimum support */
  k = tbg_itemcnt(tabag);       /* get and check the number of items */
  if (k <= 0) return isr_report(report);
  n = tbg_cnt(tabag);           /* get the number of transactions */
  a = (TAAE*)malloc((size_t)(n+2) *sizeof(TAAE));
  if (!a) return -1;            /* create initial transaction array */
  for (i = n; --i >= 0; ) {     /* traverse the transactions */
    t = tbg_tract(tabag, i);    /* initialize a new array element */
    a[i].items = ta_items(t);   /* for each transaction and */
    a[i].occ   = ta_wgt(t);     /* set the item array (pointer) */
  }                             /* and the transaction weight */
  if (a[n-1].items[0] <= TA_END)
    n--;                        /* remove an empty transaction */
  a[n].items = NULL;            /* store a sentinel at the end */
  rd.fim16   = NULL;            /* default: no 16 items machine */
  if (mode & SAM_FIM16) {       /* if to use a 16 items machine */
    rd.fim16 = m16_create(-1, rd.smin, report);
    if (!rd.fim16) { free(a); return -1; }
  }                             /* create a 16 items machine */
  rd.report = report;           /* note the item set reporter */
  r = rec_dsb(a, n, &rd);       /* execute split and merge recursion */
  if (rd.fim16)                 /* if a 16 items machine was used, */
    m16_delete(rd.fim16);       /* delete the 16 items machine */
  free(a);                      /* deallocate the transaction array */
  if (r >= 0)                   /* if no error occurred, */
    r = isr_report(report);     /* report the empty item set */
  return r;                     /* return the error status */
}  /* sam_dsb() */

/*----------------------------------------------------------------------
  Split and Merge: Transactions as Prefix Tree
----------------------------------------------------------------------*/

static int add (TTNODE **root, const ITEM *items, SUPP supp,MEMSYS *mem)
{                               /* --- add an item set to the tree */
  ITEM   i;                     /* buffer for an item */
  TTNODE **p;                   /* insertion position */
  TTNODE *node;                 /* to insert new nodes */

  assert(root                   /* check the function arguments */
  &&     items && (supp >= 0) && mem);
  p = root;                     /* start the search at the root node */
  while (1) {                   /* traverse the items of the set */
    if ((i = *items++) <= TA_END)  /* get the next item in the set */
      return 0;                    /* and abort if there is none */
    while (*p && ((*p)->item > i)) p = &(*p)->sibling;
    node = *p;                  /* find the item/insertion position */
    if (!node || (node->item != i)) break;
    node->supp += supp;         /* if the item does not exist, abort */
    p = &node->children;        /* else update the item set support */
  }                             /* and get the list of children */
  node = (TTNODE*)ms_alloc(mem);/* create a new prefix tree node */
  if (!node) return -1;         /* for the next item */
  node->item    = i;            /* store the current item and */
  node->supp    = supp;         /* the support of the item set */
  node->sibling = *p;           /* insert the created node */
  *p = node;                    /* into the sibling list */
  while (*items > TA_END) {     /* traverse the remaining items */
    node = node->children = (TTNODE*)ms_alloc(mem);
    if (!node) return -1;       /* create a new prefix tree node */
    node->item    = *items++;   /* create a new prefix tree node */
    node->supp    = supp;       /* store item and its support */
    node->sibling = NULL;       /* there are no siblings yet */
  }
  node->children = NULL;        /* last created node is a leaf */
  return 1;                     /* return that nodes were added */
}  /* add() */

/*--------------------------------------------------------------------*/

static TTNODE* merge (TTNODE *s1, TTNODE *s2)
{                               /* --- merge two node lists */
  TTNODE *out, **end = &out;    /* output node list and end pointer */

  if (!s1) return s2;           /* if there is only one node list, */
  if (!s2) return s1;           /* simply return the other list */
  end = &out;                   /* start the output list */
  while (1) {                   /* node list merge loop */
    if      (s1->item > s2->item) {
      *end = s1; end = &s1->sibling; s1 = *end; if (!s1) break; }
    else if (s2->item > s1->item) {
      *end = s2; end = &s2->sibling; s2 = *end; if (!s2) break; }
    else {                      /* copy nodes with singular items */
      s1->children = merge(s1->children, s2->children);
      s1->supp += s2->supp;     /* merge the children recursively */
      *end = s1; end = &s1->sibling; s1 = *end; s2 = s2->sibling;
      if (!s1 || !s2) break;    /* move node from the first source */
    }                           /* to the output and delete the node */
  }                             /* from the second source */
  *end = (s1) ? s1 : s2;        /* append the remaining nodes */
  return out;                   /* return the merged top-down tree */
}  /* merge() */

/*--------------------------------------------------------------------*/

static TTNODE* copy (TTNODE *src, MEMSYS *mem)
{                               /* --- copy a top-down tree */
  TTNODE *node, *c, *dst;       /* created copy of the node list */
  TTNODE **end = &dst;          /* end of the created copy */

  assert(src && mem);           /* check the function arguments */
  do {                          /* sibling list copying loop */
    c = src->children;          /* if there are children, copy them */
    if (c && ((c = copy(c, mem)) == COPYERR)) return COPYERR;
    *end = node = (TTNODE*)ms_alloc(mem);
    if (!node) return COPYERR;  /* create a copy of the current node */
    node->item  = src->item;    /* store the item */
    node->supp  = src->supp;    /* and the support */
    node->children = c;         /* set the (copied) children */
    end = &node->sibling;       /* get the new end of the output */
    src = src->sibling;         /* and the next sibling node */
  } while (src);                /* while there is another node */
  *end = NULL;                  /* terminate the copied list */
  return dst;                   /* and return the created copy */
}  /* copy() */

/*--------------------------------------------------------------------*/

static int rec_tree (TTNODE *node, RECDATA *rd)
{                               /* --- find item sets recursively */
  int    r = 0;                 /* error status */
  SUPP   pex;                   /* minimum for perfect extensions */
  TTNODE *proj = NULL;          /* created projection */

  assert(node && rd);           /* check the function arguments */
  pex = (rd->mode & SAM_PERFECT) ? isr_supp(rd->report) : SUPP_MAX;
  for ( ; node; node = merge(node->sibling, node->children)) {
    if (node->item < 0) {       /* if only packed items left */
      do { assert(node->item < 0);
        m16_add(rd->fim16, (BITTA)(node->item & ~TA_END),node->supp);
      } while ((node = node->sibling));    /* add transactions */
      r = m16_mine(rd->fim16); break;      /* to 16 items machine */
    }                           /* and mine with 16 items machine */
    if (node->supp < rd->smin) continue; /* skip infrequent items */
    if (node->supp >= pex) {    /* and collect perfect extensions */
      isr_addpex(rd->report, node->item); continue; }
    r = isr_add(rd->report, node->item, node->supp);
    if (r <  0) break;          /* add current item to the reporter */
    if (r <= 0) continue;       /* check if item needs processing */
    if (node->children          /* if current node has children and */
    &&  isr_xable(rd->report, 1)) {    /* another item can be added */
      if (ms_push(rd->mem) < 0) {
        r = -1; break; }        /* copy the subtree for the item */
      proj = copy(node->children, rd->mem);
      r = (proj != COPYERR) ? rec_tree(proj, rd) : -1;
      ms_pop(rd->mem);          /* recursively process the copy */
      if (r < 0) break;         /* then restore memory system state */
    }                           /* and check for an error */
    r = isr_report(rd->report); /* report the current item set */
    if (r < 0) break;           /* and check for an error */
    isr_remove(rd->report, 1);  /* remove the current item */
  }                             /* from the item set reporter */
  return r;                     /* return the error status */
}  /* rec_tree() */

/*--------------------------------------------------------------------*/

int sam_tree (TABAG *tabag, SUPP smin, int mode, ISREPORT *report)
{                               /* --- search for frequent item sets */
  int     r = 0;                /* result of recursion */
  ITEM    k;                    /* number of items */
  TID     i, n;                 /* loop variable, number of trans. */
  TRACT   *t;                   /* to traverse the transactions */
  TTNODE  *root;                /* root of transaction prefix tree */
  RECDATA rd;                   /* structure for recursive search */

  assert(tabag && report);      /* check the function arguments */
  rd.mode = mode;               /* make search mode consistent */
  rd.smin = (smin > 0) ? smin : 1; /* check and adapt the support */
  if (tbg_wgt(tabag) < rd.smin) /* check the total transaction weight */
    return 0;                   /* against the minimum support */
  k = tbg_itemcnt(tabag);       /* get and check the number of items */
  if (k <= 0) return isr_report(report);
  rd.mem = ms_create(sizeof(TTNODE), 65535);
  if (!rd.mem) return 0;        /* create a memory mgmt. system */
  rd.fim16 = NULL;              /* default: no 16 items machine */
  if (mode & SAM_FIM16) {       /* if to use a 16 items machine */
    rd.fim16 = m16_create(-1, rd.smin, report);
    if (!rd.fim16) { ms_delete(rd.mem); return -1; }
  }                             /* create a 16 items machine */
  root = NULL;                  /* initialize the root node */
  for (i = n = tbg_cnt(tabag); --i >= 0; ) {
    t = tbg_tract(tabag, i);    /* traverse the transactions */
    n = add(&root, ta_items(t), ta_wgt(t), rd.mem);
    if (n < 0) break;           /* add the transaction to the tree */
  }                             /* and check for a memory error */
  if (n >= 0) {                 /* if a frequent pattern tree */
    rd.report = report;         /* has successfully been built, */
    r = rec_tree(root, &rd);    /* note the item set reporter, */
    if (r >= 0)                 /* find freq. item sets recursively */
      r = isr_report(report);   /* finally report the empty item set */
  }
  if (rd.fim16)                 /* if a 16 items machine was used, */
    m16_delete(rd.fim16);       /* delete the 16 items machine */
  ms_delete(rd.mem);            /* delete the memory mgmt. system */
  return r;                     /* return the error status */
}  /* sam_tree() */

/*----------------------------------------------------------------------
  Split and Merge: Unlimited Item Insertions
----------------------------------------------------------------------*/

static int rec_ins (TXAE *a, TID n, ITEM k, RECDATA *rd)
{                               /* --- split and merge recursion */
  int    r = 0;                 /* function result, error status */
  TXAE   *proj;                 /* projected transaction database */
  TXAE   *s, *t, *d, *e;        /* to traverse the transactions */
  double pen, wgt;              /* insertion penalty and tra. weight */
  double sum;                   /* sum of transaction weights */
  SUPP   supp;                  /* support of current item */
  double pex;                   /* minimum weight for perfect exts. */

  assert(a && (n > 0) && rd);   /* check the function arguments */
  pex  = (rd->mode & SAM_PERFECT) ? isr_wgt(rd->report) : INFINITY;
  proj = (TXAE*)malloc((size_t)(n+1) *sizeof(TXAE));
  if (!proj) return -1;         /* allocate the projection array and */
  while (--k >= 0) {            /* split and merge loop */
    pen = ib_getpen(rd->base,k);/* get the item insertion penalty */
    d = t = (pen <= 0) ? proj : (TXAE*)rd->buf;
    sum = 0; supp = 0;          /* -- split the transaction array */
    for (s = a; s->items && (s->items[0] == k); s++, d++) {
      d->items = ++s->items;    /* copy all transactions that */
      sum  += d->wgt = s->wgt;  /* start with the current item */
      supp += d->occ = s->occ;  /* sum the transaction weights */
    }                           /* and the actual occurrences */
    if (d == proj) continue;    /* skip items that do not occur */
    d->items = NULL;            /* store a sentinel at the end */
    if (sum >= pex) {           /* identify perfect extensions */
      isr_addpex(rd->report, k); continue; }
    if (pen <= 0) {             /* if item insertion is not allowed */
      e = d; d = a;             /* -- merge the transaction arrays */
      while (s->items && t->items) {         /* compare transactions */
        int c = cmp(s->items, t->items);     /* from the two sources */
        if (c) *d++ = (c > 0) ? *s++ : *t++; /* and copy the smaller */
        else { *d = *s++; d->occ += t->occ; (d++)->wgt += (t++)->wgt; }
      }                         /* combine equal transactions */
      while (t->items) *d++ = *t++;   /* copy remaining transactions */
      while (s->items) *d++ = *s++; } /* from the non-empty source */
    else {                      /* if item insertion is allowed */
      e = proj; d = a;          /* -- merge the transaction arrays */
      while (s->items && t->items) {       /* compare transactions */
        int c = cmp(s->items, t->items);   /* from the two sources */
        if      (c < 0)         /* if trans. with item is smaller, */
          *d++ = *e++ = *t++;   /* copy it to both destinations */
        else if (c > 0) {       /* if trans. w/o  item is smaller */
          sum     += wgt = rd->tnorm(s->wgt, pen);
          e->wgt   = wgt;       /* copy it with a penalized weight */
          e->occ   = 0;         /* to the second destination and */
          e->items = s->items;  /* with full weight to the first */
          *d++ = *s++; e++; }   /* no occurrences (item is missing) */
        else {                  /* if the transactions are equal */
          sum     += wgt = rd->tnorm(s->wgt, pen);
          d->wgt   = t->wgt +s->wgt; e->wgt = t->wgt +wgt;
          d->occ   = t->occ +s->occ; e->occ = t->occ;
          d->items = e->items = s->items;
          s++; t++; d++; e++;   /* combine equal transactions, */
        }                       /* but weight them differently */
      }                         /* for the two destinations */
      while (t->items)          /* copy remaining transactions */
        *d++ = *e++ = *t++;     /* containing the current item */
      while (s->items) {        /* while transactions w/o item left */
        sum     += wgt = rd->tnorm(s->wgt, pen);
        e->wgt   = wgt;         /* copy it with a penalized weight */
        e->occ   = 0;           /* to the second destination and */
        e->items = s->items;    /* with full weight to the first */
        *d++ = *s++; e++;       /* clear the number of occurrences */
      }                         /* (the item is missing) */
    }
    d->items = e->items = NULL; /* store a sentinel at the end */
    if ((supp < rd->smin)       /* if the item set support */
    ||  (sum  < rd->sins))      /* or the weight are too low */
      continue;                 /* skip the recursive processing */
    r = isr_addwgt(rd->report, k, supp, sum);
    if (r <  0) break;          /* add current item to the reporter */
    if (r <= 0) continue;       /* and check if it needs processing */
    n = (TID)(e -proj);         /* get num. of items in projection */
    if ((n > 0)                 /* if the projection is not empty */
    && isr_xable(rd->report,1)){/* and another item can be added */
      r = rec_ins(proj, n, k, rd);
      if (r < 0) break;         /* search projection recursively */
    }
    r = isr_report(rd->report); /* report the current item set */
    if (r < 0) break;           /* and check for an error */
    isr_remove(rd->report, 1);  /* and remove the current item */
  }                             /* from the item set reporter */
  free(proj);                   /* deallocate the projection array */
  return r;                     /* return the error status */
}  /* rec_ins() */

/*--------------------------------------------------------------------*/

int sam_ins (TABAG *tabag, SUPP smin, double sins, int tnorm,
             int mode, ISREPORT *report)
{                               /* --- split and merge algorithm */
  int     r;                    /* result of recursion */
  ITEM    k;                    /* number of items */
  TID     i, n;                 /* loop variable, number of trans. */
  TRACT   *t;                   /* to traverse the transactions */
  TXAE    *a;                   /* initial transaction array */
  RECDATA rd;                   /* recursion data */

  assert(tabag && report);      /* check the function arguments */
  rd.mode   = mode;             /* check and adapt minimum support */
  rd.smin   = (smin > 0) ? smin : 0;
  rd.sins   = (sins > 0) ? sins : DBL_MIN;
  if ((tnorm < 0) || (tnorm >= (int)(sizeof(tnorms)/sizeof(*tnorms))))
    tnorm = T_MIN;              /* check and adapt the t-norm */
  rd.tnorm = tnorms[tnorm];     /* note the t-norm for item weights */
  if (tbg_wgt(tabag) < rd.smin) /* check the total transaction weight */
    return 0;                   /* against the minimum support */
  k = tbg_itemcnt(tabag);       /* get and check the number of items */
  if (k <= 0) return isr_report(report);
  n = tbg_cnt(tabag);           /* get the number of transactions */
  a = (TXAE*)malloc((size_t)(n+n+2) *sizeof(TXAE));
  if (!a) return -1;            /* create initial transaction array */
  for (i = n; --i >= 0; ) {     /* traverse the transactions */
    t = tbg_tract(tabag, i);    /* initialize a new array element */
    a[i].items = ta_items(t);   /* store the item array (pointer) */
    a[i].wgt   = (double)(a[i].occ = ta_wgt(t));
  }                             /* store the transaction weight */
  a[n].items = NULL;            /* store a sentinel at the end */
  rd.buf     = a+n+1;           /* use rest of array as a buffer */
  rd.base    = tbg_base(tabag); /* note the underlying item base */
  rd.report  = report;          /* and the item set reporter */
  r = rec_ins(a, n, k, &rd);    /* execute split and merge recursion */
  free(a);                      /* deallocate the transaction array */
  if (r >= 0)                   /* if no error occurred, */
    r = isr_report(report);     /* report the empty item set */
  return r;                     /* return the error status */
}  /* sam_ins() */

/*----------------------------------------------------------------------
  Split and Merge: Limited Item Insertions
----------------------------------------------------------------------*/

static int rec_lim (TZAE *a, TID n, ITEM k, RECDATA *rd)
{                               /* --- split and merge recursion */
  int    r = 0;                 /* function result, error status */
  TZAE   *proj;                 /* projected transaction database */
  TZAE   *s, *t, *d, *e;        /* to traverse the transactions */
  TZAE   *x, *y;                /* ditto, for equal trans. merging */
  double pen, wgt;              /* insertion penalty and tra. weight */
  double sum;                   /* sum of transaction weights */
  SUPP   supp;                  /* support of current item */
  double pex;                   /* minimum weight for perfect exts. */

  assert(a && (n > 0) && rd);   /* check the function arguments */
  pex  = (rd->mode & SAM_PERFECT) ? isr_wgt(rd->report) : INFINITY;
  proj = (TZAE*)malloc((size_t)(n+1) *sizeof(TZAE));
  if (!proj) return -1;         /* allocate the projection array and */
  while (--k >= 0) {            /* split and merge loop */
    pen = ib_getpen(rd->base,k);/* get the item insertion penalty */
    d = t = (pen <= 0) ? proj : (TZAE*)rd->buf;
    sum = 0; supp = 0;          /* -- split the transaction array */
    for (s = a; s->items && (s->items[0] == k); s++, d++) {
      d->items = ++s->items;    /* copy all transactions that */
      d->cnt   = s->cnt;        /* start with the current item */
      d->wgt   = s->wgt;        /* sum the transaction weights */
      sum     += (double)s->cnt * s->wgt;
      supp    += d->occ = s->occ;
    }                           /* sum the actual occurrences */
    if (d == proj) continue;    /* skip items that do not occur */
    d->items = NULL;            /* store a sentinel at the end */
    if (sum >= pex) {           /* identify perfect extensions */
      isr_addpex(rd->report, k); continue; }
    e = proj; d = a;            /* -- merge the transaction arrays */
    while (s->items && t->items) {      /* compare transactions */
      int c = cmp(s->items, t->items);  /* from the two sources */
      if      (c < 0)           /* if trans. with item is smaller, */
        *d++ = *e++ = *t++;     /* copy it to both destinations */
      else if (c > 0) {         /* if trans. w/o  item is smaller */
        wgt = rd->tnorm(s->wgt, pen); /* compute the penalized weight */
        if (wgt >= rd->twgt) {  /* if it is greater than the minimum, */
          e->wgt = wgt;         /* copy transaction with penalized */
          e->cnt = s->cnt;      /* weight to the second destination */
          sum   += (double)s->cnt *wgt;  /* sum transaction weights */
          e->occ = 0;           /* no occurrences (item is missing) */
          (e++)->items = s->items;
        }                       /* keep the item array unchanged */
        *d++ = *s++; }          /* copy with full weight to 1st dest. */
      else {                    /* if the transactions are equal */
        x = s; y = t;           /* - merge to projection */
        if (pen > 0) {          /* if an insertion is allowed */
          while ((x->items == s->items)
          &&     (y->items == t->items)) {
            wgt = rd->tnorm(x->wgt, pen); /* compute penalized weight */
            if (wgt < rd->twgt) { x++; continue; }
            if      (wgt < y->wgt) {
              e->wgt = wgt;     /* copy weight from 1st source */
              e->cnt = x->cnt;  /* sum the transaction weight */
              sum   += (double)x->cnt *wgt;
              e->occ = 0;      x++; }
            else if (wgt > y->wgt) {
              e->wgt = y->wgt;  /* copy weight from 2nd source */
              e->cnt = y->cnt;  /* weights are already summed */
              e->occ = y->occ; y++; }
            else {              /* if the trans. weights are equal */
              e->wgt = wgt;     /* sum the transactions */
              e->cnt = x->cnt +y->cnt;
              sum   += (double)x->cnt *wgt;
              e->occ = y->occ; x++; y++;
            }                   /* combine the two transactions */
            (e++)->items = s->items;
          }                     /* copy item array from 1st source */
          while (x->items == s->items) {
            wgt = rd->tnorm(x->wgt, pen); /* compute penalized weight */
            if (wgt < rd->twgt) { x++; continue; }
            e->items = s->items;/* copy the remaining transactions */
            e->wgt   = wgt;     /* with an unchanged item array */
            e->cnt   = x->cnt;  /* sum the transaction weights */
            sum     += (double)x->cnt *wgt;
            e->occ   = 0; e++; x++;
          }                     /* copy rest from first source) */
        }
        while (y->items == t->items) {
          e->items = s->items;  /* copy the remaining transactions */
          e->wgt   = y->wgt;    /* with an unchanged item array */
          e->cnt   = y->cnt;    /* and their full trans. weight */
          e->occ   = y->occ; e++; y++;
        }
        x = s; y = t;           /* - merge to remove item */
        while ((s->items == x->items)
        &&     (t->items == y->items)) {
          if      (s->wgt < t->wgt) {
            d->wgt = s->wgt;    /* copy weight and counters */
            d->cnt = s->cnt;    /* from the first  souce */
            d->occ = s->occ; s++; }
          else if (s->wgt > t->wgt) {
            d->wgt = t->wgt;    /* copy weight and counters */
            d->cnt = t->cnt;    /* from the second source */
            d->occ = t->occ; t++; }
          else {                /* if the weights are equal */
            d->wgt = t->wgt;    /* copy weight from second source */
            d->cnt = s->cnt +t->cnt;
            d->occ = s->occ +t->occ; s++; t++;
          }                     /* combine the two transactions */
          (d++)->items = x->items;
        }                       /* store the item array */
        while (s->items == x->items)
          *d++ = *s++;          /* copy the remaining transactions */
        while (t->items == y->items) {
          d->items = x->items;  /* copy the remaining transactions */
          d->wgt   = t->wgt;    /* from the second source */
          d->cnt   = t->cnt;
          d->occ   = t->occ; d++; t++;
        }                       /* (only for equal transactions */
      }                         /* the two merge processes are */
    }                           /* separated to allow combining) */
    while (t->items)            /* copy remaining transactions */
      *d++ = *e++ = *t++;       /* containing the current item */
    while (s->items) {          /* while trans. w/o current item left */
      wgt = rd->tnorm(s->wgt, pen);    /* get penalized trans. weight */
      if (wgt >= rd->twgt) {    /* if it is greater than the minimum, */
        e->items = s->items;    /* copy the transaction */
        e->cnt   = s->cnt;      /* with the penalized weight */
        e->wgt   = wgt;         /* to the second destination */
        sum     += (double)s->cnt *wgt; /* sum transaction weights */
        e->occ   = 0; e++;      /* clear the number of occurrences */
      }                         /* (the item is missing) */
      *d++ = *s++;              /* in any case copy it with full */
    }                           /* weight to the 1st destination */
    d->items = e->items = NULL; /* store a sentinel at the end */
    if ((supp < rd->smin)       /* if the item set support */
    ||  (sum  < rd->sins))      /* or the weight are too low, */
      continue;                 /* skip the recursive processing */
    r = isr_addwgt(rd->report, k, supp, sum);
    if (r <  0) { free(proj); return r; }
    if (r <= 0) continue;       /* add current item to the reporter */
    n = (TID)(e -proj);         /* and check if it needs processing */
    if ((n > 0)                 /* if the projection is not empty */
    && isr_xable(rd->report,1)){/* and another item can be added */
      n = rec_lim(proj,n,k,rd); /* search projection recursively */
      if (n < 0) { free(proj); return -1; }
    }
    r = isr_report(rd->report); /* report the current item set */
    if (r < 0) break;           /* and check for an error */
    isr_remove(rd->report, 1);  /* remove the current item */
  }                             /* from the item set reporter */
  free(proj);                   /* deallocate the projection array */
  return r;                     /* return the error status */
}  /* rec_lim() */

/*--------------------------------------------------------------------*/

int sam_lim (TABAG *tabag, SUPP smin, double sins, int tnorm,
             double twgt, int mode, ISREPORT *report)
{                               /* --- split and merge algorithm */
  int     r;                    /* result of recursion */
  ITEM    k;                    /* number of items */
  TID     i, n;                 /* loop variable, number of trans. */
  TRACT   *t;                   /* to traverse the transactions */
  TZAE    *a;                   /* initial transaction array */
  RECDATA rd;                   /* recursion data */

  assert(tabag && report);      /* check the function arguments */
  rd.mode   = mode;             /* check and adapt minimum support */
  rd.smin   = (smin > 0) ? smin : 0;
  rd.sins   = (sins > 0) ? sins : DBL_MIN;
  rd.twgt   = (twgt > 0) ? twgt : DBL_MIN;
  if ((tnorm < 0) || (tnorm >= (int)(sizeof(tnorms)/sizeof(*tnorms))))
    tnorm = T_MIN;              /* check and adapt the t-norm */
  rd.tnorm = tnorms[tnorm];     /* note the t-norm for item weights */
  if (tbg_wgt(tabag) < rd.smin) /* check the total transaction weight */
    return 0;                   /* against the minimum support */
  k = tbg_itemcnt(tabag);       /* get and check the number of items */
  if (k <= 0) return isr_report(report);
  n = tbg_cnt(tabag);           /* get the number of transactions */
  a = (TZAE*)malloc((size_t)(n+n+2) *sizeof(TZAE));
  if (!a) return -1;            /* create initial transaction array */
  for (i = 0; i < n; i++) {     /* traverse the transactions */
    t = tbg_tract(tabag, i);    /* initialize a new array element */
    a[i].items = ta_items(t);   /* for each transaction */
    a[i].occ   = a[i].cnt = ta_wgt(t);
    a[i].wgt   = 1.0;           /* copy the occurrence counter and */
  }                             /* init. the individual weight */
  a[n].items = NULL;            /* store a sentinel at the end */
  rd.buf     = a+n+1;           /* use rest of array as a buffer */
  rd.base    = tbg_base(tabag); /* note the underlying item base */
  rd.report  = report;          /* and the item set reporter */
  r = rec_lim(a, n, k, &rd);    /* execute split and merge recursion */
  free(a);                      /* deallocate the transaction array */
  if (r >= 0)                   /* if no error occurred, */
    r = isr_report(report);     /* report the empty item set */
  return r;                     /* return the error status */
}  /* sam_lim() */

/*----------------------------------------------------------------------
  SaM (generic)
----------------------------------------------------------------------*/

int sam_data (TABAG *tabag, int target, SUPP smin, ITEM zmin,
              double twgt, int eval, int algo, int mode, int sort)
{                               /* --- prepare data for SaM */
  ITEM    m;                    /* number of items */
  TID     n;                    /* number of transactions */
  SUPP    w;                    /* total transaction weight */
  int     pack;                 /* number of items to pack */
  clock_t t;                    /* timer for measurements */

  assert(tabag);                /* check the function arguments */
  pack = mode & SAM_FIM16;      /* get number of items to pack */
  if (pack > 16) pack = 16;     /* pack at most 16 items */

  /* --- sort and recode items --- */
  t = clock();                  /* start timer, print log message */
  XMSG(stderr, "filtering, sorting and recoding items ... ");
  m = tbg_recode(tabag, smin, -1, -1, -sort);
  if (m <  0) return E_NOMEM;   /* recode items and transactions */
  if (m <= 0) return E_NOITEMS; /* and check the number of items */
  XMSG(stderr, "[%"ITEM_FMT" item(s)]", m);
  XMSG(stderr, " done [%.2fs].\n", SEC_SINCE(t));

  /* --- sort and reduce transactions --- */
  t = clock();                  /* start timer, print log message */
  XMSG(stderr, "sorting and reducing transactions ... ");
  tbg_filter(tabag, (twgt >= 0) ? 0 : zmin, NULL, 0);
  tbg_itsort(tabag, -1, 0);     /* sort items in transactions and */
  tbg_sort  (tabag, -1, 0);     /* sort the trans. lexicographically */
  n = tbg_reduce(tabag, 0);     /* reduce transactions to unique ones */
  if ((twgt < 0) && (pack > 0)) /* if insertions and 16 items mach., */
    tbg_pack(tabag, pack);      /* pack the most frequent items */
  w = tbg_wgt(tabag);           /* get the new transaction weight */
  XMSG(stderr, "[%"TID_FMT, n); /* print number of transactions */
  if (w != (SUPP)n) { XMSG(stderr, "/%"SUPP_FMT, w); }
  XMSG(stderr, " transaction(s)] done [%.2fs].\n", SEC_SINCE(t));
  return 0;                     /* return 'ok' */
}  /* sam_data() */

/*--------------------------------------------------------------------*/

int sam_repo (ISREPORT *report, int target,
              int eval, double thresh, int algo, int mode)
{                               /* --- prepare reporter for SaM */
  assert(report);               /* check the function arguments */
  if (eval == SAM_LDRATIO)      /* set the evaluation function */
    isr_seteval(report, isr_logrto, NULL, +1, thresh);
  return (isr_settarg(report, target, 0, -1)) ? E_NOMEM : 0;
}  /* sam_repo() */

/*--------------------------------------------------------------------*/

int sam (TABAG *tabag, int target, SUPP smin, double sins,
         int tnorm, double twgt, int eval, double thresh,
         int algo, int mode, TID merge, ISREPORT *report)
{                               /* --- SaM algorithm (generic) */
  int     r;                    /* result of function call */
  clock_t t;                    /* timer for measurements */

  assert(tabag && report);      /* check the function arguments */
  t = clock();                  /* start timer, print log message */
  XMSG(stderr, "writing %s ... ", isr_name(report));
  if      (twgt >  0)           /*   limited item insertions */
    r = sam_lim (tabag, smin, sins, tnorm, twgt, mode, report);
  else if (twgt >= 0)           /* unlimited item insertions */
    r = sam_ins (tabag, smin, sins, tnorm,       mode, report);
  else if (algo == SAM_TREE)    /* transaction prefix tree */
    r = sam_tree(tabag, smin,                    mode, report);
  else if (algo == SAM_DOUBLE)  /* double source buffering */
    r = sam_dsb (tabag, smin, merge,             mode, report);
  else if (algo == SAM_BSEARCH) /* optimized with binary search merge */
    r = sam_opt (tabag, smin,                    mode, report);
  else                          /* basic frequent item set search */
    r = sam_base(tabag, smin,                    mode, report);
  if (r < 0) return E_NOMEM;    /* search for frequent item sets */
  XMSG(stderr, "[%"SIZE_FMT" set(s)]", isr_repcnt(report));
  XMSG(stderr, " done [%.2fs].\n", SEC_SINCE(t));
  return 0;                     /* return 'ok' */
}  /* sam() */

/*----------------------------------------------------------------------
  Main Functions
----------------------------------------------------------------------*/
#ifdef SAM_MAIN

static void help (void)
{                               /* --- print add. option information */
  #ifndef QUIET
  fprintf(stderr, "\n");        /* terminate startup message */
  printf("SaM algorithm variants (option -a#)\n");
  printf("  s   basic split and merge search\n");
  printf("  b   split and merge with binary search (default)\n");
  printf("  d   split and merge with double source buffering\n");
  printf("  t   split and merge with transaction prefix tree\n");
  printf("\n");
  printf("additional evaluation measures (option -e#)\n");
  printf("  x   no measure (default)\n");
  printf("  b   binary logarithm of support quotient\n");
  printf("\n");
  printf("t-norms (triangular norms) for combining item penalties"
         " (option -N#)\n");
  printf("  m   minimum              T(a,b) = min(a,b)\n");
  printf("  n   nil-potent minimum   T(a,b) = min(a,b) "
                                    "if a+b > 1 else 0\n");
  printf("  p   product              T(a,b) = a*b\n");
  printf("  l   Lukasiewicz          T(a,b) = max(0,a+b-1)\n");
  printf("  h   Hamacher product     T(a,b) = 0 if a = b = 0 "
                                    "else a*b/(a+b-a*b)\n");
  printf("\n");
  printf("information output format characters (option -v#)\n");
  printf("  %%%%  a percent sign\n");
  printf("  %%i  number of items (item set size)\n");
  printf("  %%a  absolute item set support\n");
  printf("  %%s  relative item set support as a fraction\n");
  printf("  %%S  relative item set support as a percentage\n");
  printf("  %%w  absolute support with insertions\n");
  printf("  %%r  relative support with insertions as a fraction\n");
  printf("  %%R  relative support with insertions as a percentage\n");
  printf("  %%e  additional evaluation measure\n");
  printf("  %%E  additional evaluation measure as a percentage\n");
  printf("All format characters can be preceded by the number\n");
  printf("of significant digits to be printed (at most 32 digits),\n");
  printf("even though this value is ignored for integer numbers.\n");
  #endif                        /* print help information */
  exit(0);                      /* abort the program */
}  /* help() */

/*--------------------------------------------------------------------*/

static ITEM getbdr (char *s, char **end, double **border)
{                               /* --- get the support border */
  ITEM   i, k;                  /* loop variables */
  double *b;                    /* support border */

  assert(s && end && border);   /* check the function arguments */
  for (i = k = 0; s[i]; i++)    /* traverse the string and */
    if (s[i] == ':') k++;       /* count the number separators */
  *border = b = (double*)malloc((size_t)++k *sizeof(double));
  if (!b) return -1;            /* allocate a support border */
  for (i = 0; i < k; i++) {     /* traverse the parameters */
    b[i] = strtod(s, end);      /* get the next parameter and */
    if (*end == s) break;       /* check for an empty parameter */
    s = *end; if (*s++ != ':') break;
  }                             /* check for a colon */
  if (++i < k)                  /* shrink support array if possible */
    *border = (double*)realloc(b, (size_t)i *sizeof(double));
  return i;                     /* return number of support values */
}  /* getbdr() */

/*--------------------------------------------------------------------*/

static int setbdr (ISREPORT *report, SUPP w, ITEM zmin,
                   double **border, ITEM n)
{                               /* --- set the support border */
  double s;                     /* to traverse the support values */

  assert(report                 /* check the function arguments */
  &&    (w > 0) && (zmin >= 0) && border && (*border || (n <= 0)));
  while (--n >= 0) {            /* traverse the support values */
    s = (*border)[n];           /* transform to absolute count */
    s = ceilsupp((s >= 0) ? 0.01 *s *(double)w *(1-DBL_EPSILON) : -s);
    if (isr_setbdr(report, n+zmin, (RSUPP)s) < 0) return -1;
  }                             /* set support in item set reporter */
  if (*border) { free(*border); *border = NULL; }
  return 0;                     /* return 'ok' */
}  /* setbdr() */

/*--------------------------------------------------------------------*/

#ifndef NDEBUG                  /* if debug version */
  #undef  CLEANUP               /* clean up memory and close files */
  #define CLEANUP \
  if (twrite) twr_delete(twrite, 1); \
  if (report) isr_delete(report, 0); \
  if (tabag)  tbg_delete(tabag,  0); \
  if (tread)  trd_delete(tread,  1); \
  if (ibase)  ib_delete (ibase);     \
  if (border) free(border);
#endif

GENERROR(error, exit)           /* generic error reporting function */

/*--------------------------------------------------------------------*/

int main (int argc, char *argv[])
{                               /* --- main function */
  int     i, k = 0;             /* loop variables, counters */
  char    *s;                   /* to traverse the options */
  CCHAR   **optarg = NULL;      /* option argument */
  CCHAR   *fn_inp  = NULL;      /* name of input  file */
  CCHAR   *fn_out  = NULL;      /* name of output file */
  CCHAR   *fn_sel  = NULL;      /* name of item selection file */
  CCHAR   *fn_psp  = NULL;      /* name of pattern spectrum file */
  CCHAR   *recseps = NULL;      /* record  separators */
  CCHAR   *fldseps = NULL;      /* field   separators */
  CCHAR   *blanks  = NULL;      /* blank   characters */
  CCHAR   *comment = NULL;      /* comment characters */
  CCHAR   *hdr     = "";        /* record header  for output */
  CCHAR   *sep     = " ";       /* item separator for output */
  CCHAR   *dflt    = " (%S)";   /* default format for check */
  CCHAR   *info    = dflt;      /* format for information output */
  int     target   = 's';       /* target type (closed/maximal) */
  double  supp     = 10;        /* minimum support of an item set */
  SUPP    smin     = 1;         /* minimum support of an item set */
  double  sins     = 10;        /* minimum support with insertions */
  ITEM    zmin     = 1;         /* minimum size of an item set */
  ITEM    zmax     = ITEM_MAX;  /* maximum size of an item set */
  int     tnorm    = 'p';       /* t-norm for combining item weights */
  double  twgt     = -1;        /* minimum transaction weight */
  int     eval     = 'x';       /* additional evaluation measure */
  double  thresh   = 10;        /* threshold for evaluation measure */
  int     sort     =  2;        /* flag for item sorting and recoding */
  int     algo     = 'b';       /* variant of SaM algorithm */
  int     mode     = SAM_DEFAULT;  /* search mode (e.g. pruning) */
  int     pack     = 16;        /* number of bit-packed items */
  TID     merge    = 8192;      /* threshold for source list merging */
  int     mtar     = 0;         /* mode for transaction reading */
  int     scan     = 0;         /* flag for scanable item output */
  int     bdrcnt   = 0;         /* number of support values in border */
  int     stats    = 0;         /* flag for item set statistics */
  PATSPEC *psp;                 /* collected pattern spectrum */
  ITEM    m;                    /* number of items */
  TID     n;                    /* number of transactions */
  SUPP    w;                    /* total transaction weight */
  clock_t t;                    /* timer for measurements */

  #ifndef QUIET                 /* if not quiet version */
  prgname = argv[0];            /* get program name for error msgs. */

  /* --- print usage message --- */
  if (argc > 1) {               /* if arguments are given */
    fprintf(stderr, "%s - %s\n", argv[0], DESCRIPTION);
    fprintf(stderr, VERSION); } /* print a startup message */
  else {                        /* if no arguments given */
    printf("usage: %s [options] infile [outfile]\n", argv[0]);
    printf("%s\n", DESCRIPTION);
    printf("%s\n", VERSION);
    printf("-t#      target type                              "
                    "(default: %c)\n", target);
    printf("         (s: frequent, c: closed, m: maximal item sets)\n");
    printf("-m#      minimum number of items per item set     "
                    "(default: %"ITEM_FMT")\n", zmin);
    printf("-n#      maximum number of items per item set     "
                    "(default: no limit)\n");
    printf("-s#      minimum support of an item set           "
                    "(default: %g%%)\n", supp);
    printf("         (positive: percentage, "
                     "negative: absolute number)\n");
    printf("-i#      minimum support with item insertions     "
                    "(default: %g%%)\n", sins);
    printf("         (only with item insertions, option -u)\n");
    printf("-N#      t-norm for combining item penalties      "
                    "(default: %c)\n",   tnorm);
    printf("-u#      minimum weight of a transaction          "
                    "(default: %g)\n",   twgt);
    printf("         (a value >= 0 selects item insertions)\n");
    printf("-e#      additional evaluation measure            "
                    "(default: none)\n");
    printf("-d#      threshold for add. evaluation measure    "
                    "(default: %g%%)\n", thresh);
    printf("-q#      sort items w.r.t. their frequency        "
                    "(default: %d)\n", sort);
    printf("         (1: ascending, -1: descending, 0: do not sort,\n"
           "          2: ascending, -2: descending w.r.t. "
                    "transaction size sum)\n");
    printf("-A#      variant of the SaM algorithm to use      "
                    "(default: %c)\n", algo);
    printf("-y#      threshold for transaction source merging "
                    "(default: %"TID_FMT")\n", merge);
    printf("         (for algorithm variant 'b', option '-ab')\n");
    printf("-x       do not prune with perfect extensions     "
                    "(default: prune)\n");
    printf("-l#      number of items for k-items machine      "
                    "(default: %d)\n", pack);
    printf("-F#:#..  support border for filtering item sets   "
                    "(default: none)\n");
    printf("         (list of minimum support values, "
                    "one per item set size,\n");
    printf("         starting at the minimum size, "
                    "as given with option -m#)\n");
    printf("-R#      read item selection/insertion penalties\n");
    printf("-P#      write a pattern spectrum to a file\n");
    printf("-Z       print item set statistics "
                    "(number of item sets per size)\n");
    printf("-g       write output in scanable form "
                    "(quote certain characters)\n");
    printf("-h#      record header  for output                "
                    "(default: \"%s\")\n", hdr);
    printf("-k#      item separator for output                "
                    "(default: \"%s\")\n", sep);
    printf("-v#      output format for item set information   "
                    "(default: \"%s\")\n", info);
    printf("-w       integer transaction weight in last field "
                    "(default: only items)\n");
    printf("-r#      record/transaction separators            "
                    "(default: \"\\n\")\n");
    printf("-f#      field /item        separators            "
                    "(default: \" \\t,\")\n");
    printf("-b#      blank   characters                       "
                    "(default: \" \\t\\r\")\n");
    printf("-C#      comment characters                       "
                    "(default: \"#\")\n");
    printf("-!       print additional option information\n");
    printf("infile   file to read transactions from           "
                    "[required]\n");
    printf("outfile  file to write frequent item sets to      "
                    "[optional]\n");
    return 0;                   /* print a usage message */
  }                             /* and abort the program */
  #endif  /* #ifndef QUIET */
  /* free option characters: ijlop [A-Z]\[CFPRZ] */

  /* --- evaluate arguments --- */
  for (i = 1; i < argc; i++) {  /* traverse arguments */
    s = argv[i];                /* get option argument */
    if (optarg) { *optarg = s; optarg = NULL; continue; }
    if ((*s == '-') && *++s) {  /* -- if argument is an option */
      while (*s) {              /* traverse options */
        switch (*s++) {         /* evaluate switches */
          case '!': help();                          break;
          case 't': target = (*s) ? *s++ : 's';      break;
          case 'm': zmin   = (ITEM)strtol(s, &s, 0); break;
          case 'n': zmax   = (ITEM)strtol(s, &s, 0); break;
          case 's': supp   =       strtod(s, &s);    break;
          case 'i': sins   =       strtod(s, &s);    break;
          case 'N': tnorm  = (*s) ? *s++ : 'p';      break;
          case 'u': twgt   =       strtod(s, &s);    break;
          case 'e': eval   = (*s) ? *s++ : 0;        break;
          case 'd': thresh =       strtod(s, &s);    break;
          case 'q': sort   = (int) strtol(s, &s, 0); break;
          case 'A': algo   = (*s) ? *s++ : 0;        break;
          case 'x': mode  &= ~SAM_PERFECT;           break;
          case 'l': pack   = (int) strtol(s, &s, 0); break;
          case 'y': merge  = (TID) strtol(s, &s, 0); break;
          case 'F': bdrcnt = getbdr(s, &s, &border); break;
          case 'R': optarg = &fn_sel;                break;
          case 'P': optarg = &fn_psp;                break;
          case 'Z': stats  = 1;                      break;
          case 'g': scan   = 1;                      break;
          case 'h': optarg = &hdr;                   break;
          case 'k': optarg = &sep;                   break;
          case 'v': optarg = &info;                  break;
          case 'w': mtar  |= TA_WEIGHT;              break;
          case 'r': optarg = &recseps;               break;
          case 'f': optarg = &fldseps;               break;
          case 'b': optarg = &blanks;                break;
          case 'C': optarg = &comment;               break;
          default : error(E_OPTION, *--s);           break;
        }                       /* set option variables */
        if (optarg && *s) { *optarg = s; optarg = NULL; break; }
      } }                       /* get option argument */
    else {                      /* -- if argument is no option */
      switch (k++) {            /* evaluate non-options */
        case  0: fn_inp = s;      break;
        case  1: fn_out = s;      break;
        default: error(E_ARGCNT); break;
      }                         /* note filenames */
    }
  }
  if (optarg)       error(E_OPTARG);     /* check option arguments */
  if (k      < 1)   error(E_ARGCNT);     /* and number of arguments */
  if (zmin   < 0)   error(E_SIZE, zmin); /* check the size limits */
  if (zmax   < 0)   error(E_SIZE, zmax); /* and the minimum support */
  if (supp   > 100) error(E_SUPPORT, supp);
  if (sins   > 100) error(E_SUPPORT, sins);
  if (twgt   > 1)   error(E_WEIGHT,  twgt);
  if (bdrcnt < 0)   error(E_NOMEM);
  if ((!fn_inp || !*fn_inp) && (fn_sel && !*fn_sel))
    error(E_STDIN);             /* stdin must not be used twice */
  switch (algo) {               /* check and translate alg. variant */
    case 's': algo = SAM_BASIC;              break;
    case 'b': algo = SAM_BSEARCH;            break;
    case 'd': algo = SAM_DOUBLE;             break;
    case 't': algo = SAM_TREE;               break;
    default : error(E_VARIANT, (char)algo);  break;
  }                             /* (get eclat algorithm code) */
  switch (target) {             /* check and translate target type */
    case 's': target = ISR_ALL;              break;
    case 'c': target = ISR_CLOSED;           break;
    case 'm': target = ISR_MAXIMAL;          break;
    default : error(E_TARGET, (char)target); break;
  }                             /* (get target type code) */
  switch (tnorm) {              /* check and translate t-norm */
    case 'm': tnorm = T_MIN;                 break;
    case 'n': tnorm = T_NILP;                break;
    case 'p': tnorm = T_PROD;                break;
    case 'l': tnorm = T_LUKA;                break;
    case 'h': tnorm = T_HAMA;                break;
    default : error(E_TNORM, (char)tnorm);   break;
  }                             /* (get triangular norm) */
  switch (eval) {               /* check and translate measure */
    case 'x': eval = SAM_NONE;               break;
    case 'b': eval = SAM_LDRATIO;            break;
    default : error(E_MEASURE, (char)eval);  break;
  }                             /* (get evaluation measure code) */
  if (pack > 0)                 /* add packed items to search mode */
    mode |= (pack < 16) ? pack : 16;
  if (info == dflt)             /* adapt the default info. format */
    info = (supp < 0) ? " (%a)" : " (%S)";
  thresh *= 0.01;               /* scale the evaluation threshold */
  MSG(stderr, "\n");            /* terminate the startup message */

  /* --- read item selection/insertion penalties --- */
  ibase = ib_create(0, 0);      /* create an item base */
  if (!ibase) error(E_NOMEM);   /* to manage the items */
  tread = trd_create();         /* create a transaction reader */
  if (!tread) error(E_NOMEM);   /* and configure the characters */
  trd_allchs(tread, recseps, fldseps, blanks, "", comment);
  if (fn_sel) {                 /* if an item selection is given */
    t = clock();                /* start timer, open input file */
    if (trd_open(tread, NULL, fn_sel) != 0)
      error(E_FOPEN, trd_name(tread));
    MSG(stderr, "reading %s ... ", trd_name(tread));
    m = (twgt >= 0)             /* depending on the target type */
      ? ib_readpen(ibase,tread) /* read the insertion penalties */
      : ib_readsel(ibase,tread);/* or the given item selection */
    if (m < 0) error((int)-m, ib_errmsg(ibase, NULL, 0));
    trd_close(tread);           /* close the input file */
    MSG(stderr, "[%"ITEM_FMT" item(s)]", m);
    MSG(stderr, " done [%.2fs].\n", SEC_SINCE(t));
  }                             /* print a log message */

  /* --- read transaction database --- */
  tabag = tbg_create(ibase);    /* create a transaction bag */
  if (!tabag) error(E_NOMEM);   /* to store the transactions */
  t = clock();                  /* start timer, open input file */
  if (trd_open(tread, NULL, fn_inp) != 0)
    error(E_FOPEN, trd_name(tread));
  MSG(stderr, "reading %s ... ", trd_name(tread));
  k = tbg_read(tabag, tread, mtar);
  if (k < 0) error(-k, tbg_errmsg(tabag, NULL, 0));
  trd_delete(tread, 1);         /* read the transaction database, */
  tread = NULL;                 /* then delete the table reader */
  m = ib_cnt(ibase);            /* get the number of items, */
  n = tbg_cnt(tabag);           /* the number of transactions, */
  w = tbg_wgt(tabag);           /* the total transaction weight */
  MSG(stderr, "[%"ITEM_FMT" item(s), %"TID_FMT, m, n);
  if (w != (SUPP)n) { MSG(stderr, "/%"SUPP_FMT, w); }
  MSG(stderr, " transaction(s)] done [%.2fs].", SEC_SINCE(t));
  if ((m <= 0) || (n <= 0))     /* check for at least one item */
    error(E_NOITEMS);           /* and at least one transaction */
  MSG(stderr, "\n");            /* terminate the log message */
  supp = (supp >= 0) ? 0.01 *supp *(double)w *(1-DBL_EPSILON) : -supp;
  sins = (sins >= 0) ? 0.01 *sins *(double)w *(1-DBL_EPSILON) : -sins;
  smin = (SUPP)ceilsupp(supp);  /* compute absolute support value */

  /* --- find frequent item sets --- */
  mode |= SAM_VERBOSE|SAM_NOCLEAN;
  k = sam_data(tabag, target, smin, zmin, twgt, eval, algo, mode, sort);
  if (k) error(k);              /* prepare data for SaM */
  report = isr_create(ibase);   /* create an item set reporter */
  if (!report) error(E_NOMEM);  /* and configure it */
  isr_setsize(report,        zmin, zmax);
  isr_setsupp(report, (RSUPP)smin, RSUPP_MAX);
  if (setbdr(report, w, zmin, &border, bdrcnt) != 0)
    error(E_NOMEM);             /* set limits and support border */
  if (fn_psp && (isr_addpsp(report, NULL) < 0))
    error(E_NOMEM);             /* set a pattern spectrum if req. */
  if (isr_setfmt(report, scan, hdr, sep, NULL, info) != 0)
    error(E_NOMEM);             /* set the output format strings */
  k = isr_open(report, NULL, fn_out);
  if (k) error(k, isr_name(report)); /* open the item set file */
  if ((sam_repo(report, target, eval, thresh, algo, mode) < 0)
  ||  (isr_setup(report) < 0))  /* prepare reporter for SaM and */
    error(E_NOMEM);             /* set up the item set reporter */
  k = sam(tabag, target, smin, sins, tnorm, twgt, eval, thresh,
          algo, mode, merge, report);
  if (k) error(k);              /* find frequent item sets */
  if (stats)                    /* print item set statistics */
    isr_prstats(report, stdout, 0);
  if (isr_close(report) != 0)   /* close item set output file */
    error(E_FWRITE, isr_name(report));

  /* --- write pattern spectrum --- */
  if (fn_psp) {                 /* if to write a pattern spectrum */
    t = clock();                /* start timer, create table write */
    psp    = isr_getpsp(report);/* get the pattern spectrum */
    twrite = twr_create();      /* create a table writer and */
    if (!twrite) error(E_NOMEM);/* open the output file */
    if (twr_open(twrite, NULL, fn_psp) != 0)
      error(E_FOPEN,  twr_name(twrite));
    MSG(stderr, "writing %s ... ", twr_name(twrite));
    if (psp_report(psp, twrite, 1.0) != 0)
      error(E_FWRITE, twr_name(twrite));
    twr_delete(twrite, 1);      /* write the pattern spectrum */
    twrite = NULL;              /* and delete the table writer */
    MSG(stderr, "[%"SIZE_FMT" signature(s)]", psp_sigcnt(psp));
    MSG(stderr, " done [%.2fs].\n", SEC_SINCE(t));
  }                             /* write a log message */

  /* --- clean up --- */
  CLEANUP;                      /* clean up memory and close files */
  SHOWMEM;                      /* show (final) memory usage */
  return 0;                     /* return 'ok' */
}  /* main() */

#endif
