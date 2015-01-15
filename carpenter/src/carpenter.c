/*----------------------------------------------------------------------
  File    : carpenter.c
  Contents: carpenter algorithm for finding frequent item sets
  Author  : Christian Borgelt
  History : 2010.06.23 file created from eclat.c
            2010.07.07 first version with tid list array running
            2010.07.08 order of the items in tid list array reversed
            2010.07.09 initial tid list construction optimized
            2010.07.13 variant based on item occurrence table added
            2010.07.14 transaction database reduction bug fixed
            2010.07.22 adapted to specialized item set repository
            2010.08.19 item selection file added as optional input
            2010.08.22 adapted to modified modules tabread and tract
            2010.10.15 adapted to modified interface of module report
            2010.11.24 adapted to modified error reporting (tract)
            2010.12.11 adapted to a generic error reporting function
            2010.12.20 adapted to function tbg_icnts() (filter problem)
            2011.07.08 adapted to modified function tbg_recode()
            2011.08.28 output of item set counters per size added
            2011.11.23 output/filtering of maximal item sets added
            2012.04.27 alternative filtering of maximal item sets added
            2012.05.18 optional transaction weights and option -p added
            2012.06.13 bug resulting from exchange of n and k fixed
            2013.04.02 adapted to type changes in module tract
            2013.10.18 optional pattern spectrum collection added
            2013.11.12 item selection file changed to option -R#
            2014.05.12 option -F# added (support border for filtering)
            2014.08.02 option -c renamed to -z (maximal item set filter)
            2014.08.28 functions carp_data() and carp_repo() added
            2014.10.24 changed from LGPL license to MIT license
------------------------------------------------------------------------
  Reference for the Carpenter algorithm:
    F. Pan, G. Cong, A.K.H. Tung, J. Yang, and M. Zaki.
    Carpenter: Finding Closed Patterns in Long Biological Datasets.
    Proc. 9th ACM SIGKDD Int. Conf. on Knowledge Discovery
    and Data Mining (KDD 2003, Washington, DC), 637-642.
    ACM Press, New York, NY, USA 2003
  Reference for some improvements:
    C. Borgelt, X. Yang, R. Nogales-Cadenas,
    P. Carmona-Saez, and A. Pascual-Montano.
    Finding Closed Frequent Item Sets by Intersecting Transactions.
    Proc. 14th Int. Conf. on Extending Database Technology
    (EDBT 2011, Uppsala, Sweden), 367-376.
    ACM Press, New York, NY, USA 2011
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
#ifdef CARP_MAIN
#ifndef PSP_REPORT
#define PSP_REPORT
#endif
#ifndef TA_READ
#define TA_READ
#endif
#endif
#include "carpenter.h"
#ifdef CARP_MAIN
#include "error.h"
#endif
#ifdef STORAGE
#include "storage.h"
#endif

/*----------------------------------------------------------------------
  Preprocessor Definitions
----------------------------------------------------------------------*/
#define PRGNAME     "carpenter"
#define DESCRIPTION "find closed/maximal frequent item sets " \
                    "with the carpenter algorithm"
#define VERSION     "version 3.11 (2014.10.24)        " \
                    "(c) 2010-2014   Christian Borgelt"

/* --- error codes --- */
/* error codes   0 to  -4 defined in tract.h */
#define E_STDIN      (-5)       /* double assignment of stdin */
#define E_OPTION     (-6)       /* unknown option */
#define E_OPTARG     (-7)       /* missing option argument */
#define E_ARGCNT     (-8)       /* too few/many arguments */
#define E_TARGET     (-9)       /* invalid target type */
#define E_SIZE      (-10)       /* invalid item set size */
#define E_SUPPORT   (-11)       /* invalid item set support */
#define E_VARIANT   (-12)       /* invalid algorithm variant */
#define E_MEASURE   (-13)       /* invalid evaluation measure */
/* error codes -15 to -25 defined in tract.h */

#ifndef QUIET                   /* if not quiet version, */
#define MSG         fprintf     /* print messages */
#define XMSG        if (mode & CARP_VERBOSE) fprintf
#else                           /* if quiet version, */
#define MSG(...)                /* suppress messages */
#define XMSG(...)
#endif

#define SEC_SINCE(t)  ((double)(clock()-(t)) /(double)CLOCKS_PER_SEC)

/*----------------------------------------------------------------------
  Type Definitions
----------------------------------------------------------------------*/
typedef struct {                /* --- trans. identifier list --- */
  ITEM     item;                /* item identifier */
  SUPP     supp;                /* item support (number of trans.) */
  TID      *tids;               /* transaction identifiers */
} TIDLIST;                      /* (transaction identifier list) */

typedef struct {                /* --- recursion data --- */
  int      mode;                /* operation mode (e.g. pruning) */
  SUPP     smin;                /* minimum support of an item set */
  ITEM     zmin;                /* minimum size    of an item set */
  SUPP     **tab;               /* item occurrence counter table */
  SUPP     *muls;               /* multiplicity of transactions */
  ITEM     *set;                /* buffer for an item set */
  REPOTREE *rpt;                /* repository of item sets */
} RECDATA;                      /* (recursion data) */

/*----------------------------------------------------------------------
  Constants
----------------------------------------------------------------------*/
#if !defined QUIET && defined CARP_MAIN
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
  /* E_VARIANT -12 */  "invalid carpenter variant '%c'",
  /* E_MEASURE -13 */  "invalid evaluation measure '%c'",
  /*           -14 */  NULL,
  /* E_NOITEMS -15 */  "no (frequent) items found",
  /*           -16 */  "unknown error"
};
#endif

/*----------------------------------------------------------------------
  Global Variables
----------------------------------------------------------------------*/
#ifdef CARP_MAIN
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
  Auxiliary Functions for Debugging
----------------------------------------------------------------------*/
#if !defined NDEBUG && defined CARP_MAIN

static void indent (int k)
{ while (--k >= 0) printf("   "); }

/*--------------------------------------------------------------------*/

static void show_tab (const char *text, SUPP **tab, ITEM k, TID n)
{                               /* --- show item counter table */
  ITEM i;                       /* loop variable */
  TID  r;                       /* loop variable */

  printf("%s\n", text);         /* print the given text */
  printf("  ");                 /* skip row id/tid column */
  for (i = 0; i < k; i++)       /* print the item header */
    printf(" %s/%"ITEM_FMT, ib_name(ibase, i), i);
  printf("\n");                 /* terminate the header line */
  for (r = 0; r < n; r++) {     /* traverse the table rows */
    printf("%"TID_FMT":", r);   /* print the row number / tid */
    for (i = 0; i < k; i++) printf(" %3"SUPP_FMT, tab[r][i]);
    printf("\n");               /* print the item counters */
  }                             /* and terminate the output line */
}  /* show_tab() */

/*--------------------------------------------------------------------*/

static void show_spt (const char *text, ITEM **items, SUPP **tab, TID n)
{                               /* --- show item counter table */
  TID  r;                       /* loop variable */
  ITEM *s;                      /* to traverse items */
  SUPP *c;                      /* to traverse counters */

  printf("%s\n", text);         /* print the given text */
  for (r = 0; r < n; r++) {     /* traverse the table rows */
    printf("%"TID_FMT":", r);   /* print the row number / tid */
    for (s = items[r], c = tab[r]; *s >= 0; s++)
      printf(" %s/%"ITEM_FMT":%"SUPP_FMT, ib_name(ibase, *s), *s, *c++);
    printf("\n");               /* print the item information */
  }                             /* and terminate the output line */
}  /* show_spt() */

/*--------------------------------------------------------------------*/

static void show_set (const char *text, ITEM *set, ITEM n, int ind)
{                               /* --- show an item set */
  ITEM k, *s;                   /* to traverse the items */

  indent(ind);                  /* indent the output line */
  printf("%s:", text);          /* print the given text */
  for (k = 0, s = set; *s >= 0; s++) {
    printf(" %s/%"ITEM_FMT, ib_name(ibase, *s), *s); k++; }
  assert(k == n);               /* print the list of items */
  printf(" (%"ITEM_FMT")\n",n); /* and check their number */
}  /* show_set() */

/*--------------------------------------------------------------------*/

static void show_tid (const char *text, TIDLIST *lists, ITEM k, int ind)
{                               /* --- show a cond. trans. database */
  TID n, *s;                    /* to traverse the transaction ids */

  indent(ind);                  /* indent the output line */
  printf("%s\n", text);         /* print the given text */
  for ( ; --k >= 0; lists++) {  /* traverse the item lists */
    indent(ind);                /* indent the output line */
    printf("%s/%"ITEM_FMT":", ib_name(ibase, lists->item), lists->item);
    for (n = 0, s = lists->tids; *s >= 0; s++) {
      printf(" %"TID_FMT, *s); n++; } /* print the item and the tids */
    printf(" (%"SUPP_FMT")\n", lists->supp);
  }                             /* print the item support */
}  /* show_tid() */

#endif  /* #ifndef NDEBUG */
/*----------------------------------------------------------------------
  Carpenter based on an Item Occurrence Table
----------------------------------------------------------------------*/

static SUPP rec_tab (ITEM *set, ITEM k, TID n, SUPP supp, RECDATA *rd)
{                               /* --- carpenter recursion */
  ITEM i, m;                    /* loop variables, item counter */
  SUPP s, r;                    /* needed support, error status */
  SUPP *row;                    /* to traverse the table rows */
  ITEM *dst;                    /* intersection destination */
  ITEM pex;                     /* minimum items for perfext exts. */

  assert(set                    /* check the function arguments */
  &&    (k > 0) && (n >= 0) && (supp >= 0) && rd);
  dst = set +k;                 /* get destination for intersections */
  pex = (rd->mode & CARP_PERFECT) ? k : ITEM_MAX;
  s   = rd->smin -supp -1;      /* get minimum for perfect exts. */
  if (s < 0) s = 0;             /* and end value for trans. loop */
  while (--n >= s) {            /* traverse remaining transactions */
    row = rd->tab[n];           /* filter item set with table row */
    for (i = m = 0; i < k; i++) /* corresp. to current transaction */
      if (row[set[i]] > s) dst[m++] = set[i];
    if (m <  rd->zmin) continue;/* skip too small intersections */
    if (m <= 1) {               /* if there is only one item left */
      r = (SUPP)rpt_add(rd->rpt, dst, m, supp +row[*dst]);
      if (r < 0) return r; continue;
    }                           /* update the item set repository */
    if (m >= pex) {             /* collect perfect extensions */
      ++supp; if (s > 0) --s; continue; }
    if ((rd->mode & CARP_MAXONLY) /* if to find maximal item sets */
    &&  rpt_super(rd->rpt, dst, m, rd->smin))
      continue;                 /* check for a frequent superset */
    r = (SUPP)rpt_add(rd->rpt, dst, m, supp+1);
    if (r <  0) return r;       /* add item set to the repository */
    if (r <= 0) continue;       /* find closed item sets recursively */
    r = rec_tab(dst, m, n, supp+1, rd);
    if (r > supp+1)             /* if there are perfect extensions */
      r = (SUPP)rpt_add(rd->rpt, dst, m, r);
    if (r < 0) return r;        /* update the item set repository */
  }                             /* and finally check for an error */
  return supp;                  /* return the item set support */
}  /* rec_tab() */

/*--------------------------------------------------------------------*/

static SUPP rec_mtb (ITEM *set, ITEM k, TID n, SUPP supp, RECDATA *rd)
{                               /* --- carpenter recursion */
  ITEM i, m;                    /* loop variables, item counter */
  SUPP s, r;                    /* needed support, error status */
  SUPP *row;                    /* to traverse the table rows */
  ITEM *dst;                    /* intersection destination */
  ITEM pex;                     /* minimum items for perfext exts. */

  assert(set                    /* check the function arguments */
  &&    (k > 0) && (n >= 0) && (supp >= 0) && rd);
  dst = set +k;                 /* get destination for intersections */
  pex = (rd->mode & CARP_PERFECT) ? k : ITEM_MAX;
  while (--n >= 0) {            /* traverse the transaction ids */
    s = rd->smin -supp -1;      /* compute the minimum support -1 */
    if (s < 0) s = 0;           /* needed in remaining transactions */
    row = rd->tab[n];           /* filter item set with table row */
    for (i = m = 0; i < k; i++) /* corresp. to current transaction */
      if (row[set[i]] > s) dst[m++] = set[i];
    if (m <  rd->zmin) continue;/* skip too small intersections */
    if (m <= 1) {               /* if there is only one item left */
      r = (SUPP)rpt_add(rd->rpt, dst, m, supp +row[*dst]);
      if (r < 0) return r; continue;
    }                           /* update the item set repository */
    if (m >= pex) {             /* collect perfect extensions */
      supp += rd->muls[n]; continue; }
    if ((rd->mode & CARP_MAXONLY) /* if to find maximal item sets */
    &&  rpt_super(rd->rpt, dst, m, rd->smin))
      continue;                 /* check for a frequent superset */
    s = supp +rd->muls[n];      /* compute support with transaction */
    r = (SUPP)rpt_add(rd->rpt, dst, m, s);
    if (r <  0) return r;       /* add item set to the repository */
    if (r <= 0) continue;       /* find closed item sets recursively */
    r = rec_mtb(dst, m, n, s, rd);
    if (r > s)                  /* if there are perfect extensions */
      r = (SUPP)rpt_add(rd->rpt, dst, m, r);
    if (r < 0) return r;        /* update the item set repository */
  }                             /* and finally check for an error */
  return supp;                  /* return the item set support */
}  /* rec_mtb() */

/*----------------------------------------------------------------------
Note that no memory is allocated directly in the above functions; all
processing is done in the single memory block that is allocated in the
function below. The size of this memory block is O(n*k), where n is the
number of items and k the number of transactions. Additional memory is
only allocated in rpt_add() for the closed item set repository rd->rpt.
----------------------------------------------------------------------*/

int carp_tab (TABAG *tabag, SUPP smin, ITEM zmin, int mode,
              REPOTREE *rpt)
{                               /* --- search for frequent item sets */
  ITEM       i, k;              /* loop variable, number of items */
  TID        j, n, m;           /* loop variable, number of trans. */
  size_t     x, z;              /* number of item instances */
  SUPP       w;                 /* buffer for support */
  TRACT      *t;                /* to traverse transactions */
  ITEM       *set;              /* to traverse items */
  SUPP       *row, *frqs;       /* to traverse table rows */
  const ITEM *p;                /* to traverse transaction items */
  RECDATA    rd;                /* recursion data */

  assert(tabag && rpt);         /* check the function arguments */
  rd.mode = mode;               /* note the search mode (pruning) */
  rd.smin = (smin > 0) ? smin : 1; /* check and adapt the support */
  rd.zmin = (zmin > 0) ? zmin : 1; /* and the minimum item set size */
  if ((tbg_wgt(tabag) < rd.smin)   /* check total transaction weight */
  ||  (tbg_max(tabag) < rd.zmin))  /* and maximum transaction size */
    return 0;                   /* against the minimum support */
  k = tbg_itemcnt(tabag);       /* get the number of items */
  n = tbg_cnt(tabag);           /* and the number of transactions */
  x = tbg_extent(tabag);        /* and the number of item instances */
  w = tbg_wgt(tabag);           /* as well as the total weight */
  rpt_add(rpt, NULL, 0, w);     /* add empty set to the repository */
  if (k <= 0) return 0;         /* check whether there are items */
  for (m = j = 0; j < n; j++)   /* traverse the transactions */
    if (ta_wgt(tbg_tract(tabag, j)) != 1) {
      m = n; break; }           /* check for unit weights only */
  z = (size_t)n*(size_t)k;      /* compute the table size */
  rd.tab = (SUPP**)malloc((size_t) n       *sizeof(SUPP*)
                        +((size_t)(m+k)+z) *sizeof(SUPP)
                        +((size_t) k   +x) *sizeof(ITEM));
  if (!rd.tab) return -1;       /* allocate memory for occ. table */
  rd.muls = (SUPP*)(rd.tab +n); /* split off multiplicity array, */
  frqs    = rd.muls +m;         /* the transaction counter array, */
  row     = frqs    +k;         /* the item counter table, */
  set     = (ITEM*)(row +z);    /* and the initial set buffer */
  memset(frqs, 0, ((size_t)k +z) *sizeof(SUPP));
  if (m <= 0) {                 /* if no multiplicity array needed */
    for (j = 0; j < n; j++) {   /* traverse the transactions and */
      rd.tab[j] = row;          /* organize/set the table rows */
      for (p = ta_items(tbg_tract(tabag, j)); *p >= 0; p++)
        row[*p] = frqs[*p] += 1;/* update the base counters and */
      row += k;                 /* set the table row counters, */
    } }                         /* then get the next table row */
  else {                        /* if to use a multiplicty array */
    for (j = 0; j < n; j++) {   /* traverse the transactions */
      rd.tab [j] = row;         /* organize/set the table rows */
      rd.muls[j] = w = ta_wgt(t = tbg_tract(tabag, j));
      for (p = ta_items(t); *p >= 0; p++)
        row[*p] = frqs[*p] += w;/* update the base counters and */
      row += k;                 /* set the table row counters, */
    }                           /* then get the next table row */
  }
  if (rpt_dir(rpt) > 0)         /* set the initial (full) item set */
       for (i = 0; i < k; i++) set[i] = i;
  else for (i = 0; i < k; i++) set[i] = k-1-i;
  rd.rpt = rpt;                 /* search for closed freq. item sets */
  w = (m) ? rec_mtb(set, k, n, 0, &rd)
          : rec_tab(set, k, n, 0, &rd);
  if (w > 0)                    /* if there are perfect extensions, */
    rpt_add(rpt, set, k, w);    /* update the item set repository */
  free(rd.tab);                 /* delete the allocated table/array */
  return (w < 0) ? (int)w : 0;  /* return the error status */
}  /* carp_tab() */

/*----------------------------------------------------------------------
  Carpenter based on Transaction Identifier Lists
----------------------------------------------------------------------*/

static SUPP rec_tid (TIDLIST *lists, ITEM k, TID n, SUPP supp,
                     RECDATA *rd)
{                               /* --- carpenter recursion */
  ITEM    i, m;                 /* loop variables, item counter */
  SUPP    s, r;                 /* needed support, error status */
  TIDLIST *dst;                 /* tid lists of projected database */
  ITEM    pex;                  /* minimum items for perfext exts. */

  assert(lists                  /* check the function arguments */
  &&    (k > 0) && (n >= 0) && (supp >= 0) && rd);
  dst = lists +k;               /* get destination for intersections */
  pex = (rd->mode & CARP_PERFECT) ? k : ITEM_MAX;
  s   = rd->smin -supp -1;      /* compute end value for trans. loop */
  if (s < 0) s = 0;             /* (based on min. item set support) */
  while (--n >= s) {            /* traverse the transaction ids */
    for (i = m = 0; i < k; i++) {
      if (*lists[i].tids == n){ /* traverse items in current trans. */
        ++lists[i].tids;        /* remove the current transaction id */
        if (--lists[i].supp >= s) dst[m++] = lists[i];
      }                         /* collect the frequent items in */
    }                           /* the conditional tid list array */
    if (m <  rd->zmin) continue;/* skip too small intersections */
    if (m <= 1) {               /* if there is only one item left */
      r = (SUPP)rpt_add(rd->rpt, &dst->item, 1, supp+1 +dst->supp);
      if (r < 0) return r; continue;
    }                           /* update the item set repository */
    if (m >= pex) {             /* collect perfect extensions */
      supp++; if (s > 0) s--; continue; }
    for (i = 0; i < m; i++)     /* collect the items */
      rd->set[i] = dst[i].item; /* in the current set */
    if ((rd->mode & CARP_MAXONLY) /* if to find maximal item sets */
    &&  rpt_super(rd->rpt, rd->set, m, rd->smin))
      continue;                 /* check for a frequent superset */
    r = (SUPP)rpt_add(rd->rpt, rd->set, m, supp+1);
    if (r <  0) return r;       /* add item set to the repository */
    if (r <= 0) continue;       /* check whether recursion is needed */
    r = rec_tid(dst, m, n, supp+1, rd);
    if (r > supp+1) {           /* find closed item sets recursively */
      for (i = 0; i < m; i++)   /* if there are perfect extensions */
        rd->set[i] = dst[i].item;
      r = (SUPP)rpt_add(rd->rpt, rd->set, m, r);
    }                           /* update the item set repository */
    if (r < 0) return r;        /* check for an error */
  }
  return supp;                  /* return the item set support */
}  /* rec_tid() */

/*--------------------------------------------------------------------*/

static SUPP rec_mti (TIDLIST *lists, ITEM k, TID n, SUPP supp,
                     RECDATA *rd)
{                               /* --- carpenter recursion */
  ITEM    i, m;                 /* loop variables, item counter */
  SUPP    s, r;                 /* needed support, error status */
  TIDLIST *dst;                 /* tid lists of projected database */
  ITEM    pex;                  /* minimum items for perfext exts. */

  assert(lists                  /* check the function arguments */
  &&    (k > 0) && (n >= 0) && (supp >= 0) && rd);
  dst = lists +k;               /* get destination for intersections */
  pex = (rd->mode & CARP_PERFECT) ? k : ITEM_MAX;
  while (--n >= 0) {            /* traverse the transaction ids */
    s = rd->smin -supp -rd->muls[n]; /* compute the minimum support */
    if (s < 0) s = 0;           /* needed in remaining transactions */
    for (i = m = 0; i < k; i++) {
      if (*lists[i].tids == n){ /* traverse items in current trans. */
        lists[i].tids += 1;     /* remove the current transaction id */
        lists[i].supp -= rd->muls[n];    /* and the corresp. support */
        if (lists[i].supp >= s) dst[m++] = lists[i];
      }                         /* collect the frequent items in */
    }                           /* the conditional tid list array */
    if (m <  rd->zmin) continue;/* skip too small intersections */
    if (m <= 1) {               /* if there is only one item left */
      s = supp +rd->muls[n] +dst->supp;
      r = (SUPP)rpt_add(rd->rpt, &dst->item, 1, s);
      if (r < 0) return r; continue;
    }                           /* update the item set repository */
    if (m >= pex) {             /* collect perfect extensions */
      supp += rd->muls[n]; continue; }
    for (i = 0; i < m; i++)     /* collect the items */
      rd->set[i] = dst[i].item; /* in the current set */
    if ((rd->mode & CARP_MAXONLY)/* if to find maximal item sets */
    &&  rpt_super(rd->rpt, rd->set, m, rd->smin))
      continue;                 /* check for a frequent superset */
    s = supp +rd->muls[n];      /* compute support with transaction */
    r = (SUPP)rpt_add(rd->rpt, rd->set, m, s);
    if (r <  0) return r;       /* add item set to the repository */
    if (r <= 0) continue;       /* check whether recursion is needed */
    r = rec_mti(dst, m, n, s, rd);
    if (r > s) {                /* find closed item sets recursively */
      for (i = 0; i < m; i++)   /* if there are perfect extensions */
        rd->set[i] = dst[i].item;
      r = (SUPP)rpt_add(rd->rpt, rd->set, m, r);
    }                           /* update the item set repository */
    if (r < 0) return r;        /* check for an error */
  }
  return supp;                  /* return the item set support */
}  /* rec_mti() */

/*----------------------------------------------------------------------
Note that no memory is allocated directly in the above functions; all
processing is done in the single memory block that is allocated in the
function below. The size of this memory block is O(n*k), where n is the
number of items and k the number of transactions. Additional memory is
only allocated in rpt_add() for the closed item set repository rd->rpt.
----------------------------------------------------------------------*/

int carp_tid (TABAG *tabag, SUPP smin, ITEM zmin, int mode,
              REPOTREE *rpt)
{                               /* --- search for frequent item sets */
  ITEM       i, k;              /* loop variable, number of items */
  TID        j, n, m;           /* loop variable, number of trans. */
  size_t     x;                 /* number of item instances */
  SUPP       w;                 /* buffer for support */
  int        dir;               /* direction of item order */
  TRACT      *t;                /* to traverse transactions */
  TIDLIST    *lists, *l;        /* to traverse the tid lists */
  TID        *p, **next;        /* to traverse transaction ids */
  const ITEM *s;                /* to traverse transaction items */
  const TID  *c;                /* item occurrence counters */
  RECDATA    rd;                /* recursion data */

  assert(tabag && rpt);         /* check the function arguments */
  rd.mode = mode;               /* note the search mode (pruning) */
  rd.smin = (smin > 0) ? smin : 1; /* check and adapt the support */
  rd.zmin = (zmin > 0) ? zmin : 1; /* and the minimum item set size */
  if ((tbg_wgt(tabag) < rd.smin)   /* check total transaction weight */
  ||  (tbg_max(tabag) < rd.zmin))  /* and maximum transaction size */
    return 0;                   /* against the minimum support */
  k = tbg_itemcnt(tabag);       /* get the number of items */
  n = tbg_cnt(tabag);           /* and the number of transactions */
  x = tbg_extent(tabag);        /* and number of item instances */
  w = tbg_wgt(tabag);           /* as well as the total weight */
  rpt_add(rpt, NULL, 0, w);     /* add empty set to the repository */
  if (k <= 0) return 0;         /* check whether there are items */
  for (m = j = 0; j < n; j++)   /* traverse the transactions */
    if (ta_wgt(tbg_tract(tabag, j)) != 1) {
      m = n; break; }           /* check for unit weights only */
  c = tbg_icnts(tabag, 0);      /* get the number of containing */
  if (!c) return -1;            /* transactions per item */
  lists = (TIDLIST*)malloc(((size_t)k +x) *sizeof(TIDLIST)
                          + (size_t)m     *sizeof(SUPP)
                          + (size_t)k     *sizeof(ITEM)
                          + (size_t)k     *sizeof(TID*)
                          +((size_t)k +x) *sizeof(TID));
  if (!lists) return -1;        /* create initial tid list array and */
  rd.muls = (SUPP*)(lists+k+x); /* split off multiplicity array, */
  next    = (TID**)(rd.muls+m); /* the next position array, */
  rd.set  = (ITEM*)(next   +k); /* the item set buffer, and */
  p       = (TID*) (rd.set +k); /* the transaction id arrays */
  dir     = rpt_dir(rpt);       /* get the item order direction */
  for (i = 0; i < k; i++) {     /* traverse the items/tid lists */
    l = lists +((dir < 0) ? k-1-i : i);
    l->item = i;                /* initialize the list item */
    l->supp = 0;                /* and the support counter */
    l->tids = next[i] = p;      /* store array for the trans. ids */
    p += c[i]; *p++ = (TID)-1;  /* store a sentinel at the end */
  }
  for (j = n; --j >= 0; ) {     /* traverse the transactions */
    t = tbg_tract(tabag, j);    /* get the next transaction */
    w = ta_wgt(t);              /* and its weight to store */
    if (m > 0) rd.muls[j] = w;  /* in the multiplicity array */
    for (s = ta_items(tbg_tract(tabag, j)); *s >= 0; s++) {
      lists[(dir < 0) ? k-1-*s : *s].supp += w;
      *next[*s]++ = j;          /* traverse the transaction's items */
    }                           /* sum the transaction weight and */
  }                             /* collect the transaction ids */
  rd.rpt = rpt;                 /* search for closed freq. item sets */
  w = (m) ? rec_mti(lists, k, n, 0, &rd)
          : rec_tid(lists, k, n, 0, &rd);
  if (w > 0)                    /* if there are perfect extensions, */
    rpt_add(rpt, rd.set, k, w); /* update the item set repository */
  free(lists);                  /* delete the allocated arrays */
  return (w < 0) ? (int)w : 0;  /* return the error status */
}  /* carp_tid() */

/*----------------------------------------------------------------------
  Carpenter (generic)
----------------------------------------------------------------------*/

int carp_data (TABAG *tabag, int target, SUPP smin, ITEM zmin,
               int eval, int algo, int mode, int sort)
{                               /* --- prepare data for Carpenter */
  ITEM    m;                    /* number of items */
  TID     n;                    /* number of transactions */
  SUPP    w;                    /* total transaction weight */
  clock_t t;                    /* timer for measurements */

  assert(tabag);                /* check the function arguments */
  /* --- sort and recode items --- */
  t = clock();                  /* start timer, print log message */
  XMSG(stderr, "filtering, sorting and recoding items ... ");
  m = tbg_recode(tabag, smin, -1, -1, -sort);
  if (m <  0) return E_NOMEM;   /* recode items and transactions */
  if (m <= 0) return E_NOITEMS; /* and check the number of items */
  XMSG(stderr, "[%"ITEM_FMT" item(s)]", m);
  XMSG(stderr, " done [%.2fs].\n", SEC_SINCE(t));

  /* --- filter and sort transactions --- */
  t = clock();                  /* start timer, print log message */
  XMSG(stderr, "filtering and sorting transactions ... ");
  tbg_filter(tabag,zmin,NULL,0);/* remove items of short transactions */
  tbg_itsort(tabag, -1, 0);     /* sort items in transactions and */
  tbg_sortsz(tabag, -1, 0);     /* transactions lexicographically */
  /* The sorting direction is inverted here, because the transaction */
  /* identifiers are traversed backwards in both algorithm variants. */
  if (mode & CARP_COLLATE)      /* if to collate equal transactions, */
    tbg_reduce(tabag, 0);       /* reduce transactions to unique ones */
  n = tbg_cnt(tabag);           /* get the number of transactions */
  w = tbg_wgt(tabag);           /* and the transaction weight */
  XMSG(stderr, "[%"TID_FMT, n); /* print number of transactions */
  if (w != (SUPP)n) { XMSG(stderr, "/%"SUPP_FMT, w); }
  XMSG(stderr, " transaction(s)] done [%.2fs].\n", SEC_SINCE(t));
  return 0;                     /* return 'ok' */
}  /* carp_data() */

/*--------------------------------------------------------------------*/

int carp_repo (ISREPORT *report, int target,
               int eval, double thresh, int algo, int mode)
{                               /* --- prepare reporter for Carpenter */
  int mrep;                     /* mode for item set reporter */

  assert(report);               /* check the function arguments */
  if      (target & ISR_RULES)   target = ISR_RULES;
  else if (target & ISR_GENERAS) target = ISR_GENERAS;
  else if (target & ISR_MAXIMAL) target = ISR_MAXIMAL;
  else if (target & ISR_CLOSED)  target = ISR_CLOSED;
  else                           target = ISR_ALL;
  if (eval == CARP_LDRATIO)     /* set additional evaluation function */
    isr_seteval(report, isr_logrto, NULL, +1, thresh);
  mrep = ((target & ISR_MAXIMAL) && !(mode & CARP_FILTER))
       ? ISR_MAXIMAL|ISR_MAXONLY : 0; /* build reporting mode */
  return (isr_settarg(report, target, mrep, -1)) ? E_NOMEM : 0;
}  /* carp_repo() */

/*--------------------------------------------------------------------*/

int carpenter (TABAG *tabag, int target, SUPP smin,
               int eval, double thresh, int algo, int mode,
               ISREPORT *report)
{                               /* --- run carpenter algorithm */
  int      r;                   /* result of carpenter algorithm */
  ITEM     m;                   /* number of items */
  TID      n;                   /* number of transactions */
  ITEM     zmin;                /* minimum size of an item set */
  REPOTREE *rpt;                /* repository of item sets */
  clock_t  t;                   /* timer for measurements */

  assert(tabag && report);      /* check the function arguments */

  /* --- make parameters consistent --- */
  if      (target & ISR_RULES)   target = ISR_RULES;
  else if (target & ISR_GENERAS) target = ISR_GENERAS;
  else if (target & ISR_MAXIMAL) target = ISR_MAXIMAL;
  else if (target & ISR_CLOSED)  target = ISR_CLOSED;
  else                           target = ISR_ALL;
  if (mode & CARP_MAXONLY) mode |= CARP_PERFECT;
  m = tbg_itemcnt(tabag);       /* get the number of items */
  n = tbg_cnt(tabag);           /* and the number of transactions */
  if (algo == CARP_AUTO)        /* choose the algorithm variant */
    algo = ((double)m*(double)n > 1024*1024.0)
         ? CARP_TIDLIST : CARP_TABLE;

  /* --- intersect transactions --- */
  t = clock();                  /* start timer, print log message */
  XMSG(stderr, "enumerating transaction sets ... ");
  rpt = rpt_create(NULL, m, -1);
  if (!rpt) return E_NOMEM;     /* create an item set repository */
  zmin = isr_zmin(report);      /* get the minimum item set size */
  if (algo == CARP_TIDLIST)     /* transaction identifier lists */
    r = carp_tid(tabag, smin, zmin, mode, rpt);
  else                          /* item occurrence counter table */
    r = carp_tab(tabag, smin, zmin, mode, rpt);
  if (r < 0) {                  /* intersect transaction sets */
    if (!(mode & CARP_NOCLEAN)) rpt_delete(rpt, 1); return E_NOMEM; }
  XMSG(stderr, "[%"SIZE_FMT" node(s)]", rpt_nodecnt(rpt));
  XMSG(stderr, " done [%.2fs].\n", SEC_SINCE(t));

  /* --- report found item sets --- */
  t = clock();                  /* start the timer */
  XMSG(stderr, "writing %s ... ", isr_name(report));
  r = (target & ISR_MAXIMAL) ? 1 : 0;
  if (mode & CARP_FILTER)       /* get the maximal item set mode */
    r = -r;                     /* (filter with repo. or reporter) */
  if (r < 0)                    /* if to filter with repository, */
    rpt_prune(rpt, smin);       /* prune item set repository */
  if (rpt_report(rpt, r, smin, report) < 0) { /* report item sets */
    if (!(mode & CARP_NOCLEAN)) rpt_delete(rpt, 1); return E_NOMEM; }
  XMSG(stderr, "[%"SIZE_FMT" set(s)]", isr_repcnt(report));
  XMSG(stderr, " done [%.2fs].\n", SEC_SINCE(t));
  if (!(mode & CARP_NOCLEAN)) rpt_delete(rpt, 1);
  return 0;                     /* return 'ok' */
}  /* carpenter() */

/*----------------------------------------------------------------------
  Main Functions
----------------------------------------------------------------------*/
#ifdef CARP_MAIN

static void help (void)
{                               /* --- print add. option information */
  #ifndef QUIET
  fprintf(stderr, "\n");        /* terminate startup message */
  printf("carpenter algorithm variants (option -A#)\n");
  printf("  a   automatic choice based on table size (default)\n");
  printf("  t   item occurrence counter table\n");
  printf("  l   transaction identifier lists\n");
  printf("\n");
  printf("additional evaluation measures (option -e#)\n");
  printf("  x   no measure (default)\n");
  printf("  b   binary logarithm of support quotient\n");
  printf("\n");
  printf("information output format characters (option -v#)\n");
  printf("  %%%%  a percent sign\n");
  printf("  %%i  number of items (item set size)\n");
  printf("  %%a  absolute item set support\n");
  printf("  %%s  relative item set support as a fraction\n");
  printf("  %%S  relative item set support as a percentage\n");
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
  int     target   = 'c';       /* target type (closed/maximal) */
  double  supp     = 10;        /* minimum support (in percent) */
  SUPP    smin     = 1;         /* minimum support of an item set */
  ITEM    zmin     = 1;         /* minimum size of an item set */
  ITEM    zmax     = ITEM_MAX;  /* maximum size of an item set */
  int     eval     = 'x';       /* additional evaluation measure */
  double  thresh   = 10;        /* threshold for evaluation measure */
  int     sort     = -2;        /* flag for item sorting and recoding */
  int     algo     = 'a';       /* variant of carpenter algorithm */
  int     mode     = CARP_DEFAULT;  /* search mode (e.g. pruning) */
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
    printf("         (c: closed item sets, m: maximal item sets)\n");
    printf("-m#      minimum number of items per item set     "
                    "(default: %"ITEM_FMT")\n", zmin);
    printf("-n#      maximum number of items per item set     "
                    "(default: no limit)\n");
    printf("-s#      minimum support of an item set           "
                    "(default: %g%%)\n", supp);
    printf("         (positive: percentage, "
                     "negative: absolute number)\n");
    printf("-e#      additional evaluation measure            "
                    "(default: none)\n");
    printf("-d#      threshold for add. evaluation measure    "
                    "(default: %g%%)\n", thresh);
    printf("-q#      sort items w.r.t. their frequency        "
                    "(default: %d)\n", sort);
    printf("         (1: ascending, -1: descending, 0: do not sort,\n"
           "          2: ascending, -2: descending w.r.t. "
                    "transaction size sum)\n");
    printf("-p       do not collate equal transactions        "
                    "(default: collate)\n");
    printf("-A#      variant of the carpenter algorithm       "
                    "(default: auto)\n");
    printf("-x       do not prune with perfect extensions     "
                    "(default: prune)\n");
    printf("-z       filter maximal item sets with repository "
                    "(default: extra)\n");
    printf("-y       add only maximal item sets to repository "
                    "(default: all closed)\n");
    printf("         (options -z and -y need less memory, "
                    "but are usually slower)\n");
    printf("-F#:#..  support border for filtering item sets   "
                    "(default: none)\n");
    printf("         (list of minimum support values, "
                    "one per item set size,\n");
    printf("         starting at the minimum size, "
                    "as given with option -m#)\n");
    printf("-R#      read an item selection from a file\n");
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
  /* free option characters: cijlotuyz [A-Z]\[ACFPRZ] */

  /* --- evaluate arguments --- */
  for (i = 1; i < argc; i++) {  /* traverse arguments */
    s = argv[i];                /* get option argument */
    if (optarg) { *optarg = s; optarg = NULL; continue; }
    if ((*s == '-') && *++s) {  /* -- if argument is an option */
      while (*s) {              /* traverse options */
        switch (*s++) {         /* evaluate switches */
          case '!': help();                          break;
          case 't': target = (*s) ? *s++ : 'c';      break;
          case 'm': zmin   = (ITEM)strtol(s, &s, 0); break;
          case 'n': zmax   = (ITEM)strtol(s, &s, 0); break;
          case 's': supp   =       strtod(s, &s);    break;
          case 'e': eval   = (*s) ? *s++ : 0;        break;
          case 'd': thresh =       strtod(s, &s);    break;
          case 'q': sort   = (int) strtol(s, &s, 0); break;
          case 'p': mode  &= ~CARP_COLLATE;          break;
          case 'A': algo   = (*s) ? *s++ : 0;        break;
          case 'x': mode  &= ~CARP_PERFECT;          break;
          case 'z': mode  |=  CARP_FILTER;           break;
          case 'y': mode  |=  CARP_MAXONLY;          break;
          case 'F': bdrcnt = getbdr(s, &s, &border); break;
          case 'R': optarg = &fn_sel;                break;
          case 'P': optarg = &fn_psp;                break;
          case 'Z': stats  = -1;                     break;
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
  if (bdrcnt < 0)   error(E_NOMEM);
  if ((!fn_inp || !*fn_inp) && (fn_sel && !*fn_sel))
    error(E_STDIN);             /* stdin must not be used twice */
  switch (target) {             /* check and translate target type */
    case 'c': target = ISR_CLOSED;           break;
    case 'm': target = ISR_MAXIMAL;          break;
    default : error(E_TARGET, (char)target); break;
  }                             /* (get target type code) */
  switch (algo) {               /* check and translate alg. variant */
    case 'a': algo = CARP_AUTO;               break;
    case 't': algo = CARP_TABLE;              break;
    case 'l': algo = CARP_TIDLIST;            break;
    default : error(E_VARIANT, (char)algo);  break;
  }                             /* (get eclat algorithm code) */
  switch (eval) {               /* check and translate measure */
    case 'x': eval = CARP_NONE;              break;
    case 'b': eval = CARP_LDRATIO;           break;
    default : error(E_MEASURE, (char)eval);  break;
  }
  if (info == dflt)             /* adapt the default info. format */
    info = (supp < 0) ? " (%a)" : " (%S)";
  thresh *= 0.01;               /* scale the evaluation threshold */
  MSG(stderr, "\n");            /* terminate the startup message */

  /* --- read item selection --- */
  ibase = ib_create(0, 0);      /* create an item base */
  if (!ibase) error(E_NOMEM);   /* to manage the items */
  tread = trd_create();         /* create a transaction reader */
  if (!tread) error(E_NOMEM);   /* and configure the characters */
  trd_allchs(tread, recseps, fldseps, blanks, "", comment);
  if (fn_sel) {                 /* if item appearances are given */
    t = clock();                /* start timer, open input file */
    if (trd_open(tread, NULL, fn_sel) != 0)
      error(E_FOPEN, trd_name(tread));
    MSG(stderr, "reading %s ... ", trd_name(tread));
    m = ib_readsel(ibase,tread);/* read the given item selection */
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
  smin = (SUPP)ceilsupp(supp);  /* compute absolute support value */

  /* --- find closed/maximal item sets --- */
  mode |= CARP_VERBOSE|CARP_NOCLEAN;
  k = carp_data(tabag, target, smin, zmin, eval, algo, mode, sort);
  if (k) error(k);              /* prepare data for Carpenter */
  report = isr_create(ibase);   /* create an item set reporter */
  if (!report) error(E_NOMEM);  /* and configure it */
  isr_setsize(report,        zmin, zmax);
  isr_setsupp(report, (RSUPP)smin, RSUPP_MAX);
  if (setbdr(report, w, zmin, &border, bdrcnt) != 0)
    error(E_NOMEM);             /* set the support border */
  if (fn_psp && (isr_addpsp(report, NULL) < 0))
    error(E_NOMEM);             /* add a pattern spectrum if req. */
  if (isr_setfmt(report, scan, hdr, sep, NULL, info) != 0)
    error(E_NOMEM);             /* set the output format strings */
  k = isr_open(report, NULL, fn_out);
  if (k) error(k, isr_name(report)); /* open item set file */
  if ((carp_repo(report, target, eval, thresh, algo, mode) < 0)
  ||  (isr_setup(report) < 0))  /* prepare reporter for Carpenter */
    error(E_NOMEM);             /* and set up the item set reporter */
  k = carpenter(tabag, target, smin, eval, thresh, algo, mode, report);
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
