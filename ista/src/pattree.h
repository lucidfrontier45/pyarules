/*----------------------------------------------------------------------
  File    : pattree.h
  Contents: patricia tree management for item sets
  Author  : Christian Borgelt
  History : 2012.07.06 file created from pfxtree.h
            2012.07.09 first version completed (without pruning yet)
            2012.07.13 pruning added (implementation of pat_prunex())
            2013.04.01 adapted to type changes in module tract
----------------------------------------------------------------------*/
#ifndef __PATTREE__
#define __PATTREE__
#include "tract.h"
#ifndef ISR_CLOMAX
#define ISR_CLOMAX
#endif
#include "report.h"

/*----------------------------------------------------------------------
  Type Definitions
----------------------------------------------------------------------*/
typedef struct patnode {        /* --- a patricia tree node --- */
  TID            step;          /* last update step (for pat_isect) */
  SUPP           supp;          /* support of represented item set(s) */
  struct patnode *sibling;      /* successor node in sibling list */
  struct patnode *children;     /* list of child nodes */
  ITEM           cnt;           /* number of associated items */
  ITEM           items[1];      /* associated items */
} PATNODE;                      /* (patricia tree node) */

typedef struct {                /* --- a prefix tree --- */
  ITEM           size;          /* number of items / array size */
  size_t         cnt;           /* current number of nodes */
  size_t         max;           /* maximum number of nodes */
  int            dir;           /* direction of item order */
  TID            step;          /* last update step (for pat_isect) */
  ITEM           last;          /* last item        (for pat_isect) */
  SUPP           supp;          /* current support  (for pat_isect) */
  SUPP           min;           /* minimum support   for reporting */
  int            err;           /* error indicator  (for pat_prunex) */
  ITEM           *items;        /* item buffer for intersection */
  ISREPORT       *rep;          /* item set reporter for reporting */
  PATNODE        root;          /* root node of the prefix tree */
  SUPP           mins[1];       /* minimum support values */
} PATTREE;                      /* (prefix tree) */

/*----------------------------------------------------------------------
  Functions
----------------------------------------------------------------------*/
extern PATTREE* pat_create  (ITEM size, int dir);
extern void     pat_delete  (PATTREE *pat);
extern size_t   pat_nodecnt (PATTREE *pat);
extern size_t   pat_nodemax (PATTREE *pat);
extern int      pat_dir     (PATTREE *pat);
extern SUPP     pat_supp    (PATTREE *pat);

extern int      pat_add     (PATTREE *pat, const ITEM *items, ITEM n,
                             SUPP supp);
extern int      pat_isect   (PATTREE *pat, const ITEM *items, ITEM n,
                             SUPP supp, SUPP min, const SUPP *frqs);
extern SUPP     pat_get     (PATTREE *pat, const ITEM *items, ITEM n);
extern int      pat_super   (PATTREE *pat, const ITEM *items, ITEM n,
                             SUPP supp);
extern int      pat_prunex  (PATTREE *pat, SUPP supp, const SUPP *frqs);
extern void     pat_prune   (PATTREE *pat, SUPP supp);

extern int      pat_report  (PATTREE *pat, int max, SUPP supp,
                             ISREPORT *rep);

#ifndef NDEBUG
extern void     pat_show    (PATTREE *pat, ITEMBASE *base);
#endif

/*----------------------------------------------------------------------
  Preprocessor Definitions
----------------------------------------------------------------------*/
#define pat_nodecnt(t)   ((t)->cnt)
#define pat_nodemax(t)   ((t)->max)
#define pat_dir(t)       ((t)->dir)
#define pat_supp(t)      ((t)->root.supp)

#endif
