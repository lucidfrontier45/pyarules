/*----------------------------------------------------------------------
  File    : pfxtree.h
  Contents: prefix tree management for item sets
  Author  : Christian Borgelt
  History : 2009.10.08 file created
            2009.10.13 function pxt_isect()  added (trans. intersection)
            2009.10.26 function pxt_dir()    added (get item order)
            2009.10.30 function pxt_report() added (recursive reporting)
            2009.11.18 function pxt_prune()  added (support pruning)
            2010.02.04 reduced to functions relevant for IsTa
            2010.06.25 support-based pruning added to isect functions
            2010.08.05 update step counter added to the nodes
            2010.08.18 function pxt_nodecnt() added (number of nodes)
            2012.04.29 function pxt_super() added (check for superset)
            2013.04.01 adapted to type changes in module tract
----------------------------------------------------------------------*/
#ifndef __PFXTREE__
#define __PFXTREE__
#include "memsys.h"
#include "tract.h"
#ifndef ISR_CLOMAX
#define ISR_CLOMAX
#endif
#include "report.h"

/*----------------------------------------------------------------------
  Type Definitions
----------------------------------------------------------------------*/
typedef struct pfxnode {        /* --- a prefix tree node --- */
  ITEM           item;          /* associated item (last item in set) */
  SUPP           supp;          /* support of represented item set */
  TID            step;          /* last update step (for pxt_isect) */
  struct pfxnode *sibling;      /* successor node in sibling list */
  struct pfxnode *children;     /* list of child nodes */
} PFXNODE;                      /* (prefix tree node) */

typedef struct {                /* --- a prefix tree --- */
  MEMSYS         *mem;          /* memory management system */
  ITEM           size;          /* number of items / array size */
  int            dir;           /* direction of item order */
  TID            step;          /* last update step (for pxt_isect) */
  ITEM           last;          /* last item        (for pxt_isect) */
  SUPP           supp;          /* current support  (for pxt_isect) */
  SUPP           min;           /* minimum support   for reporting */
  ISREPORT       *rep;          /* item set reporter for reporting */
  PFXNODE        root;          /* root node of the prefix tree */
  SUPP           mins[1];       /* minimum support values */
} PFXTREE;                      /* (prefix tree) */

/*----------------------------------------------------------------------
  Functions
----------------------------------------------------------------------*/
extern PFXTREE* pxt_create  (ITEM size, int dir, MEMSYS  *mem);
extern void     pxt_delete  (PFXTREE *pxt, int delms);
extern MEMSYS*  pxt_memsys  (PFXTREE *pxt);
extern size_t   pxt_nodecnt (PFXTREE *pxt);
extern size_t   pxt_nodemax (PFXTREE *pxt);
extern int      pxt_dir     (PFXTREE *pxt);
extern SUPP     pxt_supp    (PFXTREE *pxt);

extern int      pxt_add     (PFXTREE *pxt, const ITEM *items, ITEM n,
                             SUPP supp);
extern int      pxt_isect   (PFXTREE *pxt, const ITEM *items, ITEM n,
                             SUPP supp, SUPP min, const SUPP *frqs);
extern SUPP     pxt_get     (PFXTREE *pxt, const ITEM *items, ITEM n);
extern int      pxt_super   (PFXTREE *pxt, const ITEM *items, ITEM n,
                             SUPP supp);
extern int      pxt_prunex  (PFXTREE *pxt, SUPP supp, const SUPP *frqs);
extern void     pxt_prune   (PFXTREE *pxt, SUPP supp);

extern int      pxt_report  (PFXTREE *pxt, int max, SUPP supp,
                             ISREPORT *rep);

#ifndef NDEBUG
extern void     pxt_show    (PFXTREE *pxt, ITEMBASE *base);
#endif

/*----------------------------------------------------------------------
  Preprocessor Definitions
----------------------------------------------------------------------*/
#define pxt_memsys(t)    ((t)->mem)
#define pxt_nodecnt(t)   (ms_used((t)->mem))
#define pxt_nodemax(t)   (ms_umax((t)->mem))
#define pxt_dir(t)       ((t)->dir)
#define pxt_supp(t)      ((t)->root.supp)

#endif
