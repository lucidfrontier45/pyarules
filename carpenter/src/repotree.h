/*----------------------------------------------------------------------
  File    : repotree.h
  Contents: item set repository tree management
  Author  : Christian Borgelt
  History : 2009.10.08 file created as clomax.h
            2010.06.21 generalized by introducing definition of SUPP
            2010.07.04 rpt_add() reports whether tree was changed
            2010.07.05 function rpt_report() added (closed item sets)
            2010.07.22 file and functions specialized from clomax.h
            2010.08.18 function rpt_nodecnt() added (number of nodes)
            2012.04.26 special maximal item set functions added
            2012.04.27 function rpt_prune() added (support pruning)
            2013.04.01 adapted to type changes in module tract
----------------------------------------------------------------------*/
#ifndef __REPOTREE__
#define __REPOTREE__
#include "memsys.h"
#include "tract.h"
#ifndef ISR_CLOMAX
#define ISR_CLOMAX
#endif
#include "report.h"

/*----------------------------------------------------------------------
  Type Definitions
----------------------------------------------------------------------*/
typedef struct reponode {       /* --- repository tree node --- */
  ITEM            item;         /* associated item (last item in set) */
  SUPP            supp;         /* support of represented item set */
  struct reponode *sibling;     /* successor node in sibling list */
  struct reponode *children;    /* list of child nodes */
} REPONODE;                     /* (repository tree node) */

typedef struct {                /* --- item set repository tree --- */
  MEMSYS          *mem;         /* memory management system */
  ITEM            size;         /* (maximum) number of items */
  int             dir;          /* direction of item order */
  SUPP            supp;         /* support of the empty set */
  SUPP            min;          /* minimum support   for reporting */
  ISREPORT        *rep;         /* item set reporter for reporting */
  REPONODE        tops[1];      /* top level nodes (like roots) */
} REPOTREE;                     /* (item set repository tree) */

/*----------------------------------------------------------------------
  Item Set Repository Tree Functions
----------------------------------------------------------------------*/
REPOTREE* rpt_create  (MEMSYS *mem, ITEM size, int dir);
void      rpt_delete  (REPOTREE *rpt, int delms);
MEMSYS*   rpt_memsys  (REPOTREE *rpt);
size_t    rpt_nodecnt (REPOTREE *rpt);
size_t    rpt_nodemax (REPOTREE *rpt);
int       rpt_dir     (REPOTREE *rpt);
SUPP      rpt_supp    (REPOTREE *rpt);

int       rpt_add     (REPOTREE *rpt, const ITEM *items, ITEM n,
                       SUPP supp);
SUPP      rpt_get     (REPOTREE *rpt, const ITEM *items, ITEM n);
int       rpt_super   (REPOTREE *rpt, const ITEM *items, ITEM n,
                       SUPP min);
void      rpt_prune   (REPOTREE *rpt, SUPP supp);

int       rpt_report  (REPOTREE *rpt, int max, SUPP supp,
                       ISREPORT *rep);

#ifndef NDEBUG
void      rpt_show    (REPOTREE *rpt, ITEMBASE *base);
#endif

/*----------------------------------------------------------------------
  Preprocessor Definitions
----------------------------------------------------------------------*/
#define rpt_memsys(t)    ((t)->mem)
#define rpt_nodecnt(t)   (ms_used((t)->mem) +(size_t)(t)->size)
#define rpt_nodemax(t)   (ms_umax((t)->mem) +(size_t)(t)->size)
#define rpt_dir(t)       ((t)->dir)
#define rpt_supp(t)      ((t)->supp)

#endif
