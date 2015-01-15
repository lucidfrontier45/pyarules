/*----------------------------------------------------------------------
  File    : carpenter.h
  Contents: carpenter algorithm for finding frequent item sets
  Author  : Christian Borgelt
  History : 2013.10.31 file created from carpenter.c
            2014.08.23 interface of function carpenter() changed
            2014.08.28 functions carp_data() and carp_repo() added
----------------------------------------------------------------------*/
#ifndef __CARPENTER__
#define __CARPENTER__
#include "tract.h"
#include "repotree.h"
#ifndef ISR_CLOMAX
#define ISR_CLOMAX
#endif
#include "report.h"

/*----------------------------------------------------------------------
  Preprocessor Definitions
----------------------------------------------------------------------*/
/* --- evaluation measures --- */
#define CARP_NONE     0         /* no measure */
#define CARP_LDRATIO  1         /* binary log. of support quotient */

/* --- operation modes --- */
#define CARP_PERFECT  0x0010    /* prune with perfect extensions */
#define CARP_FILTER   0x0020    /* filter maximal sets with repo. */
#define CARP_MAXONLY  0x0040    /* add only maximal sets to repo. */
#define CARP_COLLATE  0x0080    /* flag for collating transactions */
#define CARP_DEFAULT  CARP_COLLATE|CARP_PERFECT
#ifdef NDEBUG
#define CARP_NOCLEAN  0x8000
#else                           /* do not clean up memory */
#define CARP_NOCLEAN  0         /* in function carpenter() */
#endif
#define CARP_VERBOSE  INT_MIN   /* verbose message output */

/* --- carpenter variants --- */
#define CARP_AUTO     0         /* auto. choice based on table size */
#define CARP_TABLE    1         /* item occurrence counter table */
#define CARP_TIDLIST  2         /* transaction identifier lists */

/*----------------------------------------------------------------------
  Functions
----------------------------------------------------------------------*/
extern int carp_data (TABAG *tabag, int target, SUPP smin, ITEM zmin,
                      int eval, int algo, int mode, int sort);
extern int carp_repo (ISREPORT *report, int target,
                      int eval, double thresh, int algo, int mode);
extern int carpenter (TABAG *tabag, int target, SUPP smin,
                      int eval, double thresh, int algo, int mode,
                      ISREPORT *report);
#endif
