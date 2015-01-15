/*----------------------------------------------------------------------
  File    : ista.h
  Contents: finding frequent item sets by intersecting transactions
  Author  : Christian Borgelt
  History : 2014.08.24 file created
            2014.08.28 functions ista_data() and ista_repo() added
----------------------------------------------------------------------*/
#ifndef __ISTA__
#define __ISTA__
#include "tract.h"
#ifndef ISR_CLOMAX
#define ISR_CLOMAX
#endif
#include "report.h"

/*----------------------------------------------------------------------
  Preprocessor Definitions
----------------------------------------------------------------------*/
/* --- evaluation measures --- */
#define ISTA_NONE     0         /* no measure */
#define ISTA_LDRATIO  1         /* binary log. of support quotient */

/* --- operation modes --- */
#define ISTA_PRUNE    0x0010    /* prune the prefix/patricia tree */
#define ISTA_FILTER   0x0020    /* filter maximal sets with repo. */
#define ISTA_MAXONLY  0x0040    /* add only maximal sets to repo. */
#define ISTA_DEFAULT  ISTA_PRUNE
#ifdef NDEBUG
#define ISTA_NOCLEAN  0x8000    /* do not clean up memory */
#else                           /* in function ista() */
#define ISTA_NOCLEAN  0         /* in debug version */
#endif                          /* always clean up memory */
#define ISTA_VERBOSE  INT_MIN   /* verbose message output */

/* --- algorithm variant --- */
#define ISTA_PREFIX   0         /* use a prefix   tree */
#define ISTA_PATRICIA 1         /* use a patricia tree */

/*----------------------------------------------------------------------
  Functions
----------------------------------------------------------------------*/
extern int ista_data (TABAG *tabag, int target, SUPP smin, ITEM zmin,
                      int eval, int algo, int mode, int sort);
extern int ista_repo (ISREPORT *report, int target,
                      int eval, double thresh, int algo, int mode);
extern int ista      (TABAG *tabag, int target, SUPP smin,
                      int eval, double thresh, int algo, int mode,
                      ISREPORT *report);
#endif
