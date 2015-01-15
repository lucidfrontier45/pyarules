/*----------------------------------------------------------------------
  File    : sam.h
  Contents: sam algorithm for finding frequent item sets
  Author  : Christian Borgelt
  History : 2013.11.20 file created
            2014.08.22 interface of function sam() changed
            2014.08.28 functions sam_data() and sam_repo() added
----------------------------------------------------------------------*/
#ifndef __SAM__
#define __SAM__
#include "tract.h"
#ifndef ISR_CLOMAX
#define ISR_CLOMAX
#endif
#include "report.h"

/*----------------------------------------------------------------------
  Preprocessor Definitions
----------------------------------------------------------------------*/
/* --- triangular norms (t-norms) --- */
#define T_MIN       0           /* minimum */
#define T_NILP      1           /* nil-potent minimum */
#define T_PROD      2           /* product */
#define T_LUKA      3           /* Lukasiewicz */
#define T_HAMA      4           /* Hamacher product */

/* --- evaluation measures --- */
#define SAM_NONE    0           /* no measure */
#define SAM_LDRATIO 1           /* binary log. of support quotient */

/* --- sam variants --- */
#define SAM_BASIC   0           /* basic split and merge */
#define SAM_BSEARCH 1           /* allow binary search merge */
#define SAM_DOUBLE  2           /* use double source buffering */
#define SAM_TREE    3           /* transaction prefix trees */

/* --- operation modes --- */
#define SAM_FIM16   0x001f      /* use 16 items machine (bit rep.) */
#define SAM_PERFECT 0x0020      /* prune with perfect extensions */
#define SAM_DEFAULT SAM_PERFECT
#ifdef NDEBUG
#define SAM_NOCLEAN 0x8000      /* do not clean up memory */
#else                           /* in function sam() */
#define SAM_NOCLEAN 0           /* in debug version */
#endif                          /* always clean up memory */
#define SAM_VERBOSE INT_MIN     /* verbose message output */

/*----------------------------------------------------------------------
  Functions
----------------------------------------------------------------------*/
extern int sam_data (TABAG *tabag, int target, SUPP smin,
                     ITEM zmin, double twgt, int eval,
                     int algo, int mode, int sort);
extern int sam_repo (ISREPORT *report, int target,
                     int eval, double thresh, int algo, int mode);
extern int sam      (TABAG *tabag, int target, SUPP smin, double sins,
                     int tnorm, double twgt, int eval, double thresh,
                     int algo, int mode, TID merge, ISREPORT *report);
#endif
