/*----------------------------------------------------------------------
  File    : relim.h
  Contents: relim algorithm for finding frequent item sets
  Author  : Christian Borgelt
  History : 2013.11.20 file created
            2014.08.22 interface of function relim() changed
            2014.08.28 functions relim_data() and relim_repo() added
----------------------------------------------------------------------*/
#ifndef __RELIM__
#define __RELIM__
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
#define REM_NONE    0           /* no measure */
#define REM_LDRATIO 1           /* binary log. of support quotient */

/* --- algorithm variants --- */
#define REM_BASIC   0           /* basic algorithm (dummy) */

/* --- operation modes --- */
#define REM_FIM16   0x001f      /* use 16 items machine (bit rep.) */
#define REM_PERFECT 0x0020      /* prune with perfect extensions */
#define REM_DEFAULT REM_PERFECT
#ifdef NDEBUG
#define REM_NOCLEAN 0x8000      /* do not clean up memory */
#else                           /* in function relim() */
#define REM_NOCLEAN 0           /* in debug version */
#endif                          /* always clean up memory */
#define REM_VERBOSE INT_MIN     /* verbose message output */

/*----------------------------------------------------------------------
  Functions
----------------------------------------------------------------------*/
extern int relim_data (TABAG *tabag, int target, SUPP smin,
                       ITEM zmin, double twgt, int eval,
                       int algo, int mode, int sort);
extern int relim_repo (ISREPORT *report, int target,
                       int eval, double thresh, int algo, int mode);
extern int relim      (TABAG *tabag, int target, SUPP smin, double sins,
                       int tnorm, double twgt, int eval, double thresh,
                       int algo, int mode, TID merge, ISREPORT *report);
#endif
