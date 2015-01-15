/*----------------------------------------------------------------------
  File    : accretion.h
  Contents: accretion algorithm for identifying neural assemblies
  Author  : Christian Borgelt
  History : 2011.07.14 file created
            2011.07.22 adapted to new module ruleval (rule evaluation)
            2013.01.31 definition of flag ACC_INVBXS added
            2013.03.29 adapted to type changes in module tract
            2014.08.28 functions eclat_data() and eclat_repo() added
----------------------------------------------------------------------*/
#ifndef __ACCRETION__
#define __ACCRETION__
#include "tract.h"
#include "report.h"
#include "ruleval.h"

/*----------------------------------------------------------------------
  Preprocessor Definitions
----------------------------------------------------------------------*/
/* --- evaluation measures --- */
/* most definitions in ruleval.h */
#define ACC_INVBXS  INT_MIN     /* invalidate stat. below exp. supp. */

/* --- operation modes --- */
#define ACC_DEFAULT 0           /* default operation mode */
#ifdef NDEBUG
#define ACC_NOCLEAN 0x8000      /* do not clean up memory */
#else                           /* in function accretion() */
#define ACC_NOCLEAN 0           /* in debug version */
#endif                          /* always clean up memory */
#define ACC_VERBOSE INT_MIN     /* verbose message output */

/*----------------------------------------------------------------------
  Functions
----------------------------------------------------------------------*/
extern int acc_data  (TABAG *tabag, int target, SUPP smin, ITEM zmin,
                      int mode, int sort);
extern int acc_repo  (ISREPORT *report, int target, int mode);
extern int accretion (TABAG *tabag, int target, SUPP smin,
                      int stat, double siglvl,
                      int mode, ITEM maxext, ISREPORT *isrep);
#endif
