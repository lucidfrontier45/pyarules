/*----------------------------------------------------------------------
  File    : apriori.h
  Contents: apriori algorithm for finding frequent item sets
  Author  : Christian Borgelt
  History : 2011.07.18 file created
            2011.10.18 several mode flags added
            2013.03.30 adapted to type changes in module tract
            2014.08.21 parameter 'body' added to function apriori()
            2014.08.28 functions apriori_data() and apriori_repo() added
----------------------------------------------------------------------*/
#ifndef __APRIORI__
#define __APRIORI__
#ifndef TATREEFN
#define TATREEFN
#endif
#include "istree.h"

/*----------------------------------------------------------------------
  Preprocessor Definitions
----------------------------------------------------------------------*/
/* --- evaluation measures --- */
/* measure definitions in ruleval.h */
#define APR_LDRATIO RE_FNCNT    /* binary log. of support quotient */
#define APR_INVBXS  IST_INVBXS  /* inval. eval. below exp. supp. */

/* --- aggregation modes --- */
#define APR_NONE    IST_NONE    /* no aggregation (use first value) */
#define APR_FIRST   IST_FIRST   /* no aggregation (use first value) */
#define APR_MIN     IST_MIN     /* minimum of measure values */
#define APR_MAX     IST_MAX     /* maximum of measure values */
#define APR_AVG     IST_AVG     /* average of measure values */

/* --- operation modes --- */
#define APR_PERFECT IST_PERFECT /* prune with perfect extensions */
#define APR_TATREE  0x1000      /* use transaction tree */
#define APR_POST    0x2000      /* use a-posteriori pruning */
#define APR_DEFAULT (APR_PERFECT|APR_TATREE)
#ifdef NDEBUG
#define APR_NOCLEAN 0x8000      /* do not clean up memory */
#else                           /* in function apriori() */
#define APR_NOCLEAN 0           /* in debug version */
#endif                          /* always clean up memory */
#define APR_VERBOSE INT_MIN     /* verbose message output */

/* --- algorithm variants --- */
#define APR_BASIC   0           /* basic algorithm (dummy) */

/*----------------------------------------------------------------------
  Functions
----------------------------------------------------------------------*/
extern int apriori_data (TABAG *tabag, int target, SUPP smin, ITEM zmin,
                         int eval, int algo, int mode, int sort);
extern int apriori_repo (ISREPORT *report, int target,
                         int eval, double thresh, int algo, int mode);
extern int apriori      (TABAG *tabag, int target, SUPP smin, SUPP body,
                         double conf, int eval, int agg, double thresh,
                         ITEM prune, int algo, int mode, double filter,
                         int order, ISREPORT *report);
#endif
