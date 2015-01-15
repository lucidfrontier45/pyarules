/*----------------------------------------------------------------------
  File    : pyfim.c
  Contents: Frequent Item set Mining for Python
  Author  : Christian Borgelt
  History : 2011.07.13 file created
            2011.07.19 interface for apriori implementation added
            2011.07.22 adapted to new module ruleval (rule evaluation)
            2011.07.25 adapted to modified apriori() interface
            2011.07.28 translation of test statistic strings etc. added
            2011.08.03 subset filtering options added to apriacc
            2012.04.17 functions py_eclat() and py_fpgrowth() added
            2012.08.02 missing Py_DECREFs with Append/BuildValue added
            2013.01.16 bug in function apriori() fixed (isr_create)
            2013.01.24 transactions generalized to any sequence type
            2013.02.10 module initialization adapted to Python 3
            2013.02.11 items generalized to any hashable object
            2013.02.25 transactions generalized to any iterable object
            2013.03.06 dictionary allowed for transaction database
            2013.04.01 adapted to type changes in module tract
            2013.05.08 memory leak in tbg_fromPyObj() fixed (dict)
            2013.10.18 optional pattern spectrum collection added
            2013.10.31 carpenter and supporting functions added
            2013.11.21 sam and relim and supporting functions added
            2013.12.07 returned pattern spectrum as list or dictionary
            2014.05.06 pattern spectrum generation and estimation added
            2014.05.12 support border added as fim function argument
            2014.05.15 item set evaluations "cprob" and "import" added
            2014.05.27 bugs concerning new argument 'border' fixed
            2014.06.12 bug in function py_patspec() fixed (NULL init.)
            2014.08.01 minimum improvement of evaluation measure removed
            2014.08.24 adapted to modified item set reporter interface
            2014.08.25 ista algorithm and supporting functions added
            2014.08.28 adapted to new FIM function interfaces
            2014.09.19 function py_arules() and mode parameters added
            2014.10.02 rules as target added to apriori/eclat/fpgrowth
            2014.10.09 functions repinit() and repterm() added
            2014.10.15 bug in function fim() fixed (call to eclat fn.)
            2014.10.17 bug in function patspec() fixed (rem. FPG_FIM16)
            2014.10.24 changed from LGPL license to MIT license
----------------------------------------------------------------------*/
#include <assert.h>
#if defined _WIN32 && !defined HAVE_ROUND
#define HAVE_ROUND
#endif
#include <Python.h>
#include <float.h>
#include <signal.h>
#ifdef _WIN32
#include <windows.h>
#else
#include <unistd.h>
#include <pthread.h>
#endif
#ifndef TATREEFN
#define TATREEFN
#endif
#ifndef TA_SURR
#define TA_SURR
#endif
#ifndef ISR_PATSPEC
#define ISR_PATSPEC
#endif
#ifndef ISR_CLOMAX
#define ISR_CLOMAX
#endif
#ifndef PSP_ESTIM
#define PSP_ESTIM
#endif
#include "report.h"
#include "apriori.h"
#include "eclat.h"
#include "fpgrowth.h"
#include "sam.h"
#include "relim.h"
#include "carpenter.h"
#include "ista.h"
#include "accretion.h"
#ifdef STORAGE
#include "storage.h"
#endif

/*----------------------------------------------------------------------
  Preprocessor Definitions
----------------------------------------------------------------------*/
#define ERR_VALUE(s) PyErr_SetString(PyExc_ValueError, s);  return NULL
#define ERR_TYPE(s)  PyErr_SetString(PyExc_TypeError,  s);  return NULL
#define ERR_MEM()    PyErr_SetString(PyExc_MemoryError,""); return NULL

#if PY_MAJOR_VERSION >= 3
#define PyInt_Check     PyLong_Check
#define PyInt_AsLong    PyLong_AsLong
#define PyInt_FromLong  PyLong_FromLong
#else
#define Py_hash_t       long    /* type was introduced with Python 3 */
#endif

#ifdef _WIN32                   /* if Microsoft Windows system */
#define THREAD          HANDLE     /* threads identified by handles */
#define THREAD_OK       0          /* return value is DWORD */
#define WORKERDEF(n,p)  DWORD WINAPI n (LPVOID p)
#else                           /* if Linux/Unix system */
#define THREAD          pthread_t  /* use the POSIX thread type */
#define THREAD_OK       NULL       /* return value is void* */
#define WORKERDEF(n,p)  void*        n (void* p)
#endif                          /* definition of a worker function */

/*----------------------------------------------------------------------
  Type Definitions
----------------------------------------------------------------------*/
typedef struct {                /* --- item set report data --- */
  PyObject  *res;               /* constructed result */
  int       err;                /* error flag */
  int       cnt;                /* number of value indicators */
  CCHAR     *rep;               /* indicators of values to report */
} REPDATA;                      /* (item set report data) */

typedef struct {                /* --- thread worker data --- */
  TABAG     *tabag;             /* transaction bag to analyze */
  TABAG     *tasur;             /* buffer for surrogate data set */
  long      cnt;                /* number of surrogate data sets */
  TBGSURRFN *surrfn;            /* surrogate data generator function */
  RNG       *rng;               /* random number generator */
  int       target;             /* target for freq. item set mining */
  SUPP      smin;               /* minimum support of an item set */
  ISREPORT  *isrep;             /* item set reporter (one per thread) */
  int       err;                /* error indicator */
  volatile long *comp;          /* number of completed data sets */
} WORKDATA;                     /* (thread worker data) */

/*----------------------------------------------------------------------
  Global Variables
----------------------------------------------------------------------*/
static TBGSURRFN *sur_tab[] = { /* surrogate generation functions */
  /* IDENT   0 */  tbg_ident,   /* identity (keep original data) */
  /* RANDOM  1 */  tbg_random,  /* transaction randomization */
  /* SWAP    2 */  tbg_swap,    /* permutation by pair swaps */
  /* SHUFFLE 3 */  tbg_shuffle, /* shuffle table-derived data */
};

static volatile int aborted = 0;/* whether abort interrupt received */

/*----------------------------------------------------------------------
  CPUinfo Functions
----------------------------------------------------------------------*/
#ifdef _WIN32                   /* if Microsoft Windows system */
#define cpuid   __cpuid         /* map existing function */
#elif !defined _SC_NPROCESSORS_ONLN    /* if Linux/Unix system */

static void cpuid (int info[4], int type)
{                               /* --- get CPU information */
  __asm__ __volatile__ (        /* execute cpuid instruction */
    "cpuid": "=a" (info[0]), "=b" (info[1]),
             "=c" (info[2]), "=d" (info[3]) : "a" (type) );
}  /* cpuid() */

#endif  /* #ifdef _WIN32 .. #else .. */
/*--------------------------------------------------------------------*/

int cpucnt (void)
{                               /* --- get the number of processors */
  #ifdef _WIN32                 /* if Microsoft Windows system */
  SYSTEM_INFO sysinfo;          /* system information structure */
  GetSystemInfo(&sysinfo);      /* get system information */
  return sysinfo.dwNumberOfProcessors;
  #elif defined _SC_NPROCESSORS_ONLN   /* if Linux/Unix system */
  return (int)sysconf(_SC_NPROCESSORS_ONLN);
  #else                         /* if no direct function available */
  if (!cpuinfo[4]) { cpuid(cpuinfo, 1); cpuinfo[4] = -1; }
  return (cpuinfo[1] >> 16) & 0xff;
  #endif                        /* use the cpuid instruction */
}  /* cpucnt() */

/*----------------------------------------------------------------------
  Signal Handler Functions
----------------------------------------------------------------------*/
#ifdef _WIN32

static BOOL WINAPI sighandler (DWORD type)
{ if (type == SIGINT) aborted = -1; return TRUE; }

static void siginstall (void)
{ SetConsoleCtrlHandler(sighandler, TRUE); }

static void sigremove (void)
{ SetConsoleCtrlHandler(sighandler, FALSE); }

#else /*--------------------------------------------------------------*/

static void sighandler (int type)
{ if (type == SIGINT) aborted = -1; }

static void siginstall (void)
{ signal(SIGINT, sighandler); }

static void sigremove (void)
{ }

#endif
/*----------------------------------------------------------------------
  Frequent Item Set Mining Functions
----------------------------------------------------------------------*/

static size_t hashitem (const void *a)
{                               /* --- compute hash code */
  return (size_t)PyObject_Hash(*(PyObject**)a);
}  /* hashitem() */

/*--------------------------------------------------------------------*/

static int cmpitems (const void *a, const void *b, void *data)
{                               /* --- compare two Python objects */
  return PyObject_RichCompareBool(*(PyObject**)a,*(PyObject**)b,Py_NE);
}  /* cmpitems() */             /* return 0 if objects are equal */

/*--------------------------------------------------------------------*/

static void delitem (const void *a)
{                              /* --- delete an object reference */
  Py_DECREF(idm_key(a));       /*     (object in the item base) */
}  /* delitem() */

/*--------------------------------------------------------------------*/

static void cleanup (TABAG *tabag, PyObject *a,
                     PyObject *b, PyObject *c, PyObject *d)
{                              /* --- clean up after an error */
  if (a) Py_DECREF(a); if (b) Py_DECREF(b);
  if (c) Py_DECREF(c); if (d) Py_DECREF(d);
  tbg_delete(tabag, 1);        /* drop references, delete trans. bag */
}  /* cleanup() */

/*--------------------------------------------------------------------*/

static TABAG* tbg_fromPyObj (PyObject *tracts)
{                               /* --- create a transaction bag */
  PyObject  *ti, *ii;           /* transaction and item iterator */
  PyObject  *trans;             /* to traverse the transactions */
  PyObject  *item;              /* to traverse the items */
  PyObject  *mul;               /* to get transaction multiplicity */
  Py_hash_t h;                  /* hash value of item */
  ITEM      n, k;               /* number of items, buffers */
  SUPP      w;                  /* weight/support buffer */
  int       isdict;             /* flag for transaction dictionary */
  TABAG     *tabag;             /* created transaction bag */
  ITEMBASE  *ibase;             /* underlying item base */

  assert(tracts);               /* check the function argument */
  ti = PyObject_GetIter(tracts);/* get an iterator for transactions */
  if (!ti) { ERR_TYPE("transaction database must be iterable"); }
  isdict = PyDict_Check(tracts);
  ibase = ib_create(IB_OBJNAMES, 0, hashitem, cmpitems, NULL, delitem);
  if (!ibase) { ERR_MEM(); }    /* create an item base */
  tabag = tbg_create(ibase);    /* and a transaction bag */
  if (!tabag) { ib_delete(ibase); ERR_MEM(); }
  while ((trans = PyIter_Next(ti))) {
    ib_clear(ibase);            /* traverse the transactions */
    ii = PyObject_GetIter(trans);
    if (!ii) { cleanup(tabag, NULL, NULL, trans, ti);
      ERR_TYPE("transactions must be iterable"); }
    if (!isdict) w = 1;         /* default: unit transaction weight */
    else {                      /* if trans. multiplicities given */
      mul = PyDict_GetItem(tracts, trans);
      if      (PyInt_Check (mul)) w = (SUPP)PyInt_AsLong (mul);
      else if (PyLong_Check(mul)) w = (SUPP)PyLong_AsLong(mul);
      else { cleanup(tabag, mul, ii, trans, ti);
        ERR_TYPE("transaction multiplicities must be integer"); }
      Py_DECREF(mul);           /* get transaction multiplicity */
    }                           /* and drop multiplicity reference */
    Py_DECREF(trans);           /* drop the transaction reference */
    while ((item = PyIter_Next(ii))) {
      h = PyObject_Hash(item);  /* check whether item is hashable */
      if (h == -1) { cleanup(tabag, NULL, item, ii, ti);
        ERR_TYPE("items must be hashable"); }
      n = ib_cnt(ibase);        /* add item to transaction */
      k = ib_add2ta(ibase, &item);
      if (ib_cnt(ibase) <= n) Py_DECREF(item);
      if (k < 0) { cleanup(tabag, NULL, NULL, ii, ti); ERR_MEM(); }
    }                           /* check for an error */
    Py_DECREF(ii);              /* drop the item iterator and */
    ib_finta(ibase, w);         /* set the transaction weight */
    if (PyErr_Occurred()) {     /* check for an iteration error */
      cleanup(tabag, NULL, NULL, NULL, ti); return NULL; }
    if (tbg_addib(tabag) < 0) { /* add the transaction to the bag */
      cleanup(tabag, NULL, NULL, NULL, ti); ERR_MEM(); }
  }
  Py_DECREF(ti);                /* drop the transaction iterator */
  if (PyErr_Occurred()) { tbg_delete(tabag, 1); return NULL; }
  return tabag;                 /* return the created trans. bag */
}  /* tbg_fromPyObj() */

/*--------------------------------------------------------------------*/

static void* isr_pyborder (ISREPORT *rep, PyObject *border)
{                               /* --- set reporter filtering border */
  int        e = 0;             /* error flag for number conversion */
  Py_ssize_t n;                 /* loop variable, sequence length */
  PyObject   *o;                /* to traverse the sequence elements */
  double     s;                 /* support threshold as a double */
  RSUPP      supp;              /* support threshold (scaled) */

  assert(rep && border);        /* check the function arguments */
  if (!PySequence_Check(border)) {
    ERR_TYPE("border must be a list or tuple of numbers"); }
  n = PySequence_Length(border);/* check for a sequence */
  if (n <= 0) return 0;         /* empty sequences need no processing */
  while (n > 0) { --n;          /* traverse the sequence elements */
    o = PySequence_GetItem(border, n);
    if      (PyLong_Check(o))   /* if element is a long integer */
      supp = (RSUPP)PyLong_AsLong(o);
    else if (PyInt_Check(o))    /* if element is an integer */
      supp = (RSUPP)PyInt_AsLong(o);
    else if (PyFloat_Check(o)){ /* if element is a float */
      s = PyFloat_AsDouble(o);
      supp = (s >= (double)SUPP_MAX) ? RSUPP_MAX : (RSUPP)s; }
    else e = 1;                 /* get the element value as support */
    Py_DECREF(o);               /* drop the element reference */
    if (e) { ERR_TYPE("border must be a list or tuple of numbers"); }
    if (isr_setbdr(rep, (ITEM)n, supp) < 0) { ERR_MEM(); }
  }                             /* set the border support */
  return (void*)1;              /* return 'ok' */
}  /* isr_pyborder() */

/*--------------------------------------------------------------------*/

static void isr_iset2PyObj (ISREPORT *rep, void *data)
{                               /* --- report an item set */
  int      i;                   /* loop variable */
  ITEM     k, n;                /* loop variables, number of items */
  RSUPP    supp, base;          /* item set support and base support */
  SUPP     s;                   /* support value for reporting */
  double   e, x;                /* evaluation and scaling factor */
  PyObject *pair;               /* pair of item set and values */
  PyObject *iset;               /* found item set (as a tuple) */
  PyObject *obj;                /* current item, to create objects */
  PyObject *vals;               /* values associated to item set */
  REPDATA  *rd = data;          /* report data structure */

  assert(rep && data);          /* check the function arguments */
  n    = isr_cnt(rep);          /* get the size of the item set */
  iset = PyTuple_New(n);        /* create an item set tuple */
  if (!iset) { rd->err = -1; return; }
  for (k = 0; k < n; k++) {     /* traverse the items */
    obj = (PyObject*)isr_itemobj(rep, isr_itemx(rep, k));
    Py_INCREF(obj);             /* get the corresp. Python object */
    PyTuple_SET_ITEM(iset, k, obj);
  }                             /* store the item in the set */
  vals = PyTuple_New(rd->cnt);  /* create a value tuple */
  if (!vals) { Py_DECREF(iset); rd->err = -1; return; }
  supp = isr_supp(rep);         /* get the item set support */
  base = isr_suppx(rep, 0);     /* and the total transaction weight */
  s = 0; e = 0;                 /* initialize the report variables */
  for (i = 0; i < rd->cnt; i++){/* traverse the values to store */
    switch (rd->rep[i]) {       /* evaluate the value indicator */
      case 'a': s = (SUPP)supp;                 x =   0.0; break;
      case 's': e = (double)supp /(double)base; x =   1.0; break;
      case 'S': e = (double)supp /(double)base; x = 100.0; break;
      case 'p': e = isr_eval(rep);              x =   1.0; break;
      case 'P': e = isr_eval(rep);              x = 100.0; break;
      case 'e': e = isr_eval(rep);              x =   1.0; break;
      case 'E': e = isr_eval(rep);              x = 100.0; break;
      default : s = 0;                          x =   0.0; break;
    }                           /* get the requested value */
    if (x == 0) obj = PyInt_FromLong((long)s);
    else        obj = PyFloat_FromDouble(e *x);
    if (!obj) { Py_DECREF(iset); Py_DECREF(vals); rd->err = -1; return;}
    PyTuple_SET_ITEM(vals, i, obj);
  }                             /* store the created value */
  pair = PyTuple_New(2);        /* create a pair of set and values */
  if (!pair) { Py_DECREF(iset); Py_DECREF(vals); rd->err = -1; return; }
  PyTuple_SET_ITEM(pair, 0, iset);
  PyTuple_SET_ITEM(pair, 1, vals);
  if (PyList_Append(rd->res, pair) != 0)
    rd->err = -1;               /* append the pair to the result list */
  Py_DECREF(pair);              /* remove internal reference to pair */
}  /* isr_iset2PyObj() */

/*--------------------------------------------------------------------*/

static double lift (RSUPP supp, RSUPP body, RSUPP head, RSUPP base)
{                               /* --- compute lift value of a rule */
  return ((body <= 0) || (head <= 0)) ? 0
       : ((double)supp*(double)base) /((double)body*(double)head);
}  /* lift() */

/*--------------------------------------------------------------------*/

static void isr_rule2PyObj (ISREPORT *rep, void *data,
                            ITEM item, RSUPP body, RSUPP head)
{                               /* --- report an association rule */
  int      v;                   /* loop variable for values */
  ITEM     i, k, n;             /* loop variables, number of items */
  ITEM     z;                   /* to traverse the items */
  RSUPP    supp, base;          /* item set support and base support */
  SUPP     s;                   /* support value for reporting */
  double   e, x;                /* evaluation and scaling factor */
  PyObject *rule;               /* triplet of head, body and values */
  PyObject *ante;               /* antecedent of association rule */
  PyObject *cons;               /* consequent of association rule */
  PyObject *obj;                /* current item, to create objects */
  PyObject *vals;               /* values associated to rule */
  REPDATA  *rd = data;          /* report data structure */

  assert(rep && data            /* check the function arguments */
  &&    (body > 0) && (head > 0));
  assert(isr_uses(rep, item));  /* head item must be in item set */
  n    = isr_cnt(rep);          /* get the size of the item set */
  ante = PyTuple_New(n-1);      /* create a tuple for the rule body */
  if (!ante) { rd->err = -1; return; }
  for (i = k = 0; i < n; i++) { /* traverse the items */
    z = isr_itemx(rep, i);      /* get the next item and skip it */
    if (z == item) continue;    /* if it is the head of the rule */
    obj = (PyObject*)isr_itemobj(rep, z);
    Py_INCREF(obj);             /* get the corresp. Python object */
    PyTuple_SET_ITEM(ante, k, obj);
    k++;                        /* store the item in the rule body */
  }                             /* and advance the item position */
  vals = PyTuple_New(rd->cnt);  /* create a value tuple */
  if (!vals) { Py_DECREF(ante); rd->err = -1; return; }
  supp = isr_supp(rep);         /* get the item set support */
  base = isr_suppx(rep, 0);     /* and the total transaction weight */
  s = 0; e = 0;                 /* initialize the report variables */
  for (v = 0; v < rd->cnt; v++){/* traverse the values to store */
    switch (rd->rep[v]) {       /* evaluate the value indicator */
      case 'a': s = (SUPP)supp;                   x =   0.0; break;
      case 'b': s = (SUPP)body;                   x =   0.0; break;
      case 'h': s = (SUPP)head;                   x =   0.0; break;
      case 's': e = (double)supp /(double)base;   x =   1.0; break;
      case 'S': e = (double)supp /(double)base;   x = 100.0; break;
      case 'x': e = (double)body /(double)base;   x =   1.0; break;
      case 'X': e = (double)body /(double)base;   x = 100.0; break;
      case 'y': e = (double)head /(double)base;   x =   1.0; break;
      case 'Y': e = (double)head /(double)base;   x = 100.0; break;
      case 'c': e = (double)supp /(double)body;   x =   1.0; break;
      case 'C': e = (double)supp /(double)body;   x = 100.0; break;
      case 'l': e = lift(supp, body, head, base); x =   1.0; break;
      case 'L': e = lift(supp, body, head, base); x = 100.0; break;
      case 'e': e = isr_eval(rep);                x =   1.0; break;
      case 'E': e = isr_eval(rep);                x = 100.0; break;
      default : s = 0;                            x =   0.0; break;
    }                           /* get the requested value */
    if (x == 0) obj = PyInt_FromLong((long)s);
    else        obj = PyFloat_FromDouble(e *x);
    if (!obj) { Py_DECREF(ante); Py_DECREF(vals); rd->err = -1; return;}
    PyTuple_SET_ITEM(vals, v, obj);
  }                             /* store the created value */
  cons = (PyObject*)isr_itemobj(rep, item);
  Py_INCREF(cons);              /* get rule head as a Python object */
  rule = PyTuple_New(3);        /* create a pair of set and values */
  if (!rule) { Py_DECREF(cons); Py_DECREF(ante);
               Py_DECREF(vals); rd->err = -1; return; }
  PyTuple_SET_ITEM(rule, 0, cons);
  PyTuple_SET_ITEM(rule, 1, ante);
  PyTuple_SET_ITEM(rule, 2, vals);
  if (PyList_Append(rd->res, rule) != 0)
    rd->err = -1;               /* append the pair to the result list */
  Py_DECREF(rule);              /* remove internal reference to rule */
}  /* isr_rule2PyObj() */

/*--------------------------------------------------------------------*/

static PyObject* psp_toPyObj (PATSPEC *psp, double scale, int format)
{                               /* --- report pattern spectrum */
  int      e = 0;               /* error indicator */
  size_t   i = 0;               /* result list index */
  ITEM     size;                /* loop variable for sizes */
  SUPP     supp, min, max;      /* loop variable for supports */
  size_t   frq;                 /* frequency of a pattern signature */
  PyObject *res;                /* created Python list object */
  PyObject *z, *s, *f, *t;      /* (elements of) result triplet */

  assert(psp);                  /* check the function arguments */
  if (format == '=') res = PyList_New((Py_ssize_t)psp_sigcnt(psp));
  else               res = PyDict_New();
  if (!res) return NULL;        /* create the pattern spectrum object */
  for (size = psp_min(psp); size <= psp_max(psp); size++) {
    min = psp_min4sz(psp,size); /* traverse the pattern sizes */
    max = psp_max4sz(psp,size); /* get range of support values */
    if (max < min) continue;    /* skip size with empty support range */
    z = PyInt_FromLong((long)size);
    if (!z) { e = -1; break; }  /* create size object for storing */
    for (supp = min; supp <= max; supp++) {
      if ((frq = psp_getfrq(psp, size, supp)) <= 0)
        continue;               /* traverse the support values */
      s = PyInt_FromLong((long)supp);
      if (!s) {                             e = -1; break; }
      f = PyFloat_FromDouble((double)frq *scale);
      if (!f) {               Py_DECREF(s); e = -1; break; }
      if (format == '=') {      /* if pattern spectrum is list */
        t = PyTuple_New(3);     /* create a result list element */
        if (!t) { Py_DECREF(f); Py_DECREF(s); e = -1; break; }
        PyTuple_SET_ITEM(t,0,z); Py_INCREF(z);
        PyTuple_SET_ITEM(t,1,s);/* fill element with the triplet */
        PyTuple_SET_ITEM(t,2,f);/* (size, support, frequency) */
        PyList_SET_ITEM(res, i, t);
        i++; }                  /* set element in the result list */
      else {                    /* if pattern spectrum is dictionary */
        t = PyTuple_New(2);     /* create a tuple as a key */
        if (!t) { Py_DECREF(f); Py_DECREF(s); e = -1; break; }
        PyTuple_SET_ITEM(t,0,z); Py_INCREF(z);
        PyTuple_SET_ITEM(t,1,s);/* create key as pair (size, support) */
        PyDict_SetItem(res, t, f);
      }                         /* map key to the counter */
    }
    Py_DECREF(z);               /* initial reference no longer needed */
    if (e) break;               /* if an error occurred, abort loop */
  }
  if (e) { Py_DECREF(res); return NULL; }
  return res;                   /* return created pattern spectrum */
}  /* psp_toPyObj() */

/*--------------------------------------------------------------------*/

static int get_target (const char *s, const char *targets)
{                               /* --- get target */
  if      (strcmp(s, "sets")       == 0) s = "s";
  else if (strcmp(s, "all")        == 0) s = "s";
  else if (strcmp(s, "frequent")   == 0) s = "s";
  else if (strcmp(s, "cls")        == 0) s = "c";
  else if (strcmp(s, "clsd")       == 0) s = "c";
  else if (strcmp(s, "closed")     == 0) s = "c";
  else if (strcmp(s, "max")        == 0) s = "m";
  else if (strcmp(s, "maxi")       == 0) s = "m";
  else if (strcmp(s, "maximal")    == 0) s = "m";
  else if (strcmp(s, "gen")        == 0) s = "g";
  else if (strcmp(s, "gens")       == 0) s = "g";
  else if (strcmp(s, "generas")    == 0) s = "g";
  else if (strcmp(s, "generators") == 0) s = "g";
  else if (strcmp(s, "rule")       == 0) s = "r";
  else if (strcmp(s, "rules")      == 0) s = "r";
  else if (strcmp(s, "arule")      == 0) s = "r";
  else if (strcmp(s, "arules")     == 0) s = "r";
  if ((strlen(s) == 1)          /* translate the target string */
  &&  (strchr(targets, s[0]) != NULL)) {
    switch (s[0]) {             /* evaluate the target code */
      case 'a': return ISR_SETS;
      case 's': return ISR_SETS;
      case 'c': return ISR_CLOSED;
      case 'm': return ISR_MAXIMAL;
      case 'g': return ISR_GENERAS;
      case 'r': return ISR_RULES;
    }
  }
  PyErr_SetString(PyExc_ValueError, "invalid target type");
  return -1;                    /* return an error code */
}  /* get_target() */

/*--------------------------------------------------------------------*/

static int get_stat (const char *s)
{                               /* --- get statistic code */
  if      (strcmp(s, "none")      == 0) s = "x";
  else if (strcmp(s, "chi2")      == 0) s = "p";
  else if (strcmp(s, "chi2pval")  == 0) s = "p";
  else if (strcmp(s, "yates")     == 0) s = "t";
  else if (strcmp(s, "yatespval") == 0) s = "t";
  else if (strcmp(s, "info")      == 0) s = "g";
  else if (strcmp(s, "infopval")  == 0) s = "g";
  else if (strcmp(s, "fetprob")   == 0) s = "f";
  else if (strcmp(s, "fetchi2")   == 0) s = "h";
  else if (strcmp(s, "fetinfo")   == 0) s = "m";
  else if (strcmp(s, "fetsupp")   == 0) s = "s";
  if (strlen(s) == 1) {         /* translate the statistic string */
    switch (s[0]) {             /* evaluate the statistic code */
      case 'x': return RE_NONE;
      case 'c': return RE_CHI2PVAL;
      case 'p': return RE_CHI2PVAL;
      case 'n': return RE_CHI2PVAL;
      case 'y': return RE_YATESPVAL;
      case 't': return RE_YATESPVAL;
      case 'i': return RE_INFOPVAL;
      case 'g': return RE_INFOPVAL;
      case 'f': return RE_FETPROB;
      case 'h': return RE_FETCHI2;
      case 'm': return RE_FETINFO;
      case 's': return RE_FETSUPP;
    }
  }
  PyErr_SetString(PyExc_ValueError, "invalid statistic");
  return -1;                    /* return an error code */
}  /* get_stat() */

/*--------------------------------------------------------------------*/

static int get_eval (const char *s)
{                               /* --- get evaluation measure code */
  if (strcmp(s, "none")    == 0) return 'x';
  if (strcmp(s, "x")       == 0) return 'x';
  if (strcmp(s, "ldratio") == 0) return 'b';
  if (strcmp(s, "b")       == 0) return 'b';
  PyErr_SetString(PyExc_ValueError, "invalid evaluation measure");
  return -1;                    /* return an error code */
}  /* get_eval() */

/*--------------------------------------------------------------------*/

static int get_evalx (const char *s)
{                               /* --- get evaluation measure code */
  if (strcmp(s, "none")       == 0) s = "x";
  if (strcmp(s, "conf")       == 0) s = "c";
  if (strcmp(s, "confidence") == 0) s = "c";
  if (strcmp(s, "confdiff")   == 0) s = "d";
  if (strcmp(s, "lift")       == 0) s = "l";
  if (strcmp(s, "liftdiff")   == 0) s = "a";
  if (strcmp(s, "liftquot")   == 0) s = "q";
  if (strcmp(s, "cvct")       == 0) s = "v";
  if (strcmp(s, "conviction") == 0) s = "v";
  if (strcmp(s, "cvctdiff")   == 0) s = "e";
  if (strcmp(s, "cvctquot")   == 0) s = "r";
  if (strcmp(s, "cprob")      == 0) s = "k";
  if (strcmp(s, "import")     == 0) s = "j";
  if (strcmp(s, "importance") == 0) s = "j";
  if (strcmp(s, "cert")       == 0) s = "z";
  if (strcmp(s, "chi2")       == 0) s = "n";
  if (strcmp(s, "chi2pval")   == 0) s = "p";
  if (strcmp(s, "yates")      == 0) s = "y";
  if (strcmp(s, "yatespval")  == 0) s = "t";
  if (strcmp(s, "info")       == 0) s = "i";
  if (strcmp(s, "infopval")   == 0) s = "g";
  if (strcmp(s, "fetprob")    == 0) s = "f";
  if (strcmp(s, "fetchi2")    == 0) s = "h";
  if (strcmp(s, "fetinfo")    == 0) s = "m";
  if (strcmp(s, "fetsupp")    == 0) s = "s";
  if (strcmp(s, "ldratio")    == 0) s = "b";
  if (strlen(s) == 1) {         /* translate the measure string */
    switch (s[0]) {             /* evaluate the measure code */
      case 'x': return RE_NONE;
      case 'c': return RE_CONF;
      case 'd': return RE_CONFDIFF;
      case 'l': return RE_LIFT;
      case 'a': return RE_LIFTDIFF;
      case 'q': return RE_LIFTQUOT;
      case 'v': return RE_CVCT;
      case 'e': return RE_CVCTDIFF;
      case 'r': return RE_CVCTQUOT;
      case 'k': return RE_CPROB;
      case 'j': return RE_IMPORT;
      case 'z': return RE_CERT;
      case 'n': return RE_CHI2;
      case 'p': return RE_CHI2PVAL;
      case 'y': return RE_YATES;
      case 't': return RE_YATESPVAL;
      case 'i': return RE_INFO;
      case 'g': return RE_INFOPVAL;
      case 'f': return RE_FETPROB;
      case 'h': return RE_FETCHI2;
      case 'm': return RE_FETINFO;
      case 's': return RE_FETSUPP;
      case 'b': return RE_FNCNT;
    }
  }
  PyErr_SetString(PyExc_ValueError, "invalid evaluation measure");
  return -1;                    /* return an error code */
}  /* get_evalx() */

/*--------------------------------------------------------------------*/

static int get_agg (const char *s)
{                               /* --- get aggregation mode */
  if      (strcmp(s, "none") == 0) s = "x";
  else if (strcmp(s, "min")  == 0) s = "m";
  else if (strcmp(s, "max")  == 0) s = "n";
  else if (strcmp(s, "avg")  == 0) s = "a";
  if (strlen(s) == 1) {         /* translate the aggregation string */
    switch (s[0]) {             /* evaluate the aggregation code */
      case 'x': return IST_NONE;
      case 'm': return IST_MIN;
      case 'n': return IST_MAX;
      case 'a': return IST_AVG;
    }
  }
  PyErr_SetString(PyExc_ValueError, "invalid aggregation mode");
  return -1;                    /* return an error code */
}  /* get_agg() */

/*--------------------------------------------------------------------*/

static int get_surr (const char *s)
{                               /* --- get surrogate function code */
  if      (strcmp(s, "ident")     == 0) s = "i";
  else if (strcmp(s, "identity")  == 0) s = "i";
  else if (strcmp(s, "random")    == 0) s = "r";
  else if (strcmp(s, "randomize") == 0) s = "r";
  else if (strcmp(s, "swap")      == 0) s = "p";
  else if (strcmp(s, "perm")      == 0) s = "p";
  else if (strcmp(s, "permute")   == 0) s = "p";
  else if (strcmp(s, "shuffle")   == 0) s = "s";
  if (strlen(s) == 1) {         /* translate surrogate method string */
    switch (s[0]) {             /* evaluate the surrogate method code */
      case 'i': return 0;
      case 'r': return 1;
      case 'p': return 2;
      case 'w': return 2;
      case 's': return 3;
    }
  }
  PyErr_SetString(PyExc_ValueError,
                  "invalid surrogate generation method");
  return -1;                    /* return an error code */
}  /* get_surr() */

/*--------------------------------------------------------------------*/

static int repinit (REPDATA *data, ISREPORT *isrep, CCHAR *report,
                    int target)
{                               /* --- initialize reporting */
  assert(data && isrep && report); /* check the function arguments */
  data->err = 0;                /* initialize the error indicator */
  if ((report[0] == '#')        /* if to get a pattern spectrum */
  ||  (report[0] == '='))       /* #: dictionary, =: list of triplets */
    return isr_addpsp(isrep, NULL);
  data->cnt = (int)strlen(data->rep = report);
  data->res = PyList_New(0);    /* create an empty output list */
  if (!data->res) return -1;    /* and set the reporting function */
  if (target & ISR_RULES) isr_setrule(isrep, isr_rule2PyObj, data);
  else                    isr_setrepo(isrep, isr_iset2PyObj, data);
  return 0;                     /* return 'ok' */
}  /* repinit() */

/*--------------------------------------------------------------------*/

static int repterm (REPDATA *data, ISREPORT *isrep, CCHAR *report)
{                               /* --- terminate reporting */
  assert(data && isrep && report); /* check the function arguments */
  if ((report[0] == '#')        /* if to get a pattern spectrum */
  ||  (report[0] == '=')) {     /* #: dictionary, =: list of triplets */
    data->res = psp_toPyObj(isr_getpsp(isrep), 1.0, report[0]);
    return data->err = (data->res) ? 0 : -1;
  }                             /* make Python pattern spectrum */
  return data->err;             /* return the error status */
}  /* repterm() */

/*--------------------------------------------------------------------*/
/* fim (tracts, target='s', supp=10, zmin=1, zmax=None,               */
/*      report='a', eval='x', agg='x', thresh=10, border=None)        */
/*--------------------------------------------------------------------*/

static PyObject* py_fim (PyObject *self,
                         PyObject *args, PyObject *kwds)
{                               /* --- frequent item set mining */
  char     *ckwds[] = { "tracts", "target", "supp", "zmin", "zmax",
                        "report", "eval", "agg", "thresh", "border",
                         NULL };
  CCHAR    *starg  = "s";       /* target type as a string */
  int      target  = ISR_SETS;  /* target type */
  double   supp    = 10;        /* minimum support of an item set */
  SUPP     smin    =  1;        /* minimum support as an integer */
  long     zmin    =  1;        /* minimum size of an item set */
  long     zmax    = -1;        /* maximum size of an item set */
  CCHAR    *report = "a";       /* indicators of values to report */
  CCHAR    *seval  = "x";       /* evaluation measure as a string */
  int      eval    = 'x';       /* evaluation measure */
  CCHAR    *sagg   = "x";       /* aggregation mode as a string */
  int      agg     =  0;        /* aggregation mode */
  double   thresh  = 10;        /* threshold for evaluation measure */
  int      algo    = FPG_SIMPLE;            /* algorithm variant */
  int      mode    = FPG_DEFAULT|FPG_FIM16; /* operation mode/flags */
  long     prune   = LONG_MIN;  /* min. size for evaluation filtering */
  PyObject *border = NULL;      /* support border for filtering */
  PyObject *tracts;             /* transaction database */
  TABAG    *tabag;              /* transaction bag */
  ISREPORT *isrep;              /* item set reporter */
  REPDATA  data;                /* data for item set reporting */
  double   w;                   /* total transaction weight */
  int      r;                   /* result buffer */

  /* --- evaluate the function arguments --- */
  if (!PyArg_ParseTupleAndKeywords(args, kwds, "O|sdllsssdO", ckwds,
        &tracts, &starg, &supp, &zmin, &zmax, &report,
        &seval, &sagg, &thresh, &border))
    return NULL;                /* parse the function arguments */
  target = get_target(starg, "ascmg");
  if (target < 0) return NULL;  /* translate the target string */
  if (zmin   < 0)    { ERR_VALUE("zmin must not be negative"); }
  if (zmax   < 0)    zmax = LONG_MAX;  /* check size range */
  if (zmax   < zmin) { ERR_VALUE("zmax must not be less than zmin"); }
  eval = get_evalx(seval);      /* get evaluation measure and */
  if (eval   < 0) return NULL;  /* check whether it is valid */
  if (eval   <= RE_NONE) prune = LONG_MIN;
  agg  = get_agg(sagg);         /* get aggregation mode and */
  if (agg    < 0) return NULL;  /* check whether it is valid */
  thresh *= 0.01;               /* scale evaluation threshold */

  /* --- create transaction bag --- */
  tabag = tbg_fromPyObj(tracts);/* turn the given transactions */
  if (!tabag) return NULL;      /* into a transaction bag */
  w    = tbg_wgt(tabag);        /* get the total transaction weight */
  supp = (supp >= 0) ? 0.01 *supp *(double)w *(1-DBL_EPSILON) : -supp;
  smin = (SUPP)ceilsupp(supp);  /* compute absolute support value */
  r = fpg_data(tabag, target, smin, (ITEM)zmin, eval, algo, mode, 2);
  if (r) tbg_delete(tabag, 1);  /* prepare data for FP-growth */
  if (r == -1) { ERR_MEM(); }   /* check for error and no items */
  if (r <   0) return PyList_New(0);

  /* --- create item set reporter --- */
  isrep = isr_create(tbg_base(tabag)); /* create item set reporter */
  if (!isrep) { tbg_delete(tabag, 1); ERR_MEM(); }
  isr_setsize(isrep, (zmin > ITEM_MAX) ? ITEM_MAX : (ITEM)zmin,
                     (zmax > ITEM_MAX) ? ITEM_MAX : (ITEM)zmax);
  isr_setsupp(isrep, (RSUPP)smin, RSUPP_MAX);
  if (border && !isr_pyborder(isrep, border)) {
    isr_delete(isrep, 0); tbg_delete(tabag, 1); return NULL; }
  if ((repinit(&data, isrep, report, ISR_SETS) != 0)
  ||  (fpg_repo(isrep, target, eval, thresh, algo, mode) < 0)
  ||  (isr_setup(isrep) < 0)) { /* set up the item set reporter */
    isr_delete(isrep, 0); tbg_delete(tabag, 1); ERR_MEM(); }

  /* --- frequent item set mining --- */
  r = fpgrowth(tabag, target, smin, smin, 1, eval, agg, thresh,
              (prune < ITEM_MIN) ? ITEM_MIN :
              (prune > ITEM_MAX) ? ITEM_MAX : (ITEM)prune,
               algo, mode, 0, isrep);
  if (r >= 0) r = repterm(&data, isrep, report);

  /* --- clean up --- */
  isr_delete(isrep, 0);         /* delete the item set reporter */
  tbg_delete(tabag, 1);         /* and the transaction bag */
  if (r < 0) { Py_DECREF(data.res); ERR_MEM(); }
  return data.res;              /* return the created result */
}  /* py_fim() */

/*--------------------------------------------------------------------*/
/* arules (tracts, supp=10, conf=80, zmin=1, zmax=None, report='aC',  */
/*         eval='x', thresh=10, mode='')                              */
/*--------------------------------------------------------------------*/

static PyObject* py_arules (PyObject *self,
                            PyObject *args, PyObject *kwds)
{                               /* --- association rule mining */
  char     *ckwds[] = { "tracts", "supp", "conf", "zmin", "zmax",
                        "report", "eval", "thresh", "mode", NULL };
  double   supp    = 10;        /* minimum support    of a rule */
  SUPP     smin    =  1;        /* minimum support of an item set */
  SUPP     body    =  1;        /* minimum support of a rule body */
  double   conf    = 80;        /* minimum confidence of a rule */
  long     zmin    =  1;        /* minimum size of an item set */
  long     zmax    = -1;        /* maximum size of an item set */
  CCHAR    *report = "aC";      /* indicators of values to report */
  CCHAR    *seval  = "x";       /* evaluation measure as a string */
  int      eval    =  0;        /* evaluation measure */
  double   thresh  = 10;        /* threshold for evaluation measure */
  int      algo    = FPG_SINGLE;/* algorithm variant */
  CCHAR    *smode  = "";        /* operation mode/flags as a string */
  int      mode    = FPG_DEFAULT|FPG_FIM16; /* operation mode/flags */
  PyObject *tracts;             /* transaction database */
  TABAG    *tabag;              /* transaction bag */
  ISREPORT *isrep;              /* item set reporter */
  REPDATA  data;                /* data for item set reporting */
  double   w;                   /* total transaction weight */
  int      r;                   /* result buffer */

  /* --- evaluate the function arguments --- */
  if (!PyArg_ParseTupleAndKeywords(args, kwds, "O|ddllssds", ckwds,
        &tracts, &supp, &conf, &zmin, &zmax, &report,
        &seval, &thresh, &smode))
    return NULL;                /* parse the function arguments */
  if ((conf < 0) || (conf > 100)) { ERR_VALUE("invalid confidence"); }
  if (zmin < 0)    { ERR_VALUE("zmin must not be negative"); }
  if (zmax < 0)    zmax = LONG_MAX; /* check size range */
  if (zmax < zmin) { ERR_VALUE("zmax must not be less than zmin"); }
  eval = get_evalx(seval);      /* get evaluation measure and */
  if (eval < 0) return NULL;    /* check whether it is valid */
  thresh *= 0.01;               /* scale evaluation threshold */

  /* --- create transaction bag --- */
  tabag = tbg_fromPyObj(tracts);/* turn the given transactions */
  if (!tabag) return NULL;      /* into a transaction bag */
  conf *= 0.01;                 /* scale the minimum confidence and */
  w    = tbg_wgt(tabag);        /* get the total transaction weight */
  supp = (supp >= 0) ? 0.01 *supp *(double)w *(1-DBL_EPSILON) : -supp;
  body = (SUPP)ceilsupp(supp);  /* compute absolute support values */
  smin = (SUPP)ceilsupp(strchr(smode, 'o') ? supp
                      : ceilsupp(supp) *conf *(1-DBL_EPSILON));
  r = fpg_data(tabag, ISR_RULES, smin, (ITEM)zmin, eval, algo, mode, 2);
  if (r) tbg_delete(tabag, 1);  /* prepare data for FP-growth */
  if (r == -1) { ERR_MEM(); }   /* check for error and no items */
  if (r <   0) return PyList_New(0);

  /* --- create item set reporter --- */
  isrep = isr_create(tbg_base(tabag)); /* create item set reporter */
  if (!isrep) { tbg_delete(tabag, 1); ERR_MEM(); }
  isr_setsize(isrep, (zmin > ITEM_MAX) ? ITEM_MAX : (ITEM)zmin,
                     (zmax > ITEM_MAX) ? ITEM_MAX : (ITEM)zmax);
  isr_setsupp(isrep, (RSUPP)smin, RSUPP_MAX);
  if ((repinit(&data, isrep, report, ISR_RULES) != 0)
  ||  (fpg_repo(isrep, ISR_RULES, eval, thresh, algo, mode) < 0)
  ||  (isr_setup(isrep) < 0)) { /* set up the item set reporter */
    isr_delete(isrep, 0); tbg_delete(tabag, 1); ERR_MEM(); }

  /* --- association rule mining --- */
  r = fpgrowth(tabag, ISR_RULES, smin, body, conf,
               eval, FPG_NONE, thresh, 0, algo, mode, 0, isrep);
  if (r >= 0) r = repterm(&data, isrep, report);

  /* --- clean up --- */
  isr_delete(isrep, 0);         /* delete the item set reporter */
  tbg_delete(tabag, 1);         /* and the transaction bag */
  if (r < 0) { Py_DECREF(data.res); ERR_MEM(); }
  return data.res;              /* return the created result */
}  /* py_arules() */

/*--------------------------------------------------------------------*/
/* apriori (tracts, target='s', supp=10, conf=80, zmin=1, zmax=None,  */
/*          report='a', eval='x', agg='x', thresh=10, prune=None,     */
/*          algo='', mode='', border=None)                            */
/*--------------------------------------------------------------------*/

static PyObject* py_apriori (PyObject *self,
                             PyObject *args, PyObject *kwds)
{                               /* --- Apriori algorithm */
  char     *ckwds[] = { "tracts", "target", "supp", "conf",
                        "zmin", "zmax", "report",
                        "eval", "agg", "thresh", "prune",
                        "algo", "mode", "border", "conf", NULL };
  CCHAR    *starg  = "s";       /* target type as a string */
  int      target  = ISR_SETS;  /* target type */
  double   supp    = 10;        /* minimum support of an item set */
  SUPP     smin    =  1;        /* minimum support as an integer */
  SUPP     body    =  1;        /* minimum support of a rule body */
  double   conf    = 80;        /* minimum confidence of a rule */
  long     zmin    =  1;        /* minimum size of an item set */
  long     zmax    = -1;        /* maximum size of an item set */
  CCHAR    *report = "a";       /* indicators of values to report */
  CCHAR    *seval  = "x";       /* evaluation measure as a string */
  int      eval    =  0;        /* evaluation measure */
  CCHAR    *sagg   = "x";       /* aggregation mode as a string */
  int      agg     =  0;        /* aggregation mode */
  double   thresh  = 10;        /* threshold for evaluation measure */
  CCHAR    *salgo  = "a";       /* algorithm as a string */
  int      algo    = APR_BASIC; /* algorithm */
  CCHAR    *smode  = "", *s;    /* operation mode/flags as a string */
  int      mode    = APR_DEFAULT;  /* operation mode/flags */
  long     prune   = LONG_MIN;  /* min. size for evaluation filtering */
  PyObject *border = NULL;      /* support border for filtering */
  PyObject *tracts;             /* transaction database */
  TABAG    *tabag;              /* transaction bag */
  ISREPORT *isrep;              /* item set reporter */
  REPDATA  data;                /* data for item set reporting */
  double   w;                   /* total transaction weight */
  int      r;                   /* result of function call */

  /* --- evaluate the function arguments --- */
  if (!PyArg_ParseTupleAndKeywords(args, kwds, "O|sddllsssdlssO", ckwds,
        &tracts, &starg, &supp, &conf, &zmin, &zmax, &report,
        &seval, &sagg, &thresh, &prune, &salgo, &smode, &border))
    return NULL;                /* parse the function arguments */
  target = get_target(starg, "ascmgr");
  if (target < 0) return NULL;  /* translate the target string */
  if (zmin   < 0)    { ERR_VALUE("zmin must not be negative"); }
  if (zmax   < 0)    zmax = LONG_MAX; /* check size range */
  if (zmax   < zmin) { ERR_VALUE("zmax must not be less than zmin"); }
  eval = get_evalx(seval);      /* get evaluation measure and */
  if (eval   < 0) return NULL;  /* check whether it is valid */
  if (eval   <= RE_NONE) prune = LONG_MIN;
  if (strchr(smode, 'z')) eval |= IST_INVBXS;
  agg  = get_agg(sagg);         /* get aggregation mode and */
  if (agg    < 0) return NULL;  /* check whether it is valid */
  thresh *= 0.01;               /* scale evaluation threshold */
  if      (strcmp(salgo, "auto")   == 0) salgo = "a";
  else if (strcmp(salgo, "basic")  == 0) salgo = "b";
  if (strlen(salgo) != 1)       /* translate the algorithm string */
    algo = -1;                  /* if it failed, set error code */
  else {                        /* if translation worked, */
    switch (salgo[0]) {         /* evaluate the algorithm code */
      case 'a': algo = APR_BASIC; break;
      case 'b': algo = APR_BASIC; break;
      default : algo = -1;        break;
    }                           /* set an error code for all */
  }                             /* other algorithm indicators */
  if (algo < 0) { ERR_VALUE("invalid Apriori algorithm"); }
  for (s = smode; *s; s++) {    /* traverse the mode characters */
    if      (*s == 'x') mode &= ~APR_PERFECT;
    else if (*s == 't') mode &= ~APR_TATREE;
    else if (*s == 'T') mode &= ~APR_TATREE;
    else if (*s == 'y') mode |=  APR_POST;
  }                             /* adapt the operation mode */

  /* --- create transaction bag --- */
  tabag = tbg_fromPyObj(tracts);/* turn the given transactions */
  if (!tabag) return NULL;      /* into a transaction bag */
  conf *= 0.01;                 /* scale the minimum confidence and */
  w    = tbg_wgt(tabag);        /* get the total transaction weight */
  supp = (supp >= 0) ? 0.01 *supp *(double)w *(1-DBL_EPSILON) : -supp;
  body = (SUPP)ceilsupp(supp);  /* compute absolute support values */
  smin = (SUPP)ceilsupp(((target & ISR_RULES) && strchr(smode, 'o'))
                     ? supp : ceilsupp(supp) *conf *(1-DBL_EPSILON));
  r = apriori_data(tabag, target, smin, (ITEM)zmin, eval, algo, mode,2);
  if (r) tbg_delete(tabag, 1);  /* prepare data for Apriori */
  if (r == -1) { ERR_MEM(); }   /* check for error and no items */
  if (r <   0) return PyList_New(0);

  /* --- create item set reporter --- */
  isrep = isr_create(tbg_base(tabag)); /* create item set reporter */
  if (!isrep) { tbg_delete(tabag, 1); ERR_MEM(); }
  isr_setsize(isrep, (zmin > ITEM_MAX) ? ITEM_MAX : (ITEM)zmin,
                     (zmax > ITEM_MAX) ? ITEM_MAX : (ITEM)zmax);
  isr_setsupp(isrep, (RSUPP)smin, RSUPP_MAX);
  if (border && !isr_pyborder(isrep, border)) {
    isr_delete(isrep, 0); tbg_delete(tabag, 1); return NULL; }
  if ((repinit(&data, isrep, report, target) != 0)
  ||  (apriori_repo(isrep, target, eval, thresh, algo, mode) < 0)
  ||  (isr_setup(isrep) < 0)) { /* set up the item set reporter */
    isr_delete(isrep, 0); tbg_delete(tabag, 1); ERR_MEM(); }

  /* --- frequent item set mining --- */
  r = apriori(tabag, target, smin, body, conf, eval, agg, thresh,
              (prune < ITEM_MIN) ? ITEM_MIN :
              (prune > ITEM_MAX) ? ITEM_MAX : (ITEM)prune,
              algo, mode, 0.01, 0, isrep);
  if (r >= 0) r = repterm(&data, isrep, report);

  /* --- clean up --- */
  isr_delete(isrep, 0);         /* delete the item set reporter */
  tbg_delete(tabag, 1);         /* and the transaction bag */
  if (r < 0) { Py_DECREF(data.res); ERR_MEM(); }
  return data.res;              /* return the created result */
}  /* py_apriori() */

/*--------------------------------------------------------------------*/
/* eclat (tracts, target='s', supp=10, conf=80, zmin=1, zmax=None,    */
/*        report='a', eval='x', agg='x', thresh=10, prune=None,       */
/*        algo='a', mode='', border=None)                             */
/*--------------------------------------------------------------------*/

static PyObject* py_eclat (PyObject *self,
                           PyObject *args, PyObject *kwds)
{                               /* --- Eclat algorithm */
  char     *ckwds[] = { "tracts", "target", "supp", "conf",
                        "zmin", "zmax", "report",
                        "eval", "agg", "thresh", "prune",
                        "algo", "mode", "border", NULL };
  CCHAR    *starg  = "s";       /* target type as a string */
  int      target  = ISR_SETS;  /* target type */
  double   supp    = 10;        /* minimum support of an item set */
  SUPP     smin    =  1;        /* minimum support as an integer */
  SUPP     body    =  1;        /* minimum support of a rule body */
  double   conf    = 80;        /* minimum confidence of a rule */
  long     zmin    =  1;        /* minimum size of an item set */
  long     zmax    = -1;        /* maximum size of an item set */
  CCHAR    *report = "a";       /* indicators of values to report */
  CCHAR    *seval  = "x";       /* evaluation measure as a string */
  int      eval    =  0;        /* evaluation measure */
  CCHAR    *sagg   = "x";       /* aggregation mode as a string */
  int      agg     =  0;        /* aggregation mode */
  double   thresh  = 10;        /* threshold for evaluation measure */
  CCHAR    *salgo  = "a";       /* algorithm as a string */
  int      algo    = ECL_OCCDLV;/* algorithm */
  CCHAR    *smode  = "", *s;    /* operation mode/flags as a string */
  int      mode    = ECL_DEFAULT|ECL_FIM16; /* operation mode/flags */
  long     prune   = LONG_MIN;  /* min. size for evaluation filtering */
  PyObject *border = NULL;      /* support border for filtering */
  PyObject *tracts;             /* transaction database */
  TABAG    *tabag;              /* transaction bag */
  ISREPORT *isrep;              /* item set reporter */
  REPDATA  data;                /* data for item set reporting */
  ITEM     m;                   /* number of items */
  double   w;                   /* total transaction weight */
  int      r;                   /* result buffer */

  /* --- evaluate the function arguments --- */
  if (!PyArg_ParseTupleAndKeywords(args, kwds, "O|sddllsssdlssO", ckwds,
        &tracts, &starg, &supp, &conf, &zmin, &zmax, &report,
        &seval, &sagg, &thresh, &prune, &salgo, &smode, &border))
    return NULL;                /* parse the function arguments */
  target = get_target(starg, "ascmgr");
  if (target < 0) return NULL;  /* translate the target string */
  if (zmin   < 0)    { ERR_VALUE("zmin must not be negative"); }
  if (zmax   < 0)    zmax = ITEM_MAX; /* check size range */
  if (zmax   < zmin) { ERR_VALUE("zmax must not be less than zmin"); }
  eval = get_evalx(seval);     /* get the evaluation measure */
  if (eval   < 0)    return NULL;
  if (eval   <= RE_NONE) prune = LONG_MIN;
  if (strchr(smode, 'z')) eval |= ECL_INVBXS;
  agg  = get_agg(sagg);         /* get the aggregation mode */
  if (agg    < 0)    return NULL;
  thresh *= 0.01;               /* add eval. flag, scale threshold */
  if      (strcmp(salgo, "auto")   == 0) salgo = "a";
  else if (strcmp(salgo, "basic")  == 0) salgo = "e";
  else if (strcmp(salgo, "lists")  == 0) salgo = "i";
  else if (strcmp(salgo, "tids")   == 0) salgo = "i";
  else if (strcmp(salgo, "bits")   == 0) salgo = "b";
  else if (strcmp(salgo, "table")  == 0) salgo = "t";
  else if (strcmp(salgo, "simple") == 0) salgo = "s";
  else if (strcmp(salgo, "ranges") == 0) salgo = "r";
  else if (strcmp(salgo, "occdlv") == 0) salgo = "o";
  else if (strcmp(salgo, "diff")   == 0) salgo = "d";
  if (strlen(salgo) != 1)       /* translate the algorithm string */
    algo = -1;                  /* if it failed, set error code */
  else {                        /* if translation worked, */
    switch (salgo[0]) {         /* evaluate the algorithm code */
      case 'a': algo = ECL_AUTO;   break;
      case 'e': algo = ECL_BASIC;  break;
      case 'i': algo = ECL_LISTS;  break;
      case 'b': algo = ECL_BITS;   break;
      case 't': algo = ECL_TABLE;  break;
      case 's': algo = ECL_SIMPLE; break;
      case 'r': algo = ECL_RANGES; break;
      case 'o': algo = ECL_OCCDLV; break;
      case 'd': algo = ECL_DIFFS;  break;
      default : algo = -1;         break;
    }                           /* set an error code for all */
  }                             /* other algorithm indicators */
  if (algo < 0) { ERR_VALUE("invalid Eclat algorithm"); }
  for (s = smode; *s; s++) {    /* traverse the mode characters */
    if      (*s == 'l') mode &= ~ECL_FIM16;
    else if (*s == 'x') mode &= ~ECL_PERFECT;
    else if (*s == 'i') mode &= ~ECL_REORDER;
    else if (*s == 'u') mode &= ~ECL_TAIL;
    else if (*s == 'y') mode |=  ECL_HORZ;
    else if (*s == 'Y') mode |=  ECL_VERT;
  }                             /* adapt the operation mode */

  /* --- create transaction bag --- */
  tabag = tbg_fromPyObj(tracts);/* turn the given transactions */
  if (!tabag) return NULL;      /* into a transaction bag */
  conf *= 0.01;                 /* scale the minimum confidence and */
  w    = tbg_wgt(tabag);        /* get the total transaction weight */
  supp = (supp >= 0) ? 0.01 *supp *(double)w *(1-DBL_EPSILON) : -supp;
  body = (SUPP)ceilsupp(supp);  /* compute absolute support values */
  smin = (SUPP)ceilsupp(((target & ISR_RULES) && strchr(smode, 'o'))
                     ? supp : ceilsupp(supp) *conf *(1-DBL_EPSILON));
  if (algo == ECL_AUTO) {       /* if automatic variant choice */
    m    = ib_frqcnt(tbg_base(tabag), smin);
    algo = ((target & (ISR_CLOSED|ISR_MAXIMAL))
        && ((double)tbg_extent(tabag) /((double)m*(double)w) > 0.02))
         ? ECL_LISTS : ECL_OCCDLV;
  }                             /* choose the eclat variant */
  r = eclat_data(tabag, target, smin, (ITEM)zmin, eval, algo, mode, 2);
  if (r) tbg_delete(tabag, 1);  /* prepare data for Eclat */
  if (r == -1) { ERR_MEM(); }   /* check for error and no items */
  if (r <   0) return PyList_New(0);

  /* --- create item set reporter --- */
  isrep = isr_create(tbg_base(tabag)); /* create item set reporter */
  if (!isrep) { tbg_delete(tabag, 1); ERR_MEM(); }
  isr_setsize(isrep, (zmin > ITEM_MAX) ? ITEM_MAX : (ITEM)zmin,
                     (zmax > ITEM_MAX) ? ITEM_MAX : (ITEM)zmax);
  isr_setsupp(isrep, (RSUPP)smin, RSUPP_MAX);
  if (border && !isr_pyborder(isrep, border)) {
    isr_delete(isrep, 0); tbg_delete(tabag, 1); return NULL; }
  if ((repinit(&data, isrep, report, target) != 0)
  ||  (eclat_repo(isrep, target, eval, thresh, algo, mode) < 0)
  ||  (isr_setup(isrep) < 0)) { /* set up the item set reporter */
    isr_delete(isrep, 0); tbg_delete(tabag, 1); ERR_MEM(); }

  /* --- frequent item set mining --- */
  r = eclat(tabag, target, smin, body, conf, eval, agg, thresh,
            (prune < ITEM_MIN) ? ITEM_MIN :
            (prune > ITEM_MAX) ? ITEM_MAX : (ITEM)prune,
            algo, mode, 0, isrep);
  if (r >= 0) r = repterm(&data, isrep, report);

  /* --- clean up --- */
  isr_delete(isrep, 0);         /* delete the item set reporter */
  tbg_delete(tabag, 1);         /* and the transaction bag */
  if (r < 0) { Py_DECREF(data.res); ERR_MEM(); }
  return data.res;              /* return the created result */
}  /* py_eclat() */

/*--------------------------------------------------------------------*/
/* fpgrowth (tracts, target='s', supp=10, conf=80, zmin=1, zmax=None, */
/*           report='a', eval='x', agg='x', thresh=10, prune=None,    */
/*           algo='s', mode='', border=None)                          */
/*--------------------------------------------------------------------*/

static PyObject* py_fpgrowth (PyObject *self,
                              PyObject *args, PyObject *kwds)
{                               /* --- FP-growth algorithm */
  char     *ckwds[] = { "tracts", "target", "supp", "conf",
                        "zmin", "zmax", "report",
                        "eval", "agg", "thresh", "prune",
                        "algo", "mode", "border", NULL };
  CCHAR    *starg  = "s";       /* target type as a string */
  int      target  = ISR_SETS;  /* target type */
  double   supp    = 10;        /* minimum support of an item set */
  SUPP     smin    =  1;        /* minimum support as an integer */
  SUPP     body    =  1;        /* minimum support of a rule body */
  double   conf    = 80;        /* minimum confidence of a rule */
  long     zmin    =  1;        /* minimum size of an item set */
  long     zmax    = -1;        /* maximum size of an item set */
  CCHAR    *report = "a";       /* indicators of values to report */
  CCHAR    *seval  = "x";       /* evaluation measure as a string */
  int      eval    = 'x';       /* evaluation measure */
  CCHAR    *sagg   = "x";       /* aggregation mode as a string */
  int      agg     =  0;        /* aggregation mode */
  double   thresh  = 10;        /* threshold for evaluation measure */
  CCHAR    *salgo  = "s";       /* algorithm as a string */
  int      algo    = FPG_SIMPLE;/* algorithm */
  CCHAR    *smode  = "", *s;    /* operation mode/flags as a string */
  int      mode    = FPG_DEFAULT|FPG_FIM16; /* operation mode/flags */
  long     prune   = LONG_MIN;  /* min. size for evaluation filtering */
  PyObject *border = NULL;      /* support border for filtering */
  PyObject *tracts;             /* transaction database */
  TABAG    *tabag;              /* transaction bag */
  ISREPORT *isrep;              /* item set reporter */
  REPDATA  data;                /* data for item set reporting */
  double   w;                   /* total transaction weight */
  int      r;                   /* result buffer */

  /* --- evaluate the function arguments --- */
  if (!PyArg_ParseTupleAndKeywords(args, kwds, "O|sddllsssdlssO", ckwds,
        &tracts, &starg, &supp, &conf, &zmin, &zmax, &report,
        &seval, &sagg, &thresh, &prune, &salgo, &smode, &border))
    return NULL;                /* parse the function arguments */
  target = get_target(starg, "ascmgr");
  if (target < 0) return NULL;  /* translate the target string */
  if (zmin   < 0)    { ERR_VALUE("zmin must not be negative"); }
  if (zmax   < 0)    zmax = LONG_MAX;  /* check size range */
  if (zmax   < zmin) { ERR_VALUE("zmax must not be less than zmin"); }
  eval = get_evalx(seval);      /* get the evaluation measure */
  if (eval   < 0)    return NULL;
  if (eval   <= RE_NONE) prune = LONG_MIN;
  if (strchr(smode, 'z')) eval |= FPG_INVBXS;
  agg  = get_agg(sagg);         /* get the aggregation mode */
  if (agg    < 0)    return NULL;
  thresh *= 0.01;               /* scale evaluation threshold */
  if      (strcmp(salgo, "simple")  == 0) salgo = "s";
  else if (strcmp(salgo, "complex") == 0) salgo = "c";
  else if (strcmp(salgo, "single")  == 0) salgo = "d";
  else if (strcmp(salgo, "topdown") == 0) salgo = "t";
  if (strlen(salgo) != 1)       /* translate the algorithm string */
    algo = -1;                  /* if it failed, set error code */
  else {                        /* if translation worked, */
    switch (salgo[0]) {         /* evaluate the algorithm code */
      case 's': algo = FPG_SIMPLE;  break;
      case 'c': algo = FPG_COMPLEX; break;
      case 'd': algo = FPG_SINGLE;  break;
      case 't': algo = FPG_TOPDOWN; break;
      default : algo = -1;          break;
    }                           /* set an error code for all */
  }                             /* other algorithm indicators */
  if (algo < 0) { ERR_VALUE("invalid FP-growth algorithm"); }
  for (s = smode; *s; s++) {    /* traverse the mode characters */
    if      (*s == 'l') mode &= ~FPG_FIM16;
    else if (*s == 'x') mode &= ~FPG_PERFECT;
    else if (*s == 'i') mode &= ~FPG_REORDER;
    else if (*s == 'u') mode &= ~FPG_TAIL;
  }                             /* adapt the operation mode */

  /* --- create transaction bag --- */
  tabag = tbg_fromPyObj(tracts);/* turn the given transactions */
  if (!tabag) return NULL;      /* into a transaction bag */
  conf *= 0.01;                 /* scale the minimum confidence and */
  w    = tbg_wgt(tabag);        /* get the total transaction weight */
  supp = (supp >= 0) ? 0.01 *supp *(double)w *(1-DBL_EPSILON) : -supp;
  body = (SUPP)ceilsupp(supp);  /* compute absolute support values */
  smin = (SUPP)ceilsupp(((target & ISR_RULES) && strchr(smode, 'o'))
                     ? supp : ceilsupp(supp) *conf *(1-DBL_EPSILON));
  r = fpg_data(tabag, target, smin, (ITEM)zmin, eval, algo, mode, 2);
  if (r) tbg_delete(tabag, 1);  /* prepare data for FP-growth */
  if (r == -1) { ERR_MEM(); }   /* check for error and no items */
  if (r <   0) return PyList_New(0);

  /* --- create item set reporter --- */
  isrep = isr_create(tbg_base(tabag)); /* create item set reporter */
  if (!isrep) { tbg_delete(tabag, 1); ERR_MEM(); }
  isr_setsize(isrep, (zmin > ITEM_MAX) ? ITEM_MAX : (ITEM)zmin,
                     (zmax > ITEM_MAX) ? ITEM_MAX : (ITEM)zmax);
  isr_setsupp(isrep, (RSUPP)smin, RSUPP_MAX);
  if (border && !isr_pyborder(isrep, border)) {
    isr_delete(isrep, 0); tbg_delete(tabag, 1); return NULL; }
  if ((repinit(&data, isrep, report, target) != 0)
  ||  (fpg_repo(isrep, target, eval, thresh, algo, mode) < 0)
  ||  (isr_setup(isrep) < 0)) { /* set up the item set reporter */
    isr_delete(isrep, 0); tbg_delete(tabag, 1); ERR_MEM(); }

  /* --- frequent item set mining --- */
  r = fpgrowth(tabag, target, smin, body, conf, eval, agg, thresh,
              (prune < ITEM_MIN) ? ITEM_MIN :
              (prune > ITEM_MAX) ? ITEM_MAX : (ITEM)prune,
               algo, mode, 0, isrep);
  if (r >= 0) r = repterm(&data, isrep, report);

  /* --- clean up --- */
  isr_delete(isrep, 0);         /* delete the item set reporter */
  tbg_delete(tabag, 1);         /* and the transaction bag */
  if (r < 0) { Py_DECREF(data.res); ERR_MEM(); }
  return data.res;              /* return the created result */
}  /* py_fpgrowth() */

/*--------------------------------------------------------------------*/
/* sam (tracts, target='s', supp=10, zmin=1, zmax=None, report='a',   */
/*      eval='x', thresh=10, algo='b', mode='', border=None)          */
/*--------------------------------------------------------------------*/

static PyObject* py_sam (PyObject *self,
                         PyObject *args, PyObject *kwds)
{                               /* --- SaM algorithm */
  char     *ckwds[] = { "tracts", "target", "supp", "zmin", "zmax",
                        "report", "eval", "thresh", "algo", "mode",
                        "border", NULL };
  CCHAR    *starg  = "s";       /* target type as a string */
  int      target  = ISR_SETS;  /* target type */
  double   supp    = 10;        /* minimum support of an item set */
  SUPP     smin    =  1;        /* minimum support as an integer */
  long     zmin    =  1;        /* minimum size of an item set */
  long     zmax    = -1;        /* maximum size of an item set */
  CCHAR    *report = "a";       /* indicators of values to report */
  CCHAR    *seval  = "x";       /* evaluation measure as a string */
  int      eval    = 'x';       /* evaluation measure */
  double   thresh  = 10;        /* threshold for evaluation measure */
  CCHAR    *salgo  = "b";       /* algorithm as a string */
  int      algo    = SAM_BSEARCH;  /* algorithm */
  CCHAR    *smode  = "", *s;    /* operation mode/flags as a string */
  int      mode    = SAM_DEFAULT|SAM_FIM16; /* operation mode/flags */
  PyObject *border = NULL;      /* support border for filtering */
  PyObject *tracts;             /* transaction database */
  TABAG    *tabag;              /* transaction bag */
  ISREPORT *isrep;              /* item set reporter */
  REPDATA  data;                /* data for item set reporting */
  double   w;                   /* total transaction weight */
  int      r;                   /* result buffer */

  /* --- evaluate the function arguments --- */
  if (!PyArg_ParseTupleAndKeywords(args, kwds, "O|sdllsssdssO", ckwds,
        &tracts, &starg, &supp, &zmin, &zmax, &report,
        &seval, &thresh, &salgo, &smode, &border))
    return NULL;                /* parse the function arguments */
  target = get_target(starg, "ascm");
  if (target < 0) return NULL;  /* translate the target string */
  if (zmin   < 0)    { ERR_VALUE("zmin must not be negative"); }
  if (zmax   < 0)    zmax = ITEM_MAX;  /* check size range */
  if (zmax   < zmin) { ERR_VALUE("zmax must not be less than zmin"); }
  eval = get_eval(seval);       /* get evaluation measure and */
  if (eval   < 0) return NULL;  /* check whether it is valid */
  thresh *= 0.01;               /* scale evaluation threshold */
  if      (strcmp(salgo, "basic")   == 0) salgo = "s";
  else if (strcmp(salgo, "simple")  == 0) salgo = "s";
  else if (strcmp(salgo, "bsearch") == 0) salgo = "b";
  else if (strcmp(salgo, "double")  == 0) salgo = "d";
  else if (strcmp(salgo, "tree")    == 0) salgo = "t";
  if (strlen(salgo) != 1)       /* translate the algorithm string */
    algo = -1;                  /* if it failed, set error code */
  else {                        /* if translation worked, */
    switch (salgo[0]) {         /* evaluate the algorithm code */
      case 's': algo = SAM_BASIC;   break;
      case 'b': algo = SAM_BSEARCH; break;
      case 'd': algo = SAM_DOUBLE;  break;
      case 't': algo = SAM_TREE;    break;
      default : algo = -1;         break;
    }                           /* set an error code for all */
  }                             /* other algorithm indicators */
  if (algo < 0) { ERR_VALUE("invalid SaM algorithm"); }
  for (s = smode; *s; s++) {    /* traverse the mode characters */
    if      (*s == 'l') mode &= ~SAM_FIM16;
    else if (*s == 'x') mode &= ~SAM_PERFECT;
  }                             /* adapt the operation mode */

  /* --- create transaction bag --- */
  tabag = tbg_fromPyObj(tracts);/* turn the given transactions */
  if (!tabag) return NULL;      /* into a transaction bag */
  w    = tbg_wgt(tabag);        /* get the total transaction weight */
  supp = (supp >= 0) ? 0.01 *supp *(double)w *(1-DBL_EPSILON) : -supp;
  smin = (SUPP)ceilsupp(supp);  /* compute absolute support value */
  r = sam_data(tabag, target, smin, (ITEM)zmin, 0,
               eval, algo, mode, 2);
  if (r) tbg_delete(tabag, 1);  /* prepare data for SaM */
  if (r == -1) { ERR_MEM(); }   /* check for error and no items */
  if (r <   0) return PyList_New(0);

  /* --- create item set reporter --- */
  isrep = isr_create(tbg_base(tabag)); /* create item set reporter */
  if (!isrep) { tbg_delete(tabag, 1); ERR_MEM(); }
  isr_setsize(isrep, (zmin > ITEM_MAX) ? ITEM_MAX : (ITEM)zmin,
                     (zmax > ITEM_MAX) ? ITEM_MAX : (ITEM)zmax);
  isr_setsupp(isrep, (RSUPP)smin, RSUPP_MAX);
  if (border && !isr_pyborder(isrep, border)) {
    isr_delete(isrep, 0); tbg_delete(tabag, 1); return NULL; }
  if ((repinit(&data, isrep, report, target) != 0)
  ||  (sam_repo(isrep, target, eval, thresh, algo, mode) < 0)
  ||  (isr_setup(isrep) < 0)) { /* set up the item set reporter */
    isr_delete(isrep, 0); tbg_delete(tabag, 1); ERR_MEM(); }

  /* --- frequent item set mining --- */
  r = sam(tabag, target, smin, 0.0, 0, -1.0, eval, thresh,
          algo, mode, 8192, isrep);
  if (r >= 0) r = repterm(&data, isrep, report);

  /* --- clean up --- */
  isr_delete(isrep, 0);         /* delete the item set reporter */
  tbg_delete(tabag, 1);         /* and the transaction bag */
  if (r < 0) { Py_DECREF(data.res); ERR_MEM(); }
  return data.res;              /* return the created result */
}  /* py_sam() */

/*--------------------------------------------------------------------*/
/* relim (tracts, target='s', supp=10, zmin=1, zmax=None, report='a', */
/*        eval='x', thresh=10, algo='s', mode='', border=None)        */
/*--------------------------------------------------------------------*/

static PyObject* py_relim (PyObject *self,
                           PyObject *args, PyObject *kwds)
{                               /* --- RElim algorithm */
  char     *ckwds[] = { "tracts", "target", "supp", "zmin", "zmax",
                        "report", "eval", "thresh", "algo", "mode",
                        "border", NULL };
  CCHAR    *starg  = "s";       /* target type as a string */
  int      target  = ISR_SETS;  /* target type */
  double   supp    = 10;        /* minimum support of an item set */
  SUPP     smin    =  1;        /* minimum support as an integer */
  long     zmin    =  1;        /* minimum size of an item set */
  long     zmax    = -1;        /* maximum size of an item set */
  CCHAR    *report = "a";       /* indicators of values to report */
  CCHAR    *seval  = "x";       /* evaluation measure as a string */
  int      eval    = 'x';       /* evaluation measure */
  double   thresh  = 10;        /* threshold for evaluation measure */
  CCHAR    *salgo  = "s";       /* algorithm as a string */
  int      algo    = REM_BASIC; /* algorithm */
  CCHAR    *smode  = "", *s;    /* operation mode/flags as a string */
  int      mode    = REM_DEFAULT|REM_FIM16; /* operation mode/flags */
  PyObject *border = NULL;      /* support border for filtering */
  PyObject *tracts;             /* transaction database */
  TABAG    *tabag;              /* transaction bag */
  ISREPORT *isrep;              /* item set reporter */
  REPDATA  data;                /* data for item set reporting */
  double   w;                   /* total transaction weight */
  int      r;                   /* result buffer */

  /* --- evaluate the function arguments --- */
  if (!PyArg_ParseTupleAndKeywords(args, kwds, "O|sdllsssdssO", ckwds,
        &tracts, &starg, &supp, &zmin, &zmax, &report,
        &seval, &thresh, &salgo, &smode, &border))
    return NULL;                /* parse the function arguments */
  target = get_target(starg, "ascm");
  if (target < 0) return NULL;  /* translate the target string */
  if (zmin   < 0)    { ERR_VALUE("zmin must not be negative"); }
  if (zmax   < 0)    zmax = ITEM_MAX;  /* check size range */
  if (zmax   < zmin) { ERR_VALUE("zmax must not be less than zmin"); }
  eval = get_eval(seval);       /* get evaluation measure and */
  if (eval   < 0) return NULL;  /* check whether it is valid */
  thresh *= 0.01;               /* scale evaluation threshold */
  if      (strcmp(salgo, "basic")   == 0) salgo = "s";
  else if (strcmp(salgo, "simple")  == 0) salgo = "s";
  if (strlen(salgo) != 1)       /* translate the algorithm string */
    algo = -1;                  /* if it failed, set error code */
  else {                        /* if translation worked, */
    switch (salgo[0]) {         /* evaluate the algorithm code */
      case 's': algo = REM_BASIC; break;
      default : algo = -1;        break;
    }                           /* set an error code for all */
  }                             /* other algorithm indicators */
  if (algo < 0) { ERR_VALUE("invalid RElim algorithm"); }
  for (s = smode; *s; s++) {    /* traverse the mode characters */
    if      (*s == 'l') mode &= ~REM_FIM16;
    else if (*s == 'x') mode &= ~REM_PERFECT;
  }                             /* adapt the operation mode */

  /* --- create transaction bag --- */
  tabag = tbg_fromPyObj(tracts);/* turn the given transactions */
  if (!tabag) return NULL;      /* into a transaction bag */
  w    = tbg_wgt(tabag);        /* get the total transaction weight */
  supp = (supp >= 0) ? 0.01 *supp *(double)w *(1-DBL_EPSILON) : -supp;
  smin = (SUPP)ceilsupp(supp);  /* compute absolute support value */
  r = relim_data(tabag, target, smin, (ITEM)zmin, -1.0,
                 eval, algo, mode, 2);
  if (r) tbg_delete(tabag, 1);  /* prepare data for RElim */
  if (r == -1) { ERR_MEM(); }   /* check for error and no items */
  if (r <   0) return PyList_New(0);

  /* --- create item set reporter --- */
  isrep = isr_create(tbg_base(tabag)); /* create item set reporter */
  if (!isrep) { tbg_delete(tabag, 1); ERR_MEM(); }
  isr_setsize(isrep, (zmin > ITEM_MAX) ? ITEM_MAX : (ITEM)zmin,
                     (zmax > ITEM_MAX) ? ITEM_MAX : (ITEM)zmax);
  isr_setsupp(isrep, (RSUPP)smin, RSUPP_MAX);
  if (border && !isr_pyborder(isrep, border)) {
    isr_delete(isrep, 0); tbg_delete(tabag, 1); return NULL; }
  if ((repinit(&data, isrep, report, target) != 0)
  ||  (relim_repo(isrep, target, eval, thresh, algo, mode) < 0)
  ||  (isr_setup(isrep) < 0)) { /* set up the item set reporter */
    isr_delete(isrep, 0); tbg_delete(tabag, 1); ERR_MEM(); }

  /* --- frequent item set mining --- */
  r = relim(tabag, target, smin, 0.0, 0, -1.0, eval, thresh,
            algo, mode, 32, isrep);
  if (r >= 0) r = repterm(&data, isrep, report);

  /* --- clean up --- */
  isr_delete(isrep, 0);         /* delete the item set reporter */
  tbg_delete(tabag, 1);         /* and the transaction bag */
  if (r < 0) { Py_DECREF(data.res); ERR_MEM(); }
  return data.res;              /* return the created result */
}  /* py_relim() */

/*--------------------------------------------------------------------*/
/* carpenter (tracts, target='s', supp=10, zmin=1, zmax=None,         */
/*            report='a', eval='x', thresh=10, algo='a', mode='',     */
/*            border=None)                                            */
/*--------------------------------------------------------------------*/

static PyObject* py_carpenter (PyObject *self,
                               PyObject *args, PyObject *kwds)
{                               /* --- Carpenter algorithm */
  char     *ckwds[] = { "tracts", "target", "supp", "zmin", "zmax",
                        "report", "eval", "thresh", "algo", "mode",
                        "border", NULL };
  CCHAR    *starg  = "c";       /* target type as a string */
  int      target  = ISR_CLOSED;/* target type */
  double   supp    = 10;        /* minimum support of an item set */
  SUPP     smin    =  1;        /* minimum support as an integer */
  long     zmin    =  1;        /* minimum size of an item set */
  long     zmax    = -1;        /* maximum size of an item set */
  CCHAR    *report = "a";       /* indicators of values to report */
  CCHAR    *seval  = "x";       /* evaluation measure as a string */
  int      eval    = 'x';       /* evaluation measure */
  double   thresh  = 10;        /* threshold for evaluation measure */
  CCHAR    *salgo  = "a";       /* algorithm as a string */
  int      algo    = CARP_AUTO; /* algorithm */
  CCHAR    *smode  = "", *s;    /* operation mode/flags as a string */
  int      mode    = CARP_DEFAULT; /* operation mode/flags */
  PyObject *border = NULL;      /* support border for filtering */
  PyObject *tracts;             /* transaction database */
  TABAG    *tabag;              /* transaction bag */
  ISREPORT *isrep;              /* item set reporter */
  REPDATA  data;                /* data for item set reporting */
  double   w;                   /* total transaction weight */
  int      r;                   /* result buffer */

  /* --- evaluate the function arguments --- */
  if (!PyArg_ParseTupleAndKeywords(args, kwds, "O|sdllssdssO", ckwds,
        &tracts, &starg, &supp, &zmin, &zmax, &report,
        &seval, &thresh, &salgo, &smode, &border))
    return NULL;                /* parse the function arguments */
  target = get_target(starg, "cm");
  if (target < 0) return NULL;  /* translate the target string */
  if ((target != ISR_CLOSED) && (target != IST_MAXIMAL)) {
    PyErr_SetString(PyExc_ValueError, "invalid target type");
    return NULL;                /* carpenter only supports mining */
  }                             /* closed and maximal item sets */
  if (zmin   < 0)    { ERR_VALUE("zmin must not be negative"); }
  if (zmax   < 0)    zmax = LONG_MAX;  /* check size range */
  if (zmax   < zmin) { ERR_VALUE("zmax must not be less than zmin"); }
  eval = get_eval(seval);      /* get evaluation measure and */
  if (eval   < 0) return NULL;  /* check whether it is valid */
  thresh *= 0.01;               /* scale evaluation threshold */
  if      (strcmp(salgo, "auto")    == 0) salgo = "a";
  else if (strcmp(salgo, "table")   == 0) salgo = "t";
  else if (strcmp(salgo, "table")   == 0) salgo = "t";
  else if (strcmp(salgo, "tids")    == 0) salgo = "l";
  else if (strcmp(salgo, "tidlist") == 0) salgo = "l";
  else if (strcmp(salgo, "list")    == 0) salgo = "l";
  if (strlen(salgo) != 1)       /* translate the algorithm string */
    algo = -1;                  /* if it failed, set error code */
  else {                        /* if translation worked, */
    switch (salgo[0]) {         /* evaluate the algorithm code */
      case 'a': algo = CARP_AUTO;    break;
      case 't': algo = CARP_TABLE;   break;
      case 'l': algo = CARP_TIDLIST; break;
      default : algo = -1;           break;
    }                           /* set an error code for all */
  }                             /* other algorithm indicators */
  if (algo < 0) { ERR_VALUE("invalid Carpenter algorithm"); }
  for (s = smode; *s; s++) {    /* traverse the mode characters */
    if      (*s == 'x') mode &= ~CARP_PERFECT;
    else if (*s == 'z') mode |=  CARP_FILTER;
    else if (*s == 'y') mode &= ~CARP_MAXONLY;
    else if (*s == 'p') mode &= ~CARP_COLLATE;
  }                             /* adapt the operation mode */

  /* --- create transaction bag --- */
  tabag = tbg_fromPyObj(tracts);/* turn the given transactions */
  if (!tabag) return NULL;      /* into a transaction bag */
  w    = tbg_wgt(tabag);        /* get the total transaction weight */
  supp = (supp >= 0) ? 0.01 *supp *(double)w *(1-DBL_EPSILON) : -supp;
  smin = (SUPP)ceilsupp(supp);  /* compute absolute support value */
  r = carp_data(tabag, target, smin, (ITEM)zmin, eval, algo, mode, 2);
  if (r) tbg_delete(tabag, 1);  /* prepare data for Carpenter */
  if (r == -1) { ERR_MEM(); }   /* check for error and no items */
  if (r <   0) return PyList_New(0);

  /* --- create item set reporter --- */
  isrep = isr_create(tbg_base(tabag)); /* create item set reporter */
  if (!isrep) { tbg_delete(tabag, 1); ERR_MEM(); }
  isr_setsize(isrep, (zmin > ITEM_MAX) ? ITEM_MAX : (ITEM)zmin,
                     (zmax > ITEM_MAX) ? ITEM_MAX : (ITEM)zmax);
  isr_setsupp(isrep, (RSUPP)smin, RSUPP_MAX);
  if (border && !isr_pyborder(isrep, border)) {
    isr_delete(isrep, 0); tbg_delete(tabag, 1); return NULL; }
  if ((repinit(&data, isrep, report, target) != 0)
  ||  (carp_repo(isrep, target, eval, thresh, algo, mode) < 0)
  ||  (isr_setup(isrep) < 0)) { /* set up the item set reporter */
    isr_delete(isrep, 0); tbg_delete(tabag, 1); ERR_MEM(); }

  /* --- frequent item set mining --- */
  r = carpenter(tabag, target, smin, eval, thresh, algo, mode, isrep);
  if (r >= 0) r = repterm(&data, isrep, report);

  /* --- clean up --- */
  isr_delete(isrep, 0);         /* delete the item set reporter */
  tbg_delete(tabag, 1);         /* and the transaction bag */
  if (r < 0) { Py_DECREF(data.res); ERR_MEM(); }
  return data.res;              /* return the created result */
}  /* py_carpenter() */

/*--------------------------------------------------------------------*/
/* ista (tracts, target='s', supp=10, zmin=1, zmax=None, report='a',  */
/*       eval='x', thresh=10, algo='x', mode='', border=None)         */
/*--------------------------------------------------------------------*/

static PyObject* py_ista (PyObject *self,
                          PyObject *args, PyObject *kwds)
{                               /* --- IsTa algorithm */
  char     *ckwds[] = { "tracts", "target", "supp", "zmin", "zmax",
                        "report", "eval", "thresh", "algo", "mode",
                        "border", NULL };
  CCHAR    *starg  = "c";       /* target type as a string */
  int      target  = ISR_CLOSED;/* target type */
  double   supp    = 10;        /* minimum support of an item set */
  SUPP     smin    =  1;        /* minimum support as an integer */
  long     zmin    =  1;        /* minimum size of an item set */
  long     zmax    = -1;        /* maximum size of an item set */
  CCHAR    *report = "a";       /* indicators of values to report */
  CCHAR    *seval  = "x";       /* evaluation measure as a string */
  int      eval    = 'x';       /* evaluation measure */
  double   thresh  = 10;        /* threshold for evaluation measure */
  CCHAR    *salgo  = "x";       /* algorithm as a string */
  int      algo    = ISTA_PREFIX;  /* algorithm */
  CCHAR    *smode  = "", *s;    /* operation mode/flags as a string */
  int      mode    = ISTA_DEFAULT; /* operation mode/flags */
  PyObject *border = NULL;      /* support border for filtering */
  PyObject *tracts;             /* transaction database */
  TABAG    *tabag;              /* transaction bag */
  ISREPORT *isrep;              /* item set reporter */
  REPDATA  data;                /* data for item set reporting */
  double   w;                   /* total transaction weight */
  int      r;                   /* result buffer */

  /* --- evaluate the function arguments --- */
  if (!PyArg_ParseTupleAndKeywords(args, kwds, "O|sdllssdssO", ckwds,
        &tracts, &starg, &supp, &zmin, &zmax, &report,
        &seval, &thresh, &salgo, &smode, &border))
    return NULL;                /* parse the function arguments */
  target = get_target(starg, "cm");
  if (target < 0) return NULL;  /* translate the target string */
  if ((target != ISR_CLOSED) && (target != IST_MAXIMAL)) {
    PyErr_SetString(PyExc_ValueError, "invalid target type");
    return NULL;                /* carpenter only supports mining */
  }                             /* closed and maximal item sets */
  if (zmin   < 0)    { ERR_VALUE("zmin must not be negative"); }
  if (zmax   < 0)    zmax = LONG_MAX;  /* check size range */
  if (zmax   < zmin) { ERR_VALUE("zmax must not be less than zmin"); }
  eval = get_eval(seval);       /* get evaluation measure and */
  if (eval   < 0) return NULL;  /* check whether it is valid */
  thresh *= 0.01;               /* scale evaluation threshold */
  if      (strcmp(salgo, "pfx")      == 0) salgo = "x";
  else if (strcmp(salgo, "prefix")   == 0) salgo = "x";
  else if (strcmp(salgo, "pat")      == 0) salgo = "p";
  else if (strcmp(salgo, "patricia") == 0) salgo = "p";
  if (strlen(salgo) != 1)       /* translate the algorithm string */
    algo = -1;                  /* if it failed, set error code */
  else {                        /* if translation worked, */
    switch (salgo[0]) {         /* evaluate the algorithm code */
      case 'x': algo = ISTA_PREFIX;   break;
      case 'p': algo = ISTA_PATRICIA; break;
      default : algo = -1;            break;
    }                           /* set an error code for all */
  }                             /* other algorithm indicators */
  if (algo < 0) { ERR_VALUE("invalid IsTa algorithm"); }
  for (s = smode; *s; s++) {    /* traverse the mode characters */
    if      (*s == 'p') mode &= ~ISTA_PRUNE;
    else if (*s == 'z') mode |=  ISTA_FILTER;
  }                             /* adapt the operation mode */

  /* --- create transaction bag --- */
  tabag = tbg_fromPyObj(tracts);/* turn the given transactions */
  if (!tabag) return NULL;      /* into a transaction bag */
  w    = tbg_wgt(tabag);        /* get the total transaction weight */
  supp = (supp >= 0) ? 0.01 *supp *(double)w *(1-DBL_EPSILON) : -supp;
  smin = (SUPP)ceilsupp(supp);  /* compute absolute support value */
  r = ista_data(tabag, target, smin, (ITEM)zmin, eval, algo, mode, 2);
  if (r) tbg_delete(tabag, 1);  /* prepare data for IsTa */
  if (r == -1) { ERR_MEM(); }   /* check for error and no items */
  if (r <   0) return PyList_New(0);

  /* --- create item set reporter --- */
  isrep = isr_create(tbg_base(tabag)); /* create item set reporter */
  if (!isrep) { tbg_delete(tabag, 1); ERR_MEM(); }
  isr_setsize(isrep, (zmin > ITEM_MAX) ? ITEM_MAX : (ITEM)zmin,
                     (zmax > ITEM_MAX) ? ITEM_MAX : (ITEM)zmax);
  isr_setsupp(isrep, (RSUPP)smin, RSUPP_MAX);
  if (border && !isr_pyborder(isrep, border)) {
    isr_delete(isrep, 0); tbg_delete(tabag, 1); return NULL; }
  if ((repinit(&data, isrep, report, target) != 0)
  ||  (ista_repo(isrep, target, eval, thresh, algo, mode) < 0)
  ||  (isr_setup(isrep) < 0)) { /* set up the item set reporter */
    isr_delete(isrep, 0); tbg_delete(tabag, 1); ERR_MEM(); }

  /* --- frequent item set mining --- */
  r = ista(tabag, target, smin, eval, thresh, algo, mode, isrep);
  if (r >= 0) r = repterm(&data, isrep, report);

  /* --- clean up --- */
  isr_delete(isrep, 0);         /* delete the item set reporter */
  tbg_delete(tabag, 1);         /* and the transaction bag */
  if (r < 0) { Py_DECREF(data.res); ERR_MEM(); }
  return data.res;              /* return the created result */
}  /* py_ista() */

/*--------------------------------------------------------------------*/
/* apriacc (tracts, supp=-2, zmin=2, zmax=None, report='aP',          */
/*          stat='c', siglvl=1, prune=0, mode='', border=None)        */
/*--------------------------------------------------------------------*/

static PyObject* py_apriacc (PyObject *self,
                             PyObject *args, PyObject *kwds)
{                               /* --- Apriori algorithm */
  char     *ckwds[] = { "tracts", "supp", "zmin", "zmax", "report",
                        "stat", "siglvl", "prune", "mode", "border",
                        NULL };
  double   supp    = -2;        /* minimum support of an item set */
  SUPP     smin    =  1;        /* minimum support as an integer */
  long     zmin    =  2;        /* minimum size of an item set */
  long     zmax    = -1;        /* maximum size of an item set */
  CCHAR    *report = "aP";      /* indicators of values to report */
  int      stat    =  0;        /* test statistic */
  CCHAR    *sstat  = "c";       /* test statistic as a string (chi^2) */
  double   siglvl  =  1;        /* minimum evaluation measure value */
  CCHAR    *smode  = "";        /* operation mode/flags as a string */
  int      mode    = APR_DEFAULT;  /* operation mode/flags */
  long     prune   =  0;        /* min. size for evaluation filtering */
  PyObject *border = NULL;      /* support border for filtering */
  PyObject *tracts;             /* transaction database */
  TABAG    *tabag;              /* transaction bag */
  ISREPORT *isrep;              /* item set reporter */
  REPDATA  data;                /* data for item set reporting */
  double   w;                   /* total transaction weight */
  int      r;                   /* result buffer */

  /* --- evaluate the function arguments --- */
  if (!PyArg_ParseTupleAndKeywords(args, kwds, "O|dllssdlsO", ckwds,
      &tracts, &supp, &zmin, &zmax, &report,
      &sstat, &siglvl, &prune, &smode, &border))
    return NULL;                /* parse the function arguments */
  if (zmin   < 0)    { ERR_VALUE("zmin must not be negative"); }
  if (zmax   < 0)    zmax = LONG_MAX;  /* check size range */
  if (zmax   < zmin) { ERR_VALUE("zmax must not be less than zmin"); }
  stat = get_stat(sstat);       /* translate the statistic string */
  if (stat   < 0) return NULL;  /* and check whether it is valid */
  if (siglvl <= 0)  { ERR_VALUE("siglvl must be positive"); }
  if (strchr(smode, 'z')) stat |= IST_INVBXS;
  siglvl *= 0.01;               /* scale significance level */

  /* --- create transaction bag --- */
  tabag = tbg_fromPyObj(tracts);/* turn the given transactions */
  if (!tabag) return NULL;      /* into a transaction bag */
  w    = tbg_wgt(tabag);        /* get the total transaction weight */
  supp = (supp >= 0) ? 0.01 *supp *(double)w *(1-DBL_EPSILON) : -supp;
  smin = (SUPP)ceilsupp(supp);  /* compute absolute support value */
  r = apriori_data(tabag, ISR_MAXIMAL, smin, (ITEM)zmin,
                   stat, APR_BASIC, mode, 2);
  if (r) tbg_delete(tabag, 1);  /* prepare data for Apriori */
  if (r == -1) { ERR_MEM(); }   /* check for error and no items */
  if (r <   0) return PyList_New(0);

  /* --- create item set reporter --- */
  isrep = isr_create(tbg_base(tabag)); /* create item set reporter */
  if (!isrep) { tbg_delete(tabag, 1); ERR_MEM(); }
  isr_setsize(isrep, (zmin > ITEM_MAX) ? ITEM_MAX : (ITEM)zmin,
                     (zmax > ITEM_MAX) ? ITEM_MAX : (ITEM)zmax);
  isr_setsupp(isrep, (RSUPP)smin, RSUPP_MAX);
  if (border && !isr_pyborder(isrep, border)) {
    isr_delete(isrep, 0); tbg_delete(tabag, 1); return NULL; }
  if ((repinit(&data, isrep, report, ISR_SETS) != 0)
  ||  (apriori_repo(isrep,ISR_MAXIMAL,stat,siglvl,APR_BASIC,mode) < 0)
  ||  (isr_setup(isrep) < 0)) { /* set up the item set reporter */
    isr_delete(isrep, 0); tbg_delete(tabag, 1); ERR_MEM(); }

  /* --- frequent item set mining --- */
  r = apriori(tabag, ISR_MAXIMAL, smin, smin, 1, stat, IST_MAX, siglvl,
              (prune < ITEM_MIN) ? ITEM_MIN :
              (prune > ITEM_MAX) ? ITEM_MAX : (ITEM)prune,
              APR_BASIC, mode, 0.01, 0, isrep);
  if (r >= 0) r = repterm(&data, isrep, report);

  /* --- clean up --- */
  isr_delete(isrep, 0);         /* delete the item set reporter */
  tbg_delete(tabag, 1);         /* and the transaction bag */
  if (r < 0) { Py_DECREF(data.res); ERR_MEM(); }
  return data.res;              /* return the created result */
}  /* py_apriacc() */

/*--------------------------------------------------------------------*/
/* accretion (tracts, supp=-2, zmin=2, zmax=None, report='aP',        */
/*            stat='c', siglvl=1, maxext=2, mode='', border=None)     */
/*--------------------------------------------------------------------*/

static PyObject* py_accretion (PyObject *self,
                               PyObject *args, PyObject *kwds)
{                               /* --- Accretion algorithm */
  char     *ckwds[] = { "tracts", "supp", "zmin", "zmax", "report",
                        "stat", "siglvl", "maxext", "mode", "border",
                        NULL };
  double   supp    = -2;        /* minimum support of an item set */
  SUPP     smin    =  1;        /* minimum support as an integer */
  long     zmin    =  2;        /* minimum size of an item set */
  long     zmax    = -1;        /* maximum size of an item set */
  CCHAR    *report = "aP";      /* indicators of values to report */
  int      stat    =  0;        /* test statistic */
  CCHAR    *sstat  = "c";       /* test statistic as a string (chi^2) */
  double   siglvl  =  1;        /* significance level (max. p-value) */
  CCHAR    *smode  = "";        /* operation mode/flags as a string */
  int      mode    = ACC_DEFAULT;  /* operation mode/flags */
  long     maxext  =  2;        /* maximum number of extension items */
  PyObject *border = NULL;      /* support border for filtering */
  PyObject *tracts;             /* transaction database */
  TABAG    *tabag;              /* transaction bag */
  ISREPORT *isrep;              /* item set reporter */
  REPDATA  data;                /* data for item set reporting */
  double   w;                   /* total transaction weight */
  int      r;                   /* result buffer */

  /* --- evaluate the function arguments --- */
  if (!PyArg_ParseTupleAndKeywords(args, kwds, "O|dllssdlsO", ckwds,
        &tracts, &supp, &zmin, &zmax, &report,
        &sstat, &siglvl, &maxext, &smode, &border))
    return NULL;                /* parse the function arguments */
  if (zmin   < 0)    { ERR_VALUE("zmin must not be negative"); }
  if (zmax   < 0)    zmax = LONG_MAX; /* check size range */
  if (zmax   < zmin) { ERR_VALUE("zmax must not be less than zmin"); }
  stat = get_stat(sstat);       /* translate the statistic string */
  if (stat   < 0) return NULL;  /* and check whether it is valid */
  if (strchr(smode, 'z')) stat |= ACC_INVBXS;
  if (siglvl <= 0)  { ERR_VALUE("siglvl must be positive"); }
  siglvl *= 0.01;               /* scale the significance level */
  if (maxext < 0)               /* a negative value means that */
    maxext = LONG_MAX;          /* there is no limit on extensions */

  /* --- create transaction bag --- */
  tabag = tbg_fromPyObj(tracts);/* turn the given transactions */
  if (!tabag) return NULL;      /* into a transaction bag */
  w    = tbg_wgt(tabag);        /* get the total transaction weight */
  supp = (supp >= 0) ? 0.01 *supp *(double)w *(1-DBL_EPSILON) : -supp;
  smin = (SUPP)ceilsupp(supp);  /* compute absolute support value */
  r = acc_data(tabag, ISR_MAXIMAL, smin, (ITEM)zmin, mode, 2);
  if (r) tbg_delete(tabag, 1);  /* prepare data for Accretion */
  if (r == -1) { ERR_MEM(); }   /* check for error and no items */
  if (r <   0) return PyList_New(0);

  /* --- create item set reporter --- */
  isrep = isr_create(tbg_base(tabag)); /* create item set reporter */
  if (!isrep) { tbg_delete(tabag, 1); ERR_MEM(); }
  isr_setsize(isrep, (zmin > ITEM_MAX) ? ITEM_MAX : (ITEM)zmin,
                     (zmax > ITEM_MAX) ? ITEM_MAX : (ITEM)zmax);
  isr_setsupp(isrep, (RSUPP)smin, RSUPP_MAX);
  if (border && !isr_pyborder(isrep, border)) {
    isr_delete(isrep, 0); tbg_delete(tabag, 1); return NULL; }
  if ((repinit(&data, isrep, report, ISR_SETS) != 0)
  ||  (acc_repo(isrep, ISR_MAXIMAL, mode) < 0)
  ||  (isr_setup(isrep) < 0)) { /* set up the item set reporter */
    isr_delete(isrep, 0); tbg_delete(tabag, 1); ERR_MEM(); }

  /* --- frequent item set mining --- */
  r = accretion(tabag, ISR_MAXIMAL, smin, stat, siglvl, mode,
                (maxext > ITEM_MAX) ? ITEM_MAX : (ITEM)maxext, isrep);
  if (r >= 0) r = repterm(&data, isrep, report);

  /* --- clean up --- */
  isr_delete(isrep, 0);         /* delete the item set reporter */
  tbg_delete(tabag, 1);         /* and the transaction bag */
  if (r < 0) { Py_DECREF(data.res); ERR_MEM(); }
  return data.res;              /* return the created result */
} /* py_accretion() */

/*--------------------------------------------------------------------*/
/* patspec (tracts, target='c', supp=2, zmin=2, zmax=None,            */
/*          report='#', cnt=1000, surr='p', seed=0, cpus=0)           */
/*--------------------------------------------------------------------*/

static WORKERDEF(worker, p)
{                               /* --- worker function for a thread */
  WORKDATA *w = p;              /* type the argument pointer */
  long     i;                   /* loop variable for data sets */

  assert(p);                    /* check the function argument */
  for (i = 1; i <= w->cnt; i++){/* generate cnt surrogate data sets */
    w->tasur = w->surrfn(w->tabag, w->rng, w->tasur);
    tbg_itsort(w->tasur, +1,0); /* re-sort items in transactions */
    tbg_sort  (w->tasur, +1,0); /* sort the trans. lexicographically */
    tbg_pack  (w->tasur, 16);   /* pack the most frequent items */
    w->err |= fpgrowth(w->tasur, w->target, w->smin, w->smin, 1,
                       RE_NONE, FPG_NONE, 0, 0, FPG_SIMPLE,
                       FPG_DEFAULT|FPG_FIM16, 0, w->isrep);
    if (w->err < 0) break;      /* execute the CoCoNAD algorithm */
    if (aborted)    break;      /* check for an abort interrupt */
    w->comp[0]++;               /* count the surrogate data set */
    if ((w->comp[0] % 20) == 0) /* report the data set number */
      fprintf(stderr, "%10ld\b\b\b\b\b\b\b\b\b\b", w->comp[0]);
  }
  return THREAD_OK;             /* return a dummy result */
}  /* worker() */

/*--------------------------------------------------------------------*/

static PyObject* py_patspec (PyObject *self,
                             PyObject *args, PyObject *kwds)
{                               /* --- generate a pattern spectrum */
  char      *ckwds[] = { "tracts", "target", "supp", "zmin", "zmax",
                         "report", "cnt", "surr", "seed", "cpus", NULL};
  CCHAR     *starg  = "s";      /* target type as a string */
  int       target  = ISR_CLOSED;  /* target type identifier */
  double    supp    =  2;       /* minimum support of an item set */
  SUPP      smin    =  2;       /* minimum support as an integer */
  long      zmin    =  2;       /* minimum size of an item set */
  long      zmax    = -1;       /* maximum size of an item set */
  CCHAR     *report = "#";      /* indicators of reporting format */
  long      cnt     = 1000;     /* number of data sets to generate */
  CCHAR     *ssurr  = "p";      /* surrogate method as a string */
  int       surr    = 2;        /* surrogate method identifier */
  long      seed    = 0;        /* seed for random number generator */
  int       cpus    = 0;        /* number of cpus */
  PyObject  *pypsp  = NULL;     /* created Python pattern spectrum */
  int       r;                  /* result of function call */
  int       k;                  /* loop variable for threads */
  long      i, c, x;            /* loop variable for data sets */
  double    wgt;                /* total transaction weight */
  PyObject  *tracts;            /* transaction database */
  TABAG     *tabag;             /* transaction bag (original data) */
  TABAG     *tasur, *s;         /* surrogate data set */
  TBGSURRFN *surrfn;            /* surrogate data generation function */
  RNG       *rng;               /* random number generator */
  ISREPORT  *isrep;             /* item set reporter */
  PATSPEC   *psp;               /* created pattern spectrum */
  THREAD    *threads;           /* thread handles */
  WORKDATA  *w;                 /* data for worker thread */
  volatile long comp = 0;       /* number of completed surrogates */
  #ifdef _WIN32                 /* if Microsoft Windows system */
  DWORD     thid;               /* dummy for storing the thread id */
  #endif                        /* (not really needed here) */

  /* --- evaluate the function arguments --- */
  if (!PyArg_ParseTupleAndKeywords(args, kwds,"O|sdllslsli", ckwds,
         &tracts, &starg, &supp, &zmin, &zmax, &report,
         &cnt, &ssurr, &seed, &cpus))
    return NULL;                /* parse the function arguments */
  target = get_target(starg, "ascm");
  if (target < 0) return NULL;  /* translate the target string */
  if (zmin   < 1)    { ERR_VALUE("zmin must be positive"); }
  if (zmax   < 1)    zmax = ITEM_MAX;  /* check the size range */
  if (zmax   < zmin) { ERR_VALUE("zmax must not be less than zmin"); }
  if (cnt   <= 0) cnt = 1;      /* check the number of data sets */
  surr   = get_surr(ssurr);     /* translate the surrogate string */
  if (surr  < 0)  return NULL;  /* and check for a valid code */
  if (surr == 0)  cnt = 1;      /* only one surrogate for identity */
  if (seed  == 0) seed = (long)time(NULL);

  /* --- create transaction bag --- */
  tabag = tbg_fromPyObj(tracts);/* turn the given transactions */
  if (!tabag) return NULL;      /* into a transaction bag */
  if ((surr == 3) && !tbg_istab(tabag)) {
    tbg_delete(tabag, 1);       /* if shuffle surrogates requested */
    ERR_VALUE("for shuffle surrogates transactions must form a table");
  }                             /* check for table derived data */
  wgt  = tbg_wgt(tabag);        /* get the total transaction weight */
  supp = (supp >= 0) ? 0.01 *supp *(double)wgt *(1-DBL_EPSILON) : -supp;
  smin = (SUPP)ceilsupp(supp);  /* compute absolute support value */
  r = fpg_data(tabag, target, smin, (ITEM)zmin, RE_NONE,
               FPG_SIMPLE, FPG_DEFAULT, 2);
  if (r == E_NOMEM) { tbg_delete(tabag, 1); ERR_MEM(); }

  /* --- generate pattern spectrum --- */
  if (cpus <= 0) cpus = cpucnt();
  if ((cpus > 1) && (cnt > 1)){ /* if to use multi-threading */
    threads = calloc((size_t)cpus, sizeof(THREAD));
    if (!threads) {             /* create array of thread handles */
      tbg_delete(tabag, 1); ERR_MEM(); }
    w = calloc((size_t)cpus, sizeof(WORKDATA));
    if (!w) {                   /* create array of worker data */
      free(threads); tbg_delete(tabag, 1); ERR_MEM(); }
    siginstall();               /* install the signal handler */
    c = (cnt+cpus-1) /cpus;     /* number of data sets per thread */
    for (k = 0; k < cpus; k++){ /* traverse the threads */
      x = cnt -k*c;             /* get the number of data sets and */
      if (x <= 0) continue;     /* check whether thread is needed */
      w[k].cnt    = (x < c) ? x : c;
      w[k].tabag  = tabag;      /* note the original train set */
      w[k].tasur  = tbg_clone(tabag);    /* and create a clone */
      w[k].surrfn = sur_tab[surr];
      w[k].rng    = rng_create((unsigned int)(seed+k));
      w[k].target = target;     /* create random number generator */
      w[k].smin   = smin;       /* and store the mining parameters */
      w[k].isrep  = isr_create(tbg_base(tabag));
      w[k].err    = 0;          /* create an item set reporter */
      w[k].comp   = &comp;
      if (!w[k].tasur || !w[k].rng || !w[k].isrep) {
        w[k].err = -1; break; } /* check for successful creation */
      isr_setsize(w[k].isrep,   /* set the size limits */
                  (zmin > ITEM_MAX) ? ITEM_MAX : (ITEM)zmin,
                  (zmax > ITEM_MAX) ? ITEM_MAX : (ITEM)zmax);
      isr_setsupp(w[k].isrep, (RSUPP)smin, RSUPP_MAX);
      if ((isr_addpsp(w[k].isrep, NULL) < 0)
      ||  (fpg_repo  (w[k].isrep, target, RE_NONE, 0,
                      FPG_SIMPLE, FPG_DEFAULT) < 0)
      ||  (isr_setup (w[k].isrep) != 0)) {
        w[k].err = -1; break; } /* set up the item set reporter */
      #ifdef _WIN32             /* if Microsoft Windows system */
      threads[k] = CreateThread(NULL, 0, worker, w+k, 0, &thid);
      if (!threads[k]) { w[k].err = -1; break; }
      #else                     /* if Linux/Unix system */
      if (pthread_create(threads+k, NULL, worker, w+k) != 0) {
        w[k].err = -1; break; } /* create a thread for each data set */
      #endif                    /* to compute surrogates in parallel */
    }
    #ifdef _WIN32               /* if Microsoft Windows system */
    WaitForMultipleObjects(k, threads, TRUE, INFINITE);
    while (--k >= 0)            /* wait for threads to finish, */
      CloseHandle(threads[k]);  /* then close all thread handles */
    #else                       /* if Linux/Unix system */
    while (--k >= 0)            /* wait for threads to finish */
      pthread_join(threads[k], NULL);
    #endif                      /* (join threads with this one) */
    sigremove();                /* remove the signal handler */
    for (k = r = 0; k < cpus; k++)
      r |= w[k].err;            /* join the error indicators */
    if (r >= 0) {               /* if processing was successful */
      psp = isr_rempsp(w[0].isrep, 0);
      for (k = 1; k < cpus; k++) {
        r = psp_addpsp(psp, isr_getpsp(w[k].isrep));
        if (r < 0) break;       /* traverse and sum pattern spectrums */
      }                         /* in the one of the first thread */
      if (r >= 0)               /* create a Python pattern spectrum */
        pypsp = psp_toPyObj(psp, 1/(double)cnt, report[0]);
      psp_delete(psp);          /* delete the pattern spectrum */
    }                           /* retrieved from the 1st thread */
    for (k = cpus; --k >= 0;) { /* traverse the worker data */
      if (w[k].tasur)  tbg_delete(w[k].tasur, 0);
      if (w[k].rng)    rng_delete(w[k].rng);
      if (w[k].isrep)  isr_delete(w[k].isrep, 0);
    }                           /* delete the data structures */
    free(w);                    /* delete array of worker data */
    free(threads); }            /* and array of thread handles */
  else {                        /* if to use only one thread */
    isrep = isr_create(tbg_base(tabag));
    if (!isrep) { tbg_delete(tabag, 1); ERR_MEM(); }
    isr_setsize(isrep, (zmin > ITEM_MAX) ? ITEM_MAX : (ITEM)zmin,
                       (zmax > ITEM_MAX) ? ITEM_MAX : (ITEM)zmax);
    isr_setsupp(isrep, (RSUPP)smin, RSUPP_MAX);
    if ((isr_addpsp(isrep, NULL) < 0)
    ||  (fpg_repo  (isrep, target,RE_NONE,0,FPG_SIMPLE,FPG_DEFAULT) < 0)
    ||  (isr_setup (isrep) != 0)) {
      isr_delete(isrep, 0); tbg_delete(tabag, 1); ERR_MEM(); }
    if (seed == 0) seed = (long)time(NULL);
    rng = rng_create((unsigned int)seed);
    if (!rng) {                 /* create a random number generator */
      isr_delete(isrep, 0); tbg_delete(tabag, 1); ERR_MEM(); }
    siginstall();               /* install the signal handler */
    surrfn = sur_tab[surr];     /* get the surrogate data function */
    tasur  = NULL; r = 0;       /* init. surrogate and return code */
    for (i = 1; i <= cnt; i++){ /* generate cnt surrogate data sets */
      s = surrfn(tabag, rng, tasur);
      if (!s) { r = -1; break;} /* generate next surrogate data set */
      tasur = s;                /* note the created data set */
      tbg_itsort(tasur, +1, 0); /* re-sort items in transactions */
      tbg_sort  (tasur, +1, 0); /* sort the trans. lexicographically */
      tbg_pack  (tasur, 16);    /* pack the most frequent items */
      r = fpgrowth(tasur, target, smin, smin, 1, RE_NONE, FPG_NONE,
                   0, 0, FPG_SIMPLE, FPG_DEFAULT, 0, isrep);
      if (r < 0)   break;       /* execute the FP-growth algorithm */
      if ((i % 20) == 0) fprintf(stderr, "%10ld\b\b\b\b\b\b\b\b\b\b", i);
      if (aborted) break;       /* report the data set number */
    }                           /* and check for an interrupt */
    sigremove();                /* remove the signal handler */
    if (r >= 0)                 /* create a Python pattern spectrum */
      pypsp = psp_toPyObj(isr_getpsp(isrep), 1/(double)cnt, report[0]);
    if (tasur) tbg_delete(tasur, 0);
    rng_delete(rng);            /* delete the data structures */
    isr_delete(isrep,  0);
  }

  /* --- clean up --- */
  if (aborted) { PyErr_SetInterrupt(); aborted = 0; }
  tbg_delete(tabag, 1);         /* delete the created train set */
  if ((r < 0) || !pypsp) { ERR_MEM(); }
  return pypsp;                 /* return created pattern spectrum */
}  /* py_patspec() */

/*--------------------------------------------------------------------*/
/* estpsp (tracts, target='s', supp=2, zmin=2, zmax=None, report='#', */
/*         equiv=10000, alpha=0.5, smpls=1000, seed=0)                */
/*--------------------------------------------------------------------*/

static PyObject* py_estpsp (PyObject *self,
                            PyObject *args, PyObject *kwds)
{                               /* --- estimate a pattern spectrum */
  char    *ckwds[] = { "tracts", "target", "supp", "zmin", "zmax",
                       "report", "equiv", "alpha", "smpls", "seed",
                       NULL };
  long     equiv   = 10000;     /* equivalent number of surrogates */
  CCHAR    *starg  = "s";       /* target type as a string */
  int      target  = ISR_SETS;  /* target type identifier */
  double   supp    =  2;        /* minimum support of an item set */
  SUPP     smin    =  2;        /* minimum support as an integer */
  double   alpha   =  0.5;      /* probability dispersion factor */
  long     smpls   = 1000;      /* number of samples per set size */
  long     zmin    =  2;        /* minimum size of an item set */
  long     zmax    = -1;        /* maximum size of an item set */
  long     seed    =  0;        /* seed for random number generator */
  CCHAR    *report = "#";       /* indicators of reporting format */
  PyObject *pypsp  = NULL;      /* created Python pattern spectrum */
  PyObject *tracts;             /* transaction database (Python) */
  TABAG    *tabag;              /* transaction database (C) */
  PATSPEC  *psp;                /* created pattern spectrum */
  double   wgt;                 /* total transaction weight */
  int      r;                   /* result of function call */

  /* --- evaluate the function arguments --- */
  if (!PyArg_ParseTupleAndKeywords(args,kwds,"O|dllldlls", ckwds,
        &tracts, &starg, &supp, &zmin, &zmax, &report,
        &equiv, &alpha, &smpls, &seed))
    return NULL;                /* parse the function arguments */
  target = get_target(starg, "as");
  if (target < 0) return NULL;  /* translate the target string */
  if (zmin  <  1)    { ERR_VALUE("zmin must be positive"); }
  if (zmax  <  1)    zmax = ITEM_MAX;  /* check the size range */
  if (zmax  <  zmin) { ERR_VALUE("zmax must not be less than zmin"); }
  if (equiv <= 0)    equiv = 1; /* check the number of data sets */
  if (smpls <= 0)    { ERR_VALUE("smpls must be positive"); }
  if (seed  == 0)    seed = (long)time(NULL);

  /* --- create transaction bag --- */
  tabag = tbg_fromPyObj(tracts);/* turn the given transactions */
  if (!tabag) return NULL;      /* into a transaction bag */
  wgt  = tbg_wgt(tabag);        /* get the total transaction weight */
  supp = (supp >= 0) ? 0.01 *supp *(double)wgt *(1-DBL_EPSILON) : -supp;
  smin = (SUPP)ceilsupp(supp);  /* compute absolute support value */
  r = tbg_recode(tabag, smin, -1, -1, -2);
  if (r < 0) { tbg_delete(tabag, 1); ERR_MEM(); }
  tbg_filter(tabag, (ITEM)zmin, NULL, 0);

  /* --- estimate pattern spectrum --- */
  rseed((unsigned)seed);        /* seed random number generator */
  psp = psp_create((ITEM)zmin, (ITEM)zmax, smin, tbg_cnt(tabag));
  r = (!psp) ? -1               /* estimate a pattern spectrum */
    : psp_tbgest(tabag, psp, (size_t)equiv, alpha, (size_t)smpls);
  if (!r)                       /* turn it into a python object */
    pypsp = psp_toPyObj(psp, 1/(double)equiv, report[0]);

  /* --- clean up --- */
  psp_delete(psp);              /* delete the pattern spectrum */
  tbg_delete(tabag, 1);         /* and the train set */
  if (r) { ERR_MEM(); }         /* check for an error and return */
  return pypsp;                 /* the created pattern spectrum */
}  /* py_estpsp() */

/*--------------------------------------------------------------------*/
/* Python Function List                                               */
/*--------------------------------------------------------------------*/

static PyMethodDef fim_methods[] = {
  { "fim", (PyCFunction)py_fim, METH_VARARGS|METH_KEYWORDS,
    "fim (tracts, target='s', supp=10, zmin=1, zmax=None, report='a',\n"
    "     eval='x', agg='x', thresh=10, border=None)\n"
    "Find frequent item sets (simplified interface).\n"
    "tracts  transaction database to mine (mandatory)\n"
    "        The database must be an iterable of transactions;\n"
    "        each transaction must be an iterable of items;\n"
    "        each item must be a hashable object.\n"
    "        If the database is a dictionary, the transactions are\n"
    "        the keys, the values their (integer) multiplicities.\n"
    "target  type of frequent item sets to find     (default: s)\n"
    "        s/a   sets/all   all     frequent item sets\n"
    "        c     closed     closed  frequent item sets\n"
    "        m     maximal    maximal frequent item sets\n"
    "        g     gens       generators\n"
    "supp    minimum support of an item set         (default: 10)\n"
    "        (positive: percentage, negative: absolute number)\n"
    "zmin    minimum number of items per item set   (default: 1)\n"
    "zmax    maximum number of items per item set   (default: no limit)\n"
    "report  values to report with an item set      (default: a)\n"
    "        a     absolute item set support (number of transactions)\n"
    "        s     relative item set support as a fraction\n"
    "        S     relative item set support as a percentage\n"
    "        e     value of item set evaluation measure\n"
    "        E     value of item set evaluation measure as a percentage\n"
    "        =     pattern spectrum as a list (instead of patterns)\n"
    "        #     pattern spectrum as a dictionary\n"
    "eval    measure for item set evaluation        (default: x)\n"
    "        x     none       no measure / zero (default)\n"
    "        b     ldratio    binary logarithm of support quotient       (+)\n"
    "        c     conf       rule confidence                            (+)\n"
    "        d     confdiff   absolute confidence difference to prior    (+)\n"
    "        l     lift       lift value (confidence divided by prior)   (+)\n"
    "        a     liftdiff   absolute difference of lift value to 1     (+)\n"
    "        q     liftquot   difference of lift quotient to 1           (+)\n"
    "        v     cvct       conviction (inverse lift for negated head) (+)\n"
    "        e     cvctdiff   absolute difference of conviction to 1     (+)\n"
    "        r     cvctquot   difference of conviction quotient to 1     (+)\n"
    "        k     cprob      conditional probability ratio              (+)\n"
    "        j     import     importance (binary log. of prob. ratio)    (+)\n"
    "        z     cert       certainty factor (relative conf. change)   (+)\n"
    "        n     chi2       normalized chi^2 measure                   (+)\n"
    "        p     chi2pval   p-value from (unnormalized) chi^2 measure  (-)\n"
    "        y     yates      normalized chi^2 with Yates' correction    (+)\n"
    "        t     yatespval  p-value from Yates-corrected chi^2 measure (-)\n"
    "        i     info       information difference to prior            (+)\n"
    "        g     infopval   p-value from G statistic/info. difference  (-)\n"
    "        f     fetprob    Fisher's exact test (table probability)    (-)\n"
    "        h     fetchi2    Fisher's exact test (chi^2 measure)        (-)\n"
    "        m     fetinfo    Fisher's exact test (mutual information)   (-)\n"
    "        s     fetsupp    Fisher's exact test (support)              (-)\n"
    "        Measures marked with (+) must meet or exceed the threshold,\n"
    "        measures marked with (-) must not exceed the threshold\n"
    "        in order for the item set to be reported.\n"
    "agg     evaluation measure aggregation mode    (default: x)\n"
    "        x     none       no aggregation (use first value)\n"
    "        m     min        minimum of individual measure values\n"
    "        n     max        maximum of individual measure values\n"
    "        a     avg        average of individual measure values\n"
    "thresh  threshold for evaluation measure       (default: 10%)\n"
    "border  support border for filtering item sets (default: None)\n"
    "        Must be a list or tuple of (absolute) minimum support values\n"
    "        per item set size (by which the list/tuple is indexed).\n"
    "returns a list of pairs (i.e. tuples with two elements),\n"
    "        each consisting of a tuple with a found frequent item set\n"
    "        and a tuple listing the values selected with 'report' *or*\n"
    "        a list of triplets (size,supp,frq), i.e. a pattern spectrum."
  },
  { "arules", (PyCFunction)py_arules, METH_VARARGS|METH_KEYWORDS,
    "arules (tracts, supp=10, conf=80, zmin=1, zmax=None, report='aC',\n"
    "        eval='x', thresh=10, mode='')\n"
    "Find association rules (simplified interface).\n"
    "tracts  transaction database to mine (mandatory)\n"
    "        The database must be an iterable of transactions;\n"
    "        each transaction must be an iterable of items;\n"
    "        each item must be a hashable object.\n"
    "        If the database is a dictionary, the transactions are\n"
    "        the keys, the values their (integer) multiplicities.\n"
    "supp    minimum support    of an assoc. rule   (default: 10)\n"
    "        (positive: percentage, negative: absolute number)\n"
    "conf    minimum confidence of an assoc. rule   (default: 80%)\n"
    "zmin    minimum number of items per rule       (default: 1)\n"
    "zmax    maximum number of items per rule       (default: no limit)\n"
    "report  values to report with an item set      (default: aC)\n"
    "        a     absolute item set  support (number of transactions)\n"
    "        s     relative item set  support as a fraction\n"
    "        S     relative item set  support as a percentage\n"
    "        b     absolute body set  support (number of transactions)\n"
    "        x     relative body set  support as a fraction\n"
    "        X     relative body set  support as a percentage\n"
    "        h     absolute head item support (number of transactions)\n"
    "        y     relative head item support as a fraction\n"
    "        Y     relative head item support as a percentage\n"
    "        c     rule confidence as a fraction\n"
    "        C     rule confidence as a percentage\n"
    "        l     lift value of a rule (confidence/prior)\n"
    "        L     lift value of a rule as a percentage\n"
    "        e     value of rule evaluation measure\n"
    "        E     value of rule evaluation measure as a percentage\n"
    "        =     pattern spectrum as a list (instead of patterns)\n"
    "        #     pattern spectrum as a dictionary\n"
    "eval    measure for rule evaluation            (default: x)\n"
    "        x     none       no measure / zero (default)\n"
    "        b     ldratio    binary logarithm of support quotient       (+)\n"
    "        c     conf       rule confidence                            (+)\n"
    "        d     confdiff   absolute confidence difference to prior    (+)\n"
    "        l     lift       lift value (confidence divided by prior)   (+)\n"
    "        a     liftdiff   absolute difference of lift value to 1     (+)\n"
    "        q     liftquot   difference of lift quotient to 1           (+)\n"
    "        v     cvct       conviction (inverse lift for negated head) (+)\n"
    "        e     cvctdiff   absolute difference of conviction to 1     (+)\n"
    "        r     cvctquot   difference of conviction quotient to 1     (+)\n"
    "        k     cprob      conditional probability ratio              (+)\n"
    "        j     import     importance (binary log. of prob. ratio)    (+)\n"
    "        z     cert       certainty factor (relative conf. change)   (+)\n"
    "        n     chi2       normalized chi^2 measure                   (+)\n"
    "        p     chi2pval   p-value from (unnormalized) chi^2 measure  (-)\n"
    "        y     yates      normalized chi^2 with Yates' correction    (+)\n"
    "        t     yatespval  p-value from Yates-corrected chi^2 measure (-)\n"
    "        i     info       information difference to prior            (+)\n"
    "        g     infopval   p-value from G statistic/info. difference  (-)\n"
    "        f     fetprob    Fisher's exact test (table probability)    (-)\n"
    "        h     fetchi2    Fisher's exact test (chi^2 measure)        (-)\n"
    "        m     fetinfo    Fisher's exact test (mutual information)   (-)\n"
    "        s     fetsupp    Fisher's exact test (support)              (-)\n"
    "        Measures marked with (+) must meet or exceed the threshold,\n"
    "        measures marked with (-) must not exceed the threshold\n"
    "        in order for the item set to be reported.\n"
    "thresh  threshold for evaluation measure       (default: 10%)\n"
    "mode    operation mode indicators/flags        (default: None)\n"
    "        o     use original rule support definition (body & head)\n"
    "returns a list of triplets (i.e. tuples with three elements),\n"
    "        each consisting of a head/consequent item, a tuple with\n"
    "        a body/antecedent item set, and a tuple listing the values\n"
    "        selected with the parameter 'report'."
  },
  { "apriori", (PyCFunction)py_apriori, METH_VARARGS|METH_KEYWORDS,
    "apriori (tracts, target='s', supp=10, zmin=1, zmax=None, report='a',\n"
    "         eval='x', agg='x', thresh=10, prune=None, algo='b', mode='',\n"
    "         border=None)\n"
    "Find frequent item sets with the Apriori algorithm.\n"
    "tracts  transaction database to mine (mandatory)\n"
    "        The database must be an iterable of transactions;\n"
    "        each transaction must be an iterable of items;\n"
    "        each item must be a hashable object.\n"
    "        If the database is a dictionary, the transactions are\n"
    "        the keys, the values their (integer) multiplicities.\n"
    "target  type of frequent item sets to find     (default: s)\n"
    "        s/a   sets/all   all     frequent item sets\n"
    "        c     closed     closed  frequent item sets\n"
    "        m     maximal    maximal frequent item sets\n"
    "        g     gens       generators\n"
    "        r     rules      association rules\n"
    "supp    minimum support of an item set         (default: 10)\n"
    "        (positive: percentage, negative: absolute number)\n"
    "conf    minimum confidence of an assoc. rule   (default: 80%)\n"
    "zmin    minimum number of items per item set   (default: 1)\n"
    "zmax    maximum number of items per item set   (default: no limit)\n"
    "report  values to report with an item set      (default: a)\n"
    "        a     absolute item set support (number of transactions)\n"
    "        s     relative item set support as a fraction\n"
    "        S     relative item set support as a percentage\n"
    "        e     value of item set evaluation measure\n"
    "        E     value of item set evaluation measure as a percentage\n"
    "        =     pattern spectrum as a list (instead of patterns)\n"
    "        #     pattern spectrum as a dictionary\n"
    "eval    measure for item set evaluation        (default: x)\n"
    "        x     none       no measure / zero (default)\n"
    "        b     ldratio    binary logarithm of support quotient       (+)\n"
    "        c     conf       rule confidence                            (+)\n"
    "        d     confdiff   absolute confidence difference to prior    (+)\n"
    "        l     lift       lift value (confidence divided by prior)   (+)\n"
    "        a     liftdiff   absolute difference of lift value to 1     (+)\n"
    "        q     liftquot   difference of lift quotient to 1           (+)\n"
    "        v     cvct       conviction (inverse lift for negated head) (+)\n"
    "        e     cvctdiff   absolute difference of conviction to 1     (+)\n"
    "        r     cvctquot   difference of conviction quotient to 1     (+)\n"
    "        k     cprob      conditional probability ratio              (+)\n"
    "        j     import     importance (binary log. of prob. ratio)    (+)\n"
    "        z     cert       certainty factor (relative conf. change)   (+)\n"
    "        n     chi2       normalized chi^2 measure                   (+)\n"
    "        p     chi2pval   p-value from (unnormalized) chi^2 measure  (-)\n"
    "        y     yates      normalized chi^2 with Yates' correction    (+)\n"
    "        t     yatespval  p-value from Yates-corrected chi^2 measure (-)\n"
    "        i     info       information difference to prior            (+)\n"
    "        g     infopval   p-value from G statistic/info. difference  (-)\n"
    "        f     fetprob    Fisher's exact test (table probability)    (-)\n"
    "        h     fetchi2    Fisher's exact test (chi^2 measure)        (-)\n"
    "        m     fetinfo    Fisher's exact test (mutual information)   (-)\n"
    "        s     fetsupp    Fisher's exact test (support)              (-)\n"
    "        Measures marked with (+) must meet or exceed the threshold,\n"
    "        measures marked with (-) must not exceed the threshold\n"
    "        in order for the item set to be reported.\n"
    "agg     evaluation measure aggregation mode    (default: x)\n"
    "        x     none       no aggregation (use first value)\n"
    "        m     min        minimum of individual measure values\n"
    "        n     max        maximum of individual measure values\n"
    "        a     avg        average of individual measure values\n"
    "thresh  threshold for evaluation measure       (default: 10%)\n"
    "prune   min. size for evaluation filtering     (default: no pruning)\n"
    "        = 0   backward filtering       (no subset check)\n"
    "        < 0   weak   forward filtering (one subset  must qualify)\n"
    "        > 0   strong forward filtering (all subsets must qualify)\n"
    "algo    algorithm variant to use               (default: a)\n"
    "        b     basic      standard algorithm (only choice)\n"
    "mode    operation mode indicators/flags        (default: None)\n"
    "        x     do not use perfect extension pruning\n"
    "        t/T   do not organize transactions as a prefix tree\n"
    "        y     a-posteriori pruning of infrequent item sets\n"
    "        z     invalidate evaluation below expected support\n"
    "        o     use original rule support definition (body & head)\n"
    "border  support border for filtering item sets (default: None)\n"
    "        Must be a list or tuple of (absolute) minimum support values\n"
    "        per item set size (by which the list/tuple is indexed).\n"
    "returns if report is not in ['#','=']:\n"
    "          if the target is association rules:\n"
    "            a list of triplets (i.e. tuples with three elements),\n"
    "            each consisting of a head/consequent item, a tuple\n"
    "            with a body/antecedent item set, and a tuple listing\n"
    "            the values selected with the parameter 'report'.\n"
    "          if the target is a type of item sets:\n"
    "            a list of pairs (i.e. tuples with two elements), each\n"
    "            consisting of a tuple with a found frequent item set\n"
    "            and a tuple listing the values selected with 'report'.\n"
    "        if report in ['#','=']:\n"
    "          a list of triplets (size,supp,freq), "
              "i.e. a pattern spectrum."
  },
  { "eclat", (PyCFunction)py_eclat, METH_VARARGS|METH_KEYWORDS,
    "eclat (tracts, target='s', supp=10, zmin=1, zmax=None, report='a',\n"
    "       eval='x', agg='x', thresh=10, prune=None, algo='a', mode='',\n"
    "       border=None)\n"
    "Find frequent item sets with the Eclat algorithm.\n"
    "tracts  transaction database to mine (mandatory)\n"
    "        The database must be an iterable of transactions;\n"
    "        each transaction must be an iterable of items;\n"
    "        each item must be a hashable object.\n"
    "        If the database is a dictionary, the transactions are\n"
    "        the keys, the values their (integer) multiplicities.\n"
    "target  type of frequent item sets to find     (default: s)\n"
    "        s/a   sets/all   all     frequent item sets\n"
    "        c     closed     closed  frequent item sets\n"
    "        m     maximal    maximal frequent item sets\n"
    "        g     gens       generators\n"
    "        r     rules      association rules\n"
    "supp    minimum support of an item set         (default: 10)\n"
    "        (positive: percentage, negative: absolute number)\n"
    "conf    minimum confidence of an assoc. rule   (default: 80%)\n"
    "zmin    minimum number of items per item set   (default: 1)\n"
    "zmax    maximum number of items per item set   (default: no limit)\n"
    "report  values to report with an item set      (default: a)\n"
    "        a     absolute item set support (number of transactions)\n"
    "        s     relative item set support as a fraction\n"
    "        S     relative item set support as a percentage\n"
    "        e     value of item set evaluation measure\n"
    "        E     value of item set evaluation measure as a percentage\n"
    "        =     pattern spectrum as a list (instead of patterns)\n"
    "        #     pattern spectrum as a dictionary\n"
    "eval    measure for item set evaluation        (default: x)\n"
    "        x     none       no measure / zero (default)\n"
    "        b     ldratio    binary logarithm of support quotient       (+)\n"
    "        c     conf       rule confidence                            (+)\n"
    "        d     confdiff   absolute confidence difference to prior    (+)\n"
    "        l     lift       lift value (confidence divided by prior)   (+)\n"
    "        a     liftdiff   absolute difference of lift value to 1     (+)\n"
    "        q     liftquot   difference of lift quotient to 1           (+)\n"
    "        v     cvct       conviction (inverse lift for negated head) (+)\n"
    "        e     cvctdiff   absolute difference of conviction to 1     (+)\n"
    "        r     cvctquot   difference of conviction quotient to 1     (+)\n"
    "        k     cprob      conditional probability ratio              (+)\n"
    "        j     import     importance (binary log. of prob. ratio)    (+)\n"
    "        z     cert       certainty factor (relative conf. change)   (+)\n"
    "        n     chi2       normalized chi^2 measure                   (+)\n"
    "        p     chi2pval   p-value from (unnormalized) chi^2 measure  (-)\n"
    "        y     yates      normalized chi^2 with Yates' correction    (+)\n"
    "        t     yatespval  p-value from Yates-corrected chi^2 measure (-)\n"
    "        i     info       information difference to prior            (+)\n"
    "        g     infopval   p-value from G statistic/info. difference  (-)\n"
    "        f     fetprob    Fisher's exact test (table probability)    (-)\n"
    "        h     fetchi2    Fisher's exact test (chi^2 measure)        (-)\n"
    "        m     fetinfo    Fisher's exact test (mutual information)   (-)\n"
    "        s     fetsupp    Fisher's exact test (support)              (-)\n"
    "        Measures marked with (+) must meet or exceed the threshold,\n"
    "        measures marked with (-) must not exceed the threshold\n"
    "        in order for the item set to be reported.\n"
    "agg     evaluation measure aggregation mode    (default: x)\n"
    "        x     none       no aggregation (use first value)\n"
    "        m     min        minimum of individual measure values\n"
    "        n     max        maximum of individual measure values\n"
    "        a     avg        average of individual measure values\n"
    "thresh  threshold for evaluation measure       (default: 10%)\n"
    "prune   min. size for evaluation filtering     (default: no pruning)\n"
    "        = 0   backward filtering       (no subset check)\n"
    "        < 0   weak   forward filtering (one subset  must qualify)\n"
    "        > 0   strong forward filtering (all subsets must qualify)\n"
    "algo    algorithm variant to use               (default: a)\n"
    "        a     auto       automatic choice based on data properties\n"
    "        e     basic      transaction id lists intersection (basic)\n"
    "        i     tids       transaction id lists intersection (improved)\n"
    "        b     bits       transaction id lists as bit vectors\n"
    "        t     table      item occurrence table (standard)\n"
    "        s     simple     item occurrence table (simplified)\n"
    "        r     ranges     transaction id range lists intersection\n"
    "        o     occdlv     occurrence deliver from transaction lists\n"
    "        d     diff       transaction id difference sets (diffsets)\n"
    "mode    operation mode indicators/flags        (default: None)\n"
    "        l     do not use a 16-items machine\n"
    "        x     do not use perfect extension pruning\n"
    "        i     do not sort items w.r.t. conditional support\n"
    "        u     do not head union tail (hut) pruning (maximal)\n"
    "        y     check extensions for closed/maximal item sets\n"
    "        z     invalidate evaluation below expected support\n"
    "        o     use original rule support definition (body & head)\n"
    "border  support border for filtering item sets (default: None)\n"
    "        Must be a list or tuple of (absolute) minimum support values\n"
    "        per item set size (by which the list/tuple is indexed).\n"
    "returns if report is not in ['#','=']:\n"
    "          if the target is association rules:\n"
    "            a list of triplets (i.e. tuples with three elements),\n"
    "            each consisting of a head/consequent item, a tuple\n"
    "            with a body/antecedent item set, and a tuple listing\n"
    "            the values selected with the parameter 'report'.\n"
    "          if the target is a type of item sets:\n"
    "            a list of pairs (i.e. tuples with two elements), each\n"
    "            consisting of a tuple with a found frequent item set\n"
    "            and a tuple listing the values selected with 'report'.\n"
    "        if report in ['#','=']:\n"
    "          a list of triplets (size,supp,freq), "
              "i.e. a pattern spectrum."
  },
  { "fpgrowth", (PyCFunction)py_fpgrowth, METH_VARARGS|METH_KEYWORDS,
    "fpgrowth (tracts, target='s', supp=10, zmin=1, zmax=None, report='a',\n"
    "          eval='x', agg='x', thresh=10, prune=Nobe, algo='s', mode='',\n"
    "          border=None)\n"
    "Find frequent item sets with the FP-growth algorithm.\n"
    "tracts  transaction database to mine (mandatory)\n"
    "        The database must be an iterable of transactions;\n"
    "        each transaction must be an iterable of items;\n"
    "        each item must be a hashable object.\n"
    "        If the database is a dictionary, the transactions are\n"
    "        the keys, the values their (integer) multiplicities.\n"
    "target  type of frequent item sets to find     (default: s)\n"
    "        s/a   sets/all   all     frequent item sets\n"
    "        c     closed     closed  frequent item sets\n"
    "        m     maximal    maximal frequent item sets\n"
    "        g     gens       generators\n"
    "        r     rules      association rules\n"
    "supp    minimum support of an item set         (default: 10)\n"
    "        (positive: percentage, negative: absolute number)\n"
    "conf    minimum confidence of an assoc. rule   (default: 80%)\n"
    "zmin    minimum number of items per item set   (default: 1)\n"
    "zmax    maximum number of items per item set   (default: no limit)\n"
    "report  values to report with an item set      (default: a)\n"
    "        a     absolute item set support (number of transactions)\n"
    "        s     relative item set support as a fraction\n"
    "        S     relative item set support as a percentage\n"
    "        e     value of item set evaluation measure\n"
    "        E     value of item set evaluation measure as a percentage\n"
    "        =     pattern spectrum as a list (instead of patterns)\n"
    "        #     pattern spectrum as a dictionary\n"
    "eval    measure for item set evaluation        (default: x)\n"
    "        x     none       no measure / zero (default)\n"
    "        b     ldratio    binary logarithm of support quotient       (+)\n"
    "        c     conf       rule confidence                            (+)\n"
    "        d     confdiff   absolute confidence difference to prior    (+)\n"
    "        l     lift       lift value (confidence divided by prior)   (+)\n"
    "        a     liftdiff   absolute difference of lift value to 1     (+)\n"
    "        q     liftquot   difference of lift quotient to 1           (+)\n"
    "        v     cvct       conviction (inverse lift for negated head) (+)\n"
    "        e     cvctdiff   absolute difference of conviction to 1     (+)\n"
    "        r     cvctquot   difference of conviction quotient to 1     (+)\n"
    "        k     cprob      conditional probability ratio              (+)\n"
    "        j     import     importance (binary log. of prob. ratio)    (+)\n"
    "        z     cert       certainty factor (relative conf. change)   (+)\n"
    "        n     chi2       normalized chi^2 measure                   (+)\n"
    "        p     chi2pval   p-value from (unnormalized) chi^2 measure  (-)\n"
    "        y     yates      normalized chi^2 with Yates' correction    (+)\n"
    "        t     yatespval  p-value from Yates-corrected chi^2 measure (-)\n"
    "        i     info       information difference to prior            (+)\n"
    "        g     infopval   p-value from G statistic/info. difference  (-)\n"
    "        f     fetprob    Fisher's exact test (table probability)    (-)\n"
    "        h     fetchi2    Fisher's exact test (chi^2 measure)        (-)\n"
    "        m     fetinfo    Fisher's exact test (mutual information)   (-)\n"
    "        s     fetsupp    Fisher's exact test (support)              (-)\n"
    "        Measures marked with (+) must meet or exceed the threshold,\n"
    "        measures marked with (-) must not exceed the threshold\n"
    "        in order for the item set to be reported.\n"
    "agg     evaluation measure aggregation mode    (default: x)\n"
    "        x     none       no aggregation (use first value)\n"
    "        m     min        minimum of individual measure values\n"
    "        n     max        maximum of individual measure values\n"
    "        a     avg        average of individual measure values\n"
    "thresh  threshold for evaluation measure       (default: 10%)\n"
    "prune   min. size for evaluation filtering     (default: no pruning)\n"
    "        = 0   backward filtering       (no subset check)\n"
    "        < 0   weak   forward filtering (one subset  must qualify)\n"
    "        > 0   strong forward filtering (all subsets must qualify)\n"
    "algo    algorithm variant to use               (default: s)\n"
    "        s     simple     simple  tree nodes (only link and parent)\n"
    "        c     complex    complex tree nodes (children and siblings)\n"
    "        d     single     top-down processing on a single prefix tree\n"
    "        t     topdown    top-down processing of the prefix trees\n"
    "        Variant d does not support closed/maximal item set mining.\n"
    "mode    operation mode indicators/flags        (default: None)\n"
    "        l     do not use a 16-items machine\n"
    "        x     do not use perfect extension pruning\n"
    "        i     do not sort items w.r.t. conditional support\n"
    "        u     do not head union tail (hut) pruning (maximal)\n"
    "        z     invalidate evaluation below expected support\n"
    "        o     use original rule support definition (body & head)\n"
    "border  support border for filtering item sets (default: None)\n"
    "        Must be a list or tuple of (absolute) minimum support values\n"
    "        per item set size (by which the list/tuple is indexed).\n"
    "returns if report is not in ['#','=']:\n"
    "          if the target is association rules:\n"
    "            a list of triplets (i.e. tuples with three elements),\n"
    "            each consisting of a head/consequent item, a tuple\n"
    "            with a body/antecedent item set, and a tuple listing\n"
    "            the values selected with the parameter 'report'.\n"
    "          if the target is a type of item sets:\n"
    "            a list of pairs (i.e. tuples with two elements), each\n"
    "            consisting of a tuple with a found frequent item set\n"
    "            and a tuple listing the values selected with 'report'.\n"
    "        if report in ['#','=']:\n"
    "          a list of triplets (size,supp,freq), "
              "i.e. a pattern spectrum."
  },
  { "sam", (PyCFunction)py_sam, METH_VARARGS|METH_KEYWORDS,
    "sam (tracts, target='s', supp=10, zmin=1, zmax=None, report='a',\n"
    "     eval='x', thresh=10, algo='b', mode='', border=None)\n"
    "Find frequent item sets with the SaM algorithm.\n"
    "tracts  transaction database to mine (mandatory)\n"
    "        The database must be an iterable of transactions;\n"
    "        each transaction must be an iterable of items;\n"
    "        each item must be a hashable object.\n"
    "        If the database is a dictionary, the transactions are\n"
    "        the keys, the values their (integer) multiplicities.\n"
    "target  type of frequent item sets to find     (default: s)\n"
    "        s/a   sets/all   all     frequent item sets\n"
    "        c     closed     closed  frequent item sets\n"
    "        m     maximal    maximal frequent item sets\n"
    "supp    minimum support of an item set         (default: 10)\n"
    "        (positive: percentage, negative: absolute number)\n"
    "zmin    minimum number of items per item set   (default: 1)\n"
    "zmax    maximum number of items per item set   (default: no limit)\n"
    "report  values to report with an item set      (default: a)\n"
    "        a     absolute item set support (number of transactions)\n"
    "        s     relative item set support as a fraction\n"
    "        S     relative item set support as a percentage\n"
    "        e     value of item set evaluation measure\n"
    "        E     value of item set evaluation measure as a percentage\n"
    "        =     pattern spectrum as a list (instead of patterns)\n"
    "        #     pattern spectrum as a dictionary\n"
    "eval    measure for item set evaluation        (default: x)\n"
    "        x     none       no measure / zero (default)\n"
    "        b     ldratio    binary logarithm of support quotient\n"
    "thresh  threshold for evaluation measure       (default: 10%)\n"
    "algo    algorithm variant to use               (default: o)\n"
    "        s     simple     basic split and merge algorithm\n"
    "        b     bsearch    split and merge with binary search\n"
    "        d     double     SaM with double source buffering\n"
    "        t     tree       SaM with transaction prefix tree\n"
    "mode    operation mode indicators/flags        (default: None)\n"
    "        l     do not use a 16-items machine\n"
    "        x     do not use perfect extension pruning\n"
    "border  support border for filtering item sets (default: None)\n"
    "        Must be a list or tuple of (absolute) minimum support values\n"
    "        per item set size (by which the list/tuple is indexed).\n"
    "returns a list of pairs (i.e. tuples with two elements),\n"
    "        each consisting of a tuple with a found frequent item set\n"
    "        and a tuple listing the values selected with 'report' *or*\n"
    "        a list of triplets (size,supp,frq), i.e. a pattern spectrum."
  },
  { "relim", (PyCFunction)py_relim, METH_VARARGS|METH_KEYWORDS,
    "relim (tracts, target='s', supp=10, zmin=1, zmax=None, report='a',\n"
    "       eval='x', thresh=10, algo='s', mode='', border=None)\n"
    "Find frequent item sets with the RElim algorithm.\n"
    "tracts  transaction database to mine (mandatory)\n"
    "        The database must be an iterable of transactions;\n"
    "        each transaction must be an iterable of items;\n"
    "        each item must be a hashable object.\n"
    "        If the database is a dictionary, the transactions are\n"
    "        the keys, the values their (integer) multiplicities.\n"
    "target  type of frequent item sets to find     (default: s)\n"
    "        s/a   sets/all   all     frequent item sets\n"
    "        c     closed     closed  frequent item sets\n"
    "        m     maximal    maximal frequent item sets\n"
    "supp    minimum support of an item set         (default: 10)\n"
    "        (positive: percentage, negative: absolute number)\n"
    "zmin    minimum number of items per item set   (default: 1)\n"
    "zmax    maximum number of items per item set   (default: no limit)\n"
    "report  values to report with an item set      (default: a)\n"
    "        a     absolute item set support (number of transactions)\n"
    "        s     relative item set support as a fraction\n"
    "        S     relative item set support as a percentage\n"
    "        e     value of item set evaluation measure\n"
    "        E     value of item set evaluation measure as a percentage\n"
    "        =     pattern spectrum as a list (instead of patterns)\n"
    "        #     pattern spectrum as a dictionary\n"
    "eval    measure for item set evaluation        (default: x)\n"
    "        x     none       no measure / zero (default)\n"
    "        b     ldratio    binary logarithm of support quotient\n"
    "thresh  threshold for evaluation measure       (default: 10%)\n"
    "algo    algorithm variant to use               (default: o)\n"
    "        s     simple     basic recursive elimination algorithm\n"
    "        (this parameter is essentially a placeholder for extensions)\n"
    "mode    operation mode indicators/flags        (default: None)\n"
    "        l     do not use a 16-items machine\n"
    "        x     do not use perfect extension pruning\n"
    "border  support border for filtering item sets (default: None)\n"
    "        Must be a list or tuple of (absolute) minimum support values\n"
    "        per item set size (by which the list/tuple is indexed).\n"
    "returns a list of pairs (i.e. tuples with two elements),\n"
    "        each consisting of a tuple with a found frequent item set\n"
    "        and a tuple listing the values selected with 'report' *or*\n"
    "        a list of triplets (size,supp,frq), i.e. a pattern spectrum."
  },
  { "carpenter", (PyCFunction)py_carpenter, METH_VARARGS|METH_KEYWORDS,
    "carpenter (tracts, target='s', supp=10, zmin=1, zmax=None, report='a',\n"
    "           eval='x', thresh=10, algo='a', mode='', border=None)\n"
    "Find frequent item sets with the Carpenter algorithm.\n"
    "tracts  transaction database to mine (mandatory)\n"
    "        The database must be an iterable of transactions;\n"
    "        each transaction must be an iterable of items;\n"
    "        each item must be a hashable object.\n"
    "        If the database is a dictionary, the transactions are\n"
    "        the keys, the values their (integer) multiplicities.\n"
    "target  type of frequent item sets to find     (default: s)\n"
    "        c     closed     closed  frequent item sets\n"
    "        m     maximal    maximal frequent item sets\n"
    "supp    minimum support of an item set         (default: 10)\n"
    "        (positive: percentage, negative: absolute number)\n"
    "zmin    minimum number of items per item set   (default: 1)\n"
    "zmax    maximum number of items per item set   (default: no limit)\n"
    "report  values to report with an item set      (default: a)\n"
    "        a     absolute item set support (number of transactions)\n"
    "        s     relative item set support as a fraction\n"
    "        S     relative item set support as a percentage\n"
    "        e     value of item set evaluation measure\n"
    "        E     value of item set evaluation measure as a percentage\n"
    "        =     pattern spectrum as a list (instead of patterns)\n"
    "        #     pattern spectrum as a dictionary\n"
    "eval    measure for item set evaluation        (default: x)\n"
    "        x     none       no measure / zero (default)\n"
    "        b     ldratio    binary logarithm of support quotient\n"
    "thresh  threshold for evaluation measure       (default: 10%)\n"
    "algo    algorithm variant to use               (default: s)\n"
    "        a     auto       automatic choice based on table size\n"
    "        t     table      item occurrence counter table\n"
    "        l     tidlist    transaction identifier lists\n"
    "mode    operation mode indicators/flags        (default: None)\n"
    "        x     do not use perfect extension pruning\n"
    "        z     filter maximal item sets with repository\n"
    "        y     add only maximal item sets to repository\n"
    "        p     do not collate equal transactions\n"
    "border  support border for filtering item sets (default: None)\n"
    "        Must be a list or tuple of (absolute) minimum support values\n"
    "        per item set size (by which the list/tuple is indexed).\n"
    "returns a list of pairs (i.e. tuples with two elements),\n"
    "        each consisting of a tuple with a found frequent item set\n"
    "        and a tuple listing the values selected with 'report' *or*\n"
    "        a list of triplets (size,supp,frq), i.e. a pattern spectrum."
  },
  { "ista", (PyCFunction)py_ista, METH_VARARGS|METH_KEYWORDS,
    "ista (tracts, target='s', supp=10, zmin=1, zmax=None, report='a',\n"
    "      eval='x', thresh=10, algo='x', mode='', border=None)\n"
    "Find frequent item sets with the IsTa algorithm.\n"
    "tracts  transaction database to mine (mandatory)\n"
    "        The database must be an iterable of transactions;\n"
    "        each transaction must be an iterable of items;\n"
    "        each item must be a hashable object.\n"
    "        If the database is a dictionary, the transactions are\n"
    "        the keys, the values their (integer) multiplicities.\n"
    "target  type of frequent item sets to find     (default: s)\n"
    "        c     closed     closed  frequent item sets\n"
    "        m     maximal    maximal frequent item sets\n"
    "supp    minimum support of an item set         (default: 10)\n"
    "        (positive: percentage, negative: absolute number)\n"
    "zmin    minimum number of items per item set   (default: 1)\n"
    "zmax    maximum number of items per item set   (default: no limit)\n"
    "report  values to report with an item set      (default: a)\n"
    "        a     absolute item set support (number of transactions)\n"
    "        s     relative item set support as a fraction\n"
    "        S     relative item set support as a percentage\n"
    "        e     value of item set evaluation measure\n"
    "        E     value of item set evaluation measure as a percentage\n"
    "        =     pattern spectrum as a list (instead of patterns)\n"
    "        #     pattern spectrum as a dictionary\n"
    "eval    measure for item set evaluation        (default: x)\n"
    "        x     none       no measure / zero (default)\n"
    "        b     ldratio    binary logarithm of support quotient\n"
    "thresh  threshold for evaluation measure       (default: 10%)\n"
    "algo    algorithm variant to use               (default: x)\n"
    "        x     prefix     use a standard prefix tree\n"
    "        p     patricia   use a patricia tree\n"
    "        (a patricia tree may be faster for very many items\n"
    "        and very few transactions)\n"
    "mode    operation mode indicators/flags        (default: None)\n"
    "        p     do not prune the prefix/patricia tree\n"
    "        z     filter maximal item sets with repository\n"
    "border  support border for filtering item sets (default: None)\n"
    "        Must be a list or tuple of (absolute) minimum support values\n"
    "        per item set size (by which the list/tuple is indexed).\n"
    "returns a list of pairs (i.e. tuples with two elements),\n"
    "        each consisting of a tuple with a found frequent item set\n"
    "        and a tuple listing the values selected with 'report' *or*\n"
    "        a list of triplets (size,supp,frq), i.e. a pattern spectrum."
  },
  { "apriacc", (PyCFunction)py_apriacc, METH_VARARGS|METH_KEYWORDS,
    "apriacc (tracts, supp=-2, zmin=2, zmax=None, report='aP',\n"
    "         stat='c', siglvl=1, prune=0, mode='', border=None)\n"
    "Find frequent item sets with an accretion-style apriori algorithm\n"
    "(that is, with an interface analog to accretion()).\n"
    "tracts  transaction database to mine           (mandatory)\n"
    "        The database must be an iterable of transactions;\n"
    "        each transaction must be an iterable of items;\n"
    "        each item must be a hashable object.\n"
    "        If the database is a dictionary, the transactions are\n"
    "        the keys, the values their (integer) multiplicities.\n"
    "supp    minimum support of an item set         (default: -2)\n"
    "        (positive: percentage, negative: absolute number)\n"
    "zmin    minimum number of items per item set   (default: 2)\n"
    "zmax    maximum number of items per item set   (default: no limit)\n"
    "report  values to report with an item set      (default: aP)\n"
    "        a     absolute item set support (number of transactions)\n"
    "        s     relative item set support as a fraction\n"
    "        S     relative item set support as a percentage\n"
    "        p     p-value of item set test as a fraction\n"
    "        P     p-value of item set test as a percentage\n"
    "        =     pattern spectrum as a list (instead of patterns)\n"
    "        #     pattern spectrum as a dictionary\n"
    "stat    test statistic for item set evaluation (default: c)\n"
    "        x     none     no statistic / zero\n"
    "        c/n/p chi2     chi^2 measure (default)\n"
    "        y/t   yates    chi^2 measure with Yates' correction\n"
    "        i/g   info     mutual information / G statistic\n"
    "        f     fetprob  Fisher's exact test (table probability)\n"
    "        h     fetchi2  Fisher's exact test (chi^2 measure)\n"
    "        m     fetinfo  Fisher's exact test (mutual information)\n"
    "        s     fetsupp  Fisher's exact test (support)\n"
    "siglvl  significance level (maximum p-value)   (default: 1%)\n"
    "prune   min. size for evaluation filtering     (default: 0)\n"
    "        = 0   backward filtering       (no subset checks)\n"
    "        < 0   weak   forward filtering (one subset  must qualify)\n"
    "        > 0   strong forward filtering (all subsets must qualify)\n"
    "mode    operation mode indicators/flags        (default: None)\n"
    "        z     invalidate evaluation below expected support\n"
    "border  support border for filtering item sets (default: None)\n"
    "        Must be a list or tuple of (absolute) minimum support values\n"
    "        per item set size (by which the list/tuple is indexed).\n"
    "returns a list of pairs (i.e. tuples with two elements),\n"
    "        each consisting of a tuple with a found frequent item set\n"
    "        and a tuple listing the values selected with 'report' *or*\n"
    "        a list of triplets (size,supp,frq), i.e. a pattern spectrum."
  },
  { "accretion", (PyCFunction)py_accretion, METH_VARARGS|METH_KEYWORDS,
    "accretion (tracts, supp=-2, zmin=2, zmax=None, report='aP',\n"
    "           stat='c', siglvl=1, maxext=2, mode='', border=None)\n"
    "Find frequent item sets with the accretion algorithm.\n"
    "tracts  transaction database to mine           (mandatory)\n"
    "        The database must be an iterable of transactions;\n"
    "        each transaction must be an iterable of items;\n"
    "        each item must be a hashable object.\n"
    "        If the database is a dictionary, the transactions are\n"
    "        the keys, the values their (integer) multiplicities.\n"
    "supp    minimum support of an item set         (default: -2)\n"
    "        (positive: percentage, negative: absolute number)\n"
    "zmin    minimum number of items per item set   (default: 2)\n"
    "zmax    maximum number of items per item set   (default: no limit)\n"
    "report  values to report with an item set      (default: aP)\n"
    "        a     absolute item set support (number of transactions)\n"
    "        s     relative item set support as a fraction\n"
    "        S     relative item set support as a percentage\n"
    "        p     p-value of item set test as a fraction\n"
    "        P     p-value of item set test as a percentage\n"
    "        =     pattern spectrum as a list (instead of patterns)\n"
    "        #     pattern spectrum as a dictionary\n"
    "stat    test statistic for item set evaluation (default: c)\n"
    "        x     none     no statistic / zero\n"
    "        c/p/n chi2     chi^2 measure (default)\n"
    "        y/t   yates    chi^2 measure with Yates' correction\n"
    "        i/g   info     mutual information / G statistic\n"
    "        f     fetprob  Fisher's exact test (table probability)\n"
    "        h     fetchi2  Fisher's exact test (chi^2 measure)\n"
    "        m     fetinfo  Fisher's exact test (mutual information)\n"
    "        s     fetsupp  Fisher's exact test (support)\n"
    "siglvl  significance level (maximum p-value)   (default: 1%)\n"
    "maxext  maximum number of extension items      (default: 2)\n"
    "mode    operation mode indicators/flags        (default: None)\n"
    "        z     invalidate evaluation below expected support\n"
    "border  support border for filtering item sets (default: None)\n"
    "        Must be a list or tuple of (absolute) minimum support values\n"
    "        per item set size (by which the list/tuple is indexed).\n"
    "returns a list of pairs (i.e. tuples with two elements),\n"
    "        each consisting of a tuple with a found frequent item set\n"
    "        and a tuple listing the values selected with 'report' *or*\n"
    "        a list of triplets (size,supp,frq), i.e. a pattern spectrum."
  },
  { "patspec", (PyCFunction)py_patspec, METH_VARARGS|METH_KEYWORDS,
    "patspec (tracts, target='s', supp=2, zmin=2, zmax=None,\n"
    "         report='#', cnt=1000, surr='p', seed=0, cpus=0)\n"
    "Generate a pattern spectrum from surrogate data sets.\n"
    "tracts  transaction database to mine           (mandatory)\n"
    "        The database must be an iterable of transactions;\n"
    "        each transaction must be an iterable of items;\n"
    "        each item must be a hashable object.\n"
    "        If the database is a dictionary, the transactions are\n"
    "        the keys, the values their (integer) multiplicities.\n"
    "target  type of frequent item sets to find     (default: s)\n"
    "        s/a  sets/all   all     frequent item sets\n"
    "        c    closed     closed  frequent item sets\n"
    "        m    maximal    maximal frequent item sets\n"
    "supp    minimum support of an item set         (default: 2)\n"
    "zmin    minimum number of items per item set   (default: 2)\n"
    "zmax    maximum number of items per item set   (default: no limit)\n"
    "report  pattern spectrum reporting format      (default: #)\n"
    "        =    pattern spectrum as a list of triplets\n"
    "        #    pattern spectrum as a dictionary\n"
    "cnt     number of surrogate data sets          (default: 1000)\n"
    "surr    surrogate data generation method       (default: p)\n"
    "        i    ident      identity (keep original data)\n"
    "        r    random     random transaction generation\n"
    "        p    swap       permutation by pair swaps\n"
    "        s    shuffle    shuffle table-derived data (columns)\n"
    "seed    seed for random number generator       (default: 0)\n"
    "        (seed = 0: use system time as a seed)\n"
    "cpus    number of cpus to use                  (default: 0)\n"
    "        A value <= 0 means all cpus reported as available.\n"
    "returns a pattern spectrum as a dictionary mapping pairs\n"
    "        (size, support) to the corresponding occurrence counters\n"
    "        or as a list of triplets (size, support, count)"
  },
  { "estpsp", (PyCFunction)py_estpsp, METH_VARARGS|METH_KEYWORDS,
    "estpsp (tracts, target='s', supp=2, zmin=2, zmax=None,\n"
    "        report='#', equiv=10000, alpha=0.5, smpls=1000, seed=0)\n"
    "Estimate a pattern spectrum from data characteristics.\n"
    "tracts  transaction database to mine           (mandatory)\n"
    "        The database must be an iterable of transactions;\n"
    "        each transaction must be an iterable of items;\n"
    "        each item must be a hashable object.\n"
    "        If the database is a dictionary, the transactions are\n"
    "        the keys, the values their (integer) multiplicities.\n"
    "target  type of frequent item sets to find     (default: s)\n"
    "        s/a  sets/all   all     frequent item sets\n"
    "supp    minimum support of an item set         (default: 2)\n"
    "zmin    minimum number of items per item set   (default: 2)\n"
    "zmax    maximum number of items per item set   (default: no limit)\n"
    "report  pattern spectrum reporting format      (default: #)\n"
    "        =    pattern spectrum as a list of triplets\n"
    "        #    pattern spectrum as a dictionary\n"
    "equiv   equivalent number of surrogates        (default: 10000)\n"
    "alpha   probability dispersion factor          (default: 0.5)\n"
    "smpls   number of samples per item set size    (default: 1000)\n"
    "seed    seed for random number generator       (default: 0)\n"
    "        (seed = 0: use system time as a seed)\n"
    "returns a pattern spectrum as a dictionary mapping pairs\n"
    "        (size, support) to the corresponding occurrence counters\n"
    "        or as a list of triplets (size, support, count)"
  },
  { NULL }                      /* sentinel */
};

/*----------------------------------------------------------------------
  Initialization Function
----------------------------------------------------------------------*/
#define FIM_DESC \
  "Frequent Item Set Mining and Association Rule Induction for Python\n" \
  "version 6.9 (2015.01.03)      (c) 2011-2015   Christian Borgelt"

#if PY_MAJOR_VERSION >= 3

static struct PyModuleDef fimdef = {
  PyModuleDef_HEAD_INIT, "fim", FIM_DESC,
  -1, fim_methods, NULL, NULL, NULL, NULL
};

PyObject* PyInit_fim (void)
{ return PyModule_Create(&fimdef); }

#else

PyMODINIT_FUNC initfim (void)
{ Py_InitModule3("fim", fim_methods, FIM_DESC); }

#endif
