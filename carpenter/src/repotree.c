/*----------------------------------------------------------------------
  File    : repotree.c
  Contents: item set repository tree management
  Author  : Christian Borgelt
  History : 2009.10.08 file created as clomax.c
            2009.10.09 item order direction added
            2010.03.11 function rpt_add() improved (create new nodes)
            2010.06.21 generalized by introducing definition of SUPP
            2010.07.04 function rpt_add() reports when tree was changed
            2010.07.05 function rpt_report() added (closed item sets)
            2010.08.18 top-level repository nodes made a fixed array
            2010.12.07 added some explicit type casts (for C++)
            2011.11.24 item order used for traversal in rpt_report()
            2012.04.10 bug in functions rep_clo() fixed (maximal sets)
            2012.04.26 special maximal item set functions added
            2012.04.27 function rpt_prune() added (support pruning)
            2013.04.01 adapted to type changes in module tract
            2013.10.15 checks of return code of isr_report() added
----------------------------------------------------------------------*/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <limits.h>
#include <assert.h>
#include "arrays.h"
#include "repotree.h"
#ifdef STORAGE
#include "storage.h"
#endif

/*----------------------------------------------------------------------
  Preprocessor Definitions
----------------------------------------------------------------------*/
#define pos(x,y)   ((x) < (y))  /* macros for item comparison */
#define neg(x,y)   ((x) > (y))  /* (ascending and descending) */

/*----------------------------------------------------------------------
  Closed Item Set Repository Tree Functions
----------------------------------------------------------------------*/

REPOTREE* rpt_create (MEMSYS *mem, ITEM size, int dir)
{                               /* --- create a repository tree */
  REPOTREE *rpt;                /* created item set repository tree */
  REPONODE *node;               /* to traverse the top-level nodes */

  assert(size >= 0);            /* check the function arguments */
  rpt = (REPOTREE*)malloc(sizeof(REPOTREE)
                        +(size_t)(size-1) *sizeof(REPONODE));
  if (!rpt) return NULL;        /* create an item set repository tree */
  rpt->size = size;             /* and initialize its fields */
  rpt->dir  = (dir < 0) ? -1 : +1;
  rpt->supp = 0;                /* create a memory system for nodes */
  rpt->mem  = (mem) ? mem : ms_create(sizeof(REPONODE), 65535);
  if (!rpt->mem) { free(rpt); return NULL; }
  while (--size >= 0) {         /* initialize the top-level nodes */
    node = rpt->tops +size; node->supp = 0; node->item = size;
    node->sibling = node->children = NULL;
  }                             /* (sibling pointer is no used) */
  return rpt;                   /* return the created prefix tree */
}  /* rpt_create() */

/*--------------------------------------------------------------------*/

static void delete (REPONODE *node, MEMSYS *mem)
{                               /* --- recursively delete nodes */
  REPONODE *tmp;                /* buffer for deallocation */

  assert(mem);                  /* check the function arguments */
  while (node) {                /* sibling list deletion loop */
    delete(node->children,mem); /* recursively delete children */
    tmp = node; node = node->sibling; ms_free(mem, tmp);
  }                             /* finally delete the node itself */
}  /* delete() */

/*--------------------------------------------------------------------*/

void rpt_delete (REPOTREE *rpt, int delms)
{                               /* --- delete a repository tree */
  assert(rpt);                  /* check the function arguments */
  if (delms)                    /* if possible/requested, */
    ms_delete(rpt->mem);        /* delete the memory system */
  else                          /* if to keep the memory system, */
    while (--rpt->size >= 0)    /* delete tree nodes recursively */
      delete(rpt->tops[rpt->size].children, rpt->mem);
  free(rpt);                    /* delete the base structure */
}  /* rpt_delete() */

/*--------------------------------------------------------------------*/

int rpt_add (REPOTREE *rpt, const ITEM *items, ITEM n, SUPP supp)
{                               /* --- add item set to repository */
  int      c = 0;               /* changed flag */
  ITEM     i;                   /* buffer for an item */
  REPONODE *node;               /* to traverse the (new) nodes */
  REPONODE **p;                 /* insertion position for new node */

  assert(rpt                    /* check the function arguments */
  &&    (items || (n <= 0)) && (supp >= 0));
  if (supp > rpt->supp) {       /* adapt the empty set support */
    rpt->supp = supp; c = 1; }  /* and set the changed flag */
  if (--n < 0) return c;        /* if there are no items, abort */
  node = rpt->tops +*items++;   /* get top-level node for first item */
  do {                          /* traverse the items of the set */
    if (supp > node->supp) {    /* adapt the item set support */
      node->supp = supp; c = 1;}/* and set the changed flag */
    if (--n < 0) return c;      /* if all items are processed, abort */
    i = *items++;               /* get the next item in the set */
    p = &node->children;        /* traverse the list of children */
    if (rpt->dir < 0) while (*p && ((*p)->item > i)) p = &(*p)->sibling;
    else              while (*p && ((*p)->item < i)) p = &(*p)->sibling;
    node = *p;                  /* find the item/insertion position */
  } while (node && (node->item == i));
  node = (REPONODE*)ms_alloc(rpt->mem);
  if (!node) return -1;         /* create a new prefix tree node */
  node->item    = i;            /* store the current item and */
  node->supp    = supp;         /* the support of the item set */
  node->sibling = *p;           /* insert the created node */
  *p = node;                    /* into the sibling list */
  while (--n >= 0) {            /* traverse the rest of the items */
    node = node->children = (REPONODE*)ms_alloc(rpt->mem);
    if (!node) return -1;       /* create a new prefix tree node */
    node->item    = *items++;   /* store the current item and */
    node->supp    = supp;       /* the support of the item set */
    node->sibling = NULL;       /* there are no siblings yet */
  }
  node->children = NULL;        /* last created node is a leaf */
  return 1;                     /* return that the tree was changed */
}  /* rpt_add() */

/*--------------------------------------------------------------------*/

SUPP rpt_get (REPOTREE *rpt, const ITEM *items, ITEM n)
{                               /* --- get support of an item set */
  ITEM     i;                   /* buffer for an item */
  REPONODE *p;                  /* to traverse the nodes */

  assert(rpt && (items || (n <= 0))); /* check function arguments */
  if (--n < 0)                  /* if there are no items, */
    return rpt->supp;           /* return the empty set support */
  p = rpt->tops +*items++;      /* get top-level node for first item */
  while (--n >= 0) {            /* while not at the last item */
    p = p->children;            /* continue with the child nodes */
    i = *items++;               /* try to find a corresp. child node */
    if (rpt->dir < 0) while (p && (p->item > i)) p = p->sibling;
    else              while (p && (p->item < i)) p = p->sibling;
    if (!p || (p->item != i))   /* if a node with the next item */
      return -1;                /* does not exist in the tree, */
  }                             /* abort the search with failure */
  return p->supp;               /* return support of the item set */
}  /* rpt_get() */

/*--------------------------------------------------------------------*/

#define SUPER(dir) \
static int super_##dir (REPONODE *node, const ITEM *items, ITEM n,     \
                        SUPP supp)                                     \
{                               /* --- check for a superset */         \
  assert(items && (n > 0) && (supp > 0));   /* check arguments */      \
  while (node                   /* while there is another node */      \
  && !dir(*items, node->item)){ /* with item before next to match */   \
    if (node->item == *items) { /* if at node with matching item */    \
      if (--n <= 0) return node->supp >= supp;                         \
      items++; }                /* check for last and skip item */     \
    else if (super_##dir(node->sibling, items, n, supp))               \
      return -1;                /* process siblings of current node */ \
    if (node->supp < supp) return 0; /* check for a frequent set */    \
    node = node->children;      /* continue with the child nodes */    \
  }                             /* (match the remaining items) */      \
  return 0;                     /* return 'no superset exists' */      \
}  /* super() */

/*--------------------------------------------------------------------*/

SUPER(pos)                      /* function for ascending  item order */
SUPER(neg)                      /* function for descending item order */

/*--------------------------------------------------------------------*/

static int super (REPOTREE *rpt, const ITEM *items, ITEM n, SUPP supp)
{                               /* --- check for a superset */
  ITEM     i, k;                /* loop variables */
  REPONODE *node;               /* to traverse the top-level nodes */

  assert(rpt                    /* check the function arguments */
  &&     items && (n > 0) && (supp > 0));
  node = rpt->tops +*items;     /* get top-level node for first item */
  if (n <= 1) {                 /* check special case for one item */
    if (node->supp >= supp) return -1; }
  else if ((rpt->dir < 0)       /* check with matching first item */
  ?        super_neg(node->children, items+1, n-1, supp)
  :        super_pos(node->children, items+1, n-1, supp))
    return -1;                  /* abort if a superset was found */
  k = (rpt->dir < 0) ? rpt->size : -1;
  for (i = *items; (i -= rpt->dir) != k; ) {
    node = rpt->tops +i;        /* traverse the top-level nodes */
    if ((rpt->dir < 0)          /* check by skipping node item */
    ?   super_neg(node->children, items, n, supp)
    :   super_pos(node->children, items, n, supp))
      return -1;                /* if a superset could be found, */
  }                             /* abort with success, otherwise */
  return 0;                     /* return 'no superset exists' */
}  /* super() */

/*--------------------------------------------------------------------*/

int rpt_super (REPOTREE *rpt, const ITEM *items, ITEM n, SUPP supp)
{                               /* --- check for a superset */
  assert(rpt                    /* check the function arguments */
  &&    (items || (n <= 0)) && (supp > 0));
  if (n <= 0)                   /* if no item, check root support */
    return (rpt->supp >= supp) ? -1 : 0;
  return super(rpt, items, n, supp);
}  /* rpt_super() */            /* recursively check for a superset */

/*--------------------------------------------------------------------*/

static void prune (REPONODE **node, SUPP supp, MEMSYS *mem)
{                               /* --- recursively prune tree */
  REPONODE *t;                  /* temporary buffer for deletion */

  while (*node) {               /* traverse the sibling list */
    if ((*node)->children)      /* prune children recursively */
      prune(&(*node)->children, supp, mem);
    if ((*node)->supp >= supp){ /* keep nodes with sufficient supp. */
      node = &(*node)->sibling; continue; }
    t = *node; *node = (*node)->sibling;
    ms_free(mem, t);            /* nodes with insufficient support */
  }                             /* are removed from sibling list */
}  /* prune() */

/*--------------------------------------------------------------------*/

void rpt_prune (REPOTREE *rpt, SUPP supp)
{                               /* --- prune item set repository */
  ITEM i;                       /* loop variable */
  assert(rpt && (supp >= 0));   /* check the function arguments */
  for (i = rpt->size; --i >= 0;)/* recursively prune sibling lists */
    prune(&rpt->tops[i].children, supp, rpt->mem);
}  /* rpt_prune() */

/*--------------------------------------------------------------------*/

static int closed (REPOTREE *rpt, REPONODE *node)
{                               /* --- report item sets recursively */
  int  r;                       /* error status */
  int  x = 0;                   /* flag for a perfect extension */
  SUPP supp;                    /* support of current item set */

  assert(rpt && node);          /* check the function arguments */
  supp = isr_supp(rpt->rep);    /* get current item set support */
  if (isr_xable(rpt->rep, 1)) { /* if item set may be extended */
    for (node = node->children; node; node = node->sibling) {
      if (node->supp < rpt->min)/* traverse the child nodes, */
        continue;               /* but skip infrequent item sets */
      x |= (node->supp >= supp);/* set perfect extension flag */
      r  = isr_addnc(rpt->rep, node->item, node->supp);
      if (r < 0) return r;      /* add current item to the reporter */
      r  = closed(rpt, node);   /* recursively report item sets, */
      isr_remove(rpt->rep, 1);  /* then remove the current item */
      if (r < 0) return r;      /* from the item set reporter; */
    } }                         /* finally check for an error */
  else {                        /* if item set cannot be extended */
    for (node = node->children ; node; node = node->sibling)
      if (node->supp >= supp) { x = -1; break; }
  }                             /* check for a perfect extension */
  return (x) ? 0 : isr_report(rpt->rep);
}  /* closed() */               /* report sets without perfect exts. */

/*--------------------------------------------------------------------*/

static int maximal (REPOTREE *rpt, REPONODE *node)
{                               /* --- report item sets recursively */
  int r;                        /* error status */
  int x = 0;                    /* flag for a frequent child */

  assert(rpt && node);          /* check the function arguments */
  if (isr_xable(rpt->rep, 1)) { /* if item set may be extended */
    for (node = node->children; node; node = node->sibling) {
      if (node->supp < rpt->min)/* traverse the child nodes, */
        continue;               /* but skip infrequent item sets */
      x = -1;                   /* set flag for a frequent child */
      r = isr_addnc(rpt->rep, node->item, node->supp);
      if (r < 0) return r;      /* add current item to the reporter */
      r = maximal(rpt, node);   /* report the current item sets, */
      isr_remove(rpt->rep, 1);  /* then remove the current item */
      if (r < 0) return r;      /* from the item set reporter; */
    } }                         /* finally check for an error */
  else {                        /* if item set may not be extended */
    for (node = node->children; node; node = node->sibling)
      if (node->supp >= rpt->min) { x = -1; break; }
  }                             /* check for a frequent superset */
  return (x) ? 0 : isr_report(rpt->rep);
}  /* maximal() */              /* report sets w/o frequent supersets */

/*--------------------------------------------------------------------*/

static int maxonly (REPOTREE *rpt, REPONODE *node)
{                               /* --- report item sets recursively */
  int      r;                   /* error status */
  int      x = 0;               /* flag for a frequent child */
  REPONODE *curr = node;        /* buffer for current node */

  assert(rpt && node);          /* check the function arguments */
  if (isr_xable(rpt->rep, 1)) { /* if item set may be extended */
    for (node = node->children; node; node = node->sibling) {
      if (node->supp < rpt->min)/* traverse the child nodes, */
        continue;               /* but skip infrequent item sets */
      x = -1;                   /* set flag for a frequent child */
      r = isr_addnc(rpt->rep, node->item, node->supp);
      if (r < 0) return r;      /* add current item to the reporter */
      r = maxonly(rpt, node);   /* recursively report supersets */
      isr_remove(rpt->rep, 1);  /* report the current item sets, */
      if (r < 0) return r;      /* then remove the current item */
    } }                         /* from the item set reporter */
  else {                        /* if item set may not be extended */
    for (node = node->children; node; node = node->sibling)
      if (node->supp >= rpt->min) { x = -1; break; }
  }                             /* check for a frequent superset */
  if (x) return 0;              /* if there is a freq. child, abort */
  curr->supp = -curr->supp;     /* mark the current node as to ignore */
  x = super(rpt, isr_items(rpt->rep), isr_cnt(rpt->rep), rpt->min);
  curr->supp = -curr->supp;     /* unmark the current node */
  return (x) ? 0 : isr_report(rpt->rep);
}  /* maxonly() */              /* report sets w/o frequent supersets */

/*--------------------------------------------------------------------*/

int rpt_report (REPOTREE *rpt, int max, SUPP supp, ISREPORT *rep)
{                               /* --- report sets in repository */
  int      r;                   /* error status */
  int      x = 0;               /* flag for an extension */
  ITEM     i, end;              /* loop variables for top level */
  REPONODE *node;               /* to traverse the top-level nodes */

  assert(rpt && rep);           /* check the function arguments */
  rpt->min = supp;              /* note the minimum support and */
  rpt->rep = rep;               /* the item set reporter */
  if (rpt->dir < 0) { i = rpt->size-1; end = -1; }
  else              { i = 0; end = rpt->size; }
  if (max < 0) {                /* if maximal item sets only */
    if (isr_xable(rep, 1)) {    /* if item set may be extended */
      for ( ; i != end; i += rpt->dir) {
        node = rpt->tops +i;    /* traverse the top-level nodes */
        if (node->supp < supp) continue;/* skip infrequent items */
        x = -1;                 /* set flag for a frequent child */
        r = isr_addnc(rep, i, node->supp);
        if (r < 0) return r;    /* add current item to the reporter */
        r = maxonly(rpt, node); /* recursively report item sets, */
        isr_remove(rep, 1);     /* then remove the current item */
        if (r < 0) return r;    /* from the item set reporter; */
      } }                       /* finally check for an error */
    else {                      /* if item set may not be extended */
      for (i = 0; i < rpt->size; i++)
        if (rpt->tops[i].supp >= supp) { x = -1; break; }
    } }                         /* check for a frequent item */
  else if (max > 0) {           /* if maximal item sets (ext. filter) */
    if (isr_xable(rep, 1)) {    /* if item set may be extended */
      for ( ; i != end; i += rpt->dir) {
        node = rpt->tops +i;    /* traverse the top-level nodes */
        if (node->supp < supp) continue;/* skip infrequent items */
        x = -1;                 /* set flag for a frequent child */
        r = isr_addnc(rep, i, node->supp);
        if (r < 0) return r;    /* add current item to the reporter */
        r = maximal(rpt, node); /* recursively report item sets, */
        isr_remove(rep, 1);     /* then remove the current item */
        if (r < 0) return r;    /* from the item set reporter */
      } }                       /* finally check for an error */
    else {                      /* if item set may not be extended */
      for (i = 0; i < rpt->size; i++)
        if (rpt->tops[i].supp >= supp) { x = -1; break; }
    } }                         /* check for a frequent item */
  else {                        /* if closed item sets */
    if (isr_xable(rep, 1)) {    /* if item set may be extended */
      for ( ; i != end; i += rpt->dir) {
        node = rpt->tops +i;    /* traverse the top-level nodes */
        if (node->supp < supp) continue;/* skip infrequent items */
        x |= (node->supp >= rpt->supp); /* set perfect ext. flag */
        r  = isr_addnc(rep, i, node->supp);
        if (r < 0) return r;    /* add current item to the reporter */
        r  = closed(rpt, node); /* recursively report item sets, */
        isr_remove(rep, 1);     /* then remove the current item */
        if (r < 0) return r;    /* from the item set reporter, */
      } }                       /* finally check for an error */
    else {                      /* if item set may not be extended */
      for (i = 0; i < rpt->size; i++)
        if (rpt->tops[i].supp >= rpt->supp) { x = -1; break; }
    }                           /* check for a perfect extension */
  }                             /* (renders empty set non-closed) */
  if ((rpt->supp >= supp) && !x)/* if empty set is closed/maximal, */
    return isr_report(rep);     /* report the empty set */
  return 0;                     /* return 'ok' */
}  /* rpt_report() */

/*--------------------------------------------------------------------*/
#ifndef NDEBUG

static void indent (int k)
{ while (--k >= 0) printf("   "); }

/*--------------------------------------------------------------------*/

static void show (REPONODE *node, ITEMBASE *base, int ind)
{                               /* --- recursively show nodes */
  assert(ind >= 0);             /* check the function arguments */
  while (node) {                /* traverse the node list */
    indent(ind);                /* indent the output line */
    if (base) printf("%s/", ib_name(base,node->item));
    printf("%"ITEM_FMT":", node->item);
    printf("%g\n", (double)node->supp);/* print node information */
    show(node->children, base, ind+1);
    node = node->sibling;       /* recursively show child nodes, */
  }                             /* then go to the next node */
}  /* show() */

/*--------------------------------------------------------------------*/

void rpt_show (REPOTREE *rpt, ITEMBASE *base)
{                               /* --- print a repository tree */
  ITEM     i;                   /* loop variable */
  REPONODE *node;               /* to traverse the top-level nodes */

  assert(rpt);                  /* check the function arguments */
  if (!rpt) {                   /* check whether tree exists */
    printf("(null)\n"); return; }
  printf("*:%g\n", (double)rpt->supp);
  for (i = 0; i < rpt->size; i++) {
    node = rpt->tops +i;        /* traverse the top-level nodes */
    if (node->supp <= 0) continue;
    printf("   ");              /* indent the output line */
    if (base) printf("%s/", ib_name(base, i));
    printf("%"ITEM_FMT":", i);  /* print the node information */
    printf("%g\n", (double)node->supp);
    show(node->children, base, 2);
  }                             /* recursively show the nodes */
  printf("nodes: %"SIZE_FMT"\n", ms_used(rpt->mem) +(size_t)rpt->size);
}  /* rpt_show() */

#endif
