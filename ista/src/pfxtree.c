/*----------------------------------------------------------------------
  File    : pfxtree.c
  Contents: prefix tree management for item sets
  Author  : Christian Borgelt
  History : 2009.10.08 file created
            2009.10.09 choice of item order direction added
            2009.10.13 function pxt_isect() added  (intersect trans.)
            2009.10.30 function pxt_add() extended (support increment)
            2009.10.30 function pxt_report() added (recursive reporting)
            2009.10.31 merge functions improved  (for empty arguments)
            2009.10.08 bug in report() concerning set size limit fixed
            2009.11.12 root node made fixed element in tree structure
            2009.11.16 direct combination of intersection with tree
            2009.11.17 extended use of sentinels in isect recursion
            2009.11.18 function pxt_prune() added (support pruning)
            2010.02.04 reduced to functions relevant for IsTa
            2010.06.25 support-based pruning added to isect functions
            2010.08.05 update step counter added to the nodes
            2010.09.07 transaction represented by membership flags
            2010.12.07 added some explicit type casts (for C++)
            2012.04.09 (non-harmful) bug in function maximal() fixed
            2012.04.29 special maximal item set functions added
            2013.04.01 adapted to type changes in module tract
            2013.10.15 checks of return code of isr_report() added
----------------------------------------------------------------------*/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <limits.h>
#include <assert.h>
#include "arrays.h"
#include "pfxtree.h"
#ifdef STORAGE
#include "storage.h"
#endif

/*----------------------------------------------------------------------
  Preprocessor Definitions
----------------------------------------------------------------------*/
#define pos(x,y)  ((x) < (y))   /* macros for item comparison */
#define neg(x,y)  ((x) > (y))   /* (ascending and descending) */

/*----------------------------------------------------------------------
  Main Functions
----------------------------------------------------------------------*/

PFXTREE* pxt_create (ITEM size, int dir, MEMSYS *mem)
{                               /* --- create a prefix tree */
  PFXTREE *pxt;                 /* created prefix tree */

  assert(size >= 0);            /* check the function arguments */
  pxt = (PFXTREE*)malloc(sizeof(PFXTREE)+(size_t)(size-1)*sizeof(SUPP));
  if (!pxt) return NULL;        /* create a prefix tree and */
  pxt->mem  = (mem) ? mem : ms_create(sizeof(PFXNODE), 65535);
  if (!pxt->mem) { free(pxt); return NULL; }
  pxt->size = size;             /* create memory system for nodes */
  pxt->dir  = (dir < 0) ? -1 : +1;
  pxt->step = 0;                /* note item count and order direction */
  pxt->last = 0;                /* and initialize the other fields */
  pxt->min  = pxt->supp = 0;
  pxt->rep  = NULL;
  pxt->root.item     = -1;      /* initialize the root node */
  pxt->root.supp     = 0; pxt->root.step = 0;
  pxt->root.children = pxt->root.sibling = NULL;
  return pxt;                   /* return the created prefix tree */
}  /* pxt_create() */

/*--------------------------------------------------------------------*/

static void delete (PFXNODE *node, MEMSYS *mem)
{                               /* --- recursively delete nodes */
  PFXNODE *tmp;                 /* buffer for deallocation */

  assert(mem);                  /* check the function arguments */
  while (node) {                /* sibling list deletion loop */
    delete(node->children,mem); /* recursively delete the children */
    tmp = node; node = node->sibling; ms_free(mem, tmp);
  }                             /* finally delete the node itself */
}  /* delete() */

/*--------------------------------------------------------------------*/

void pxt_delete (PFXTREE *pxt, int delms)
{                               /* --- delete a prefix tree */
  assert(pxt);                  /* check the function arguments */
  if (delms)                    /* if possible/requested, */
    ms_delete(pxt->mem);        /* delete the memory system */
  else                          /* otherwise delete the nodes */
    delete(pxt->root.children, pxt->mem);
  free(pxt);                    /* delete the base structure */
}  /* pxt_delete() */

/*--------------------------------------------------------------------*/

int pxt_add (PFXTREE *pxt, const ITEM *items, ITEM n, SUPP supp)
{                               /* --- add item set to prefix tree */
  ITEM    i;                    /* buffer for an item */
  PFXNODE **p;                  /* pointer to insertion position */
  PFXNODE *node;                /* to insert new nodes */

  assert(pxt                    /* check the function arguments */
  &&    (items || (n <= 0)) && (supp >= 0));
  node = &pxt->root;            /* start at the root node */
  do {                          /* traverse the items of the set */
    if (supp > node->supp)      /* adapt the node support */
      node->supp = supp;        /* (root represents empty set) */
    if (--n < 0) return 0;      /* if all items are processed, abort */
    i = *items++;               /* get the next item in the set */
    assert((i >= 0) && (i < pxt->size));
    p = &node->children;        /* traverse the list of children */
    if (pxt->dir < 0) while (*p && ((*p)->item > i)) p = &(*p)->sibling;
    else              while (*p && ((*p)->item < i)) p = &(*p)->sibling;
    node = *p;                  /* find the item/insertion position */
  } while (node && (node->item == i));
  node = (PFXNODE*)ms_alloc(pxt->mem);
  if (!node) return -1;         /* create a new prefix tree node */
  node->item    = i;            /* store the current item, */
  node->step    = 0;            /* clear the update step, and */
  node->supp    = supp;         /* set support of the item set */
  node->sibling = *p;           /* insert the created node */
  *p = node;                    /* into the sibling list */
  while (--n >= 0) {            /* traverse the rest of the items */
    node = node->children = (PFXNODE*)ms_alloc(pxt->mem);
    if (!node) return -1;       /* create a new prefix tree node */
    assert((*items >= 0) && (*items < pxt->size));
    node->item    = *items++;   /* store the current item, */
    node->step    = 0;          /* clear the update step, and */
    node->supp    = supp;       /* set support of the item set */
    node->sibling = NULL;       /* there are no siblings yet */
  }
  node->children = NULL;        /* last created node is a leaf */
  return 0;                     /* return 'ok' */
}  /* pxt_add() */

/*--------------------------------------------------------------------*/

#define ISECT(dir) \
static int isect_##dir (PFXNODE *node, PFXNODE **ins, PFXTREE *pxt)    \
{                               /* --- intersect with an item set */   \
  ITEM    i;                    /* buffer for current item */          \
  PFXNODE *d;                   /* to allocate new nodes */            \
                                                                       \
  assert(node && ins && pxt);   /* check the function arguments */     \
  for ( ; node; node = node->sibling) {                                \
    i = node->item;             /* traverse the node list */           \
    if (node->step >= pxt->step) { /* if the node has been visited */  \
      if (!dir(i, pxt->last)) break;                                   \
      if (node->children        /* if the node has child nodes */      \
      && (isect_##dir(node->children, &node->children, pxt) < 0))      \
        return -1; }            /* intersect subtree with item set */  \
    else if (!pxt->mins[i]) {   /* if item is not in intersection */   \
      if (!dir(i, pxt->last)) break;                                   \
      if (node->children        /* if there are child nodes */         \
      && (isect_##dir(node->children, ins,             pxt) < 0))      \
        return -1; }            /* intersect subtree with item set */  \
    else if (node->supp < pxt->mins[i]) { /* if not enough support */  \
      if (!dir(i, pxt->last)) break; }    /* skip the item */          \
    else {                      /* if item is in the intersection */   \
      while ((d = *ins) && dir(d->item, i))                            \
        ins = &d->sibling;      /* find the insertion position */      \
      if (!d || (d->item != i)){/* if node does not exist */           \
        d = (PFXNODE*)ms_alloc(pxt->mem);                              \
        if (!d) return -1;      /* allocate a new node and */          \
        d->item = i;            /* store the matched item */           \
        d->step = pxt->step;    /* set update step and support */      \
        d->supp = pxt->supp +node->supp;                               \
        d->sibling  = *ins; *ins = d;                                  \
        d->children = NULL; }   /* insert node into sibling list */    \
      else {                    /* if a node already exists */         \
        if (d->step >= pxt->step) d->supp -= pxt->supp;                \
        if (d->supp < node->supp) d->supp = node->supp;                \
        d->supp += pxt->supp;   /* update the intersection support */  \
        d->step  = pxt->step;   /* and set the current update step */  \
      }                         /* (update current intersection) */    \
      if (!dir(i, pxt->last)) break;                                   \
      if (node->children        /* if there are child nodes */         \
      && (isect_##dir(node->children, &d->children, pxt) < 0))         \
        return -1;              /* recursively intersect subtree */    \
    }                           /* with the rest of the item set */    \
  }                                                                    \
  return 0;                     /* return 'ok' */                      \
}  /* isect() */

/*--------------------------------------------------------------------*/

ISECT(pos)                      /* function for ascending  item order */
ISECT(neg)                      /* function for descending item order */

/*--------------------------------------------------------------------*/

int pxt_isect (PFXTREE *pxt, const ITEM *items, ITEM n, SUPP supp,
               SUPP min, const SUPP *frqs)
{                               /* --- intersect with an item set */
  ITEM i;                       /* to traverse the items */
  SUPP s;                       /* to compute limiting support */

  assert(pxt && (items || (n <= 0)));  /* check function arguments */
  pxt->root.supp += supp;       /* update the empty set support */
  if (n <= 0) return 0;         /* handle empty item set specially */
  if (pxt_add(pxt, items, n, 0) < 0)
    return -1;                  /* add the item set to the tree */
  pxt->last = items[n-1];       /* note last item in transaction */
  pxt->supp = supp;             /* and the transaction support */
  memset(pxt->mins, 0, (size_t)pxt->size *sizeof(SUPP));
  if (!frqs) min = 0;           /* clear all item flags/supports */
  for (s = 0; --n >= 0; ) {     /* traverse items in transaction */
    i = items[n];               /* get the next item */
    if (frqs && (frqs[i] > s)) s = frqs[i];
    pxt->mins[i] = (min > s) ? min-s : (SUPP)-1;
  }                             /* compute minimum support values */
  /* For the pruning, it is not sufficient to use only the frequency  */
  /* of an item in the future transactions. Rather the frequency of   */
  /* all items with codes following the code of the item have to be   */
  /* taken into account (find maximum frequency). Otherwise items in  */
  /* subtrees rooted at nodes with the item, which could still become */
  /* become frequent, may not be processed, thus losing results.      */
  pxt->step++;                  /* increment the update step */
  return (pxt->dir < 0)         /* intersect tree with item set */
    ? isect_neg(pxt->root.children, &pxt->root.children, pxt)
    : isect_pos(pxt->root.children, &pxt->root.children, pxt);
}  /* pxt_isect() */

/*--------------------------------------------------------------------*/

SUPP pxt_get (PFXTREE *pxt, const ITEM *items, ITEM n)
{                               /* --- get support of an item set */
  ITEM    i;                    /* buffer for an item */
  PFXNODE *p;                   /* to traverse the nodes */

  assert(pxt && (items || (n <= 0))); /* check function arguments */
  p = &pxt->root;               /* start at the root node */
  while (--n >= 0) {            /* while not at the last item */
    p = p->children;            /* continue with the child nodes */
    i = *items++;               /* try to find a corresp. child node */
    if (pxt->dir < 0) while (p && (p->item > i)) p = p->sibling;
    else              while (p && (p->item < i)) p = p->sibling;
    if (!p || (p->item != i))   /* if a node with the next item */
      return -1;                /* does not exist in the tree, */
  }                             /* abort the search with failure */
  return p->supp;               /* return support of the item set */
}  /* pxt_get() */

/*--------------------------------------------------------------------*/
#if 1                           /* check siblings before children */

#define SUPER(dir) \
static int super_##dir (PFXNODE *node,                                 \
                        const ITEM *items, ITEM n, SUPP supp)          \
{                               /* --- check for a superset */         \
  assert(items && (n > 0) && (supp > 0));   /* check arguments */      \
  while (node                   /* while there is another node with */ \
  && !dir(*items, node->item)){ /* an item not after next to match */  \
    if (*items == node->item) { /* if at node with matching item */    \
      if (--n <= 0) return (node->supp >= supp);                       \
      items++; }                /* check for last and skip item */     \
    else if (super_##dir(node->sibling, items, n, supp))               \
      return -1;                /* process siblings of current node */ \
    if (node->supp < supp)      /* if the support is too low, */       \
      return 0;                 /* no superset will be found */        \
    node = node->children;      /* continue with the child nodes */    \
  }                             /* (match the remaining items) */      \
  return 0;                     /* return 'no superset exists' */      \
}  /* super() */

/*--------------------------------------------------------------------*/
#else                           /* check children before siblings */

#define SUPER(dir) \
static int super_##dir (PFXNODE *node,                                 \
                        const ITEM *items, ITEM n, SUPP supp)          \
{                               /* --- check for a superset */         \
  assert(items && (n > 0) && (supp > 0));   /* check arguments */      \
  for ( ; node && !dir(*items, node->item); node = node->sibling) {    \
    if (node->supp < supp)      /* traverse the sibling list, but */   \
      continue;                 /* skip nodes with insuff. support */  \
    if (node->item == *items) { /* if at node with matching item */    \
      if (n <= 1) return -1;    /* if last item was matched, abort */  \
      if (super_##dir(node->children, items+1, n-1, supp))             \
        return -1; }            /* check children for rem. items */    \
    else {                      /* if at node with different item */   \
      if (super_##dir(node->children, items,   n,   supp))             \
        return -1;              /* check children for rem. items */    \
    }                           /* and abort if a superset is found */ \
  }                                                                    \
  return 0;                     /* return 'no superset exists' */      \
}  /* super() */

#endif
/*--------------------------------------------------------------------*/

SUPER(pos)                      /* function for ascending  item order */
SUPER(neg)                      /* function for descending item order */

/*--------------------------------------------------------------------*/

int pxt_super (PFXTREE *pxt, const ITEM *items, ITEM n, SUPP supp)
{                               /* --- check for a superset */
  assert(pxt                    /* check the function arguments */
  &&    (items || (n <= 0)) && (supp > 0));
  if (n <= 0)                   /* if no item, check root support */
    return (pxt->root.supp >= supp) ? -1 : 0;
  return (pxt->dir < 0)         /* recursively check for a superset */
       ? super_neg(pxt->root.children, items, n, supp)
       : super_pos(pxt->root.children, items, n, supp);
}  /* pxt_super() */

/*--------------------------------------------------------------------*/

#define MERGE(dir) \
static PFXNODE* merge_##dir (PFXNODE *s1, PFXNODE *s2, MEMSYS *mem)    \
{                               /* --- merge two node list */          \
  PFXNODE *out, **end, *node;   /* output node list and end pointer */ \
                                                                       \
  assert(mem);                  /* check the function arguments */     \
  if (!s1) return s2;           /* if there is only one node list, */  \
  if (!s2) return s1;           /* simply return the other list */     \
  end = &out;                   /* start the output list */            \
  while (1) {                   /* node list merge loop */             \
    if      (dir(s1->item, s2->item)) {                                \
      *end = s1; end = &s1->sibling; if (!(s1 = *end)) break; }        \
    else if (dir(s2->item, s1->item)) {                                \
      *end = s2; end = &s2->sibling; if (!(s2 = *end)) break; }        \
    else {                      /* copy nodes with singular items */   \
      if (s1->supp < s2->supp)  /* same item: keep the maximum */      \
        s1->supp = s2->supp;    /* of the support values */            \
      s1->children = merge_##dir(s1->children, s2->children, mem);     \
      node = s2; s2  =  s2->sibling; ms_free(mem, node);               \
      *end = s1; end = &s1->sibling; s1 = *end;                        \
      if (!s1 || !s2) break;    /* if an item occurs in both lists, */ \
    }                           /* rec. merge the lists of children */ \
  }                             /* and delete one of the two nodes */  \
  *end = (s1) ? s1 : s2;        /* append the remaining nodes */       \
  return out;                   /* return the merged prefix tree */    \
}  /* merge() */

/*--------------------------------------------------------------------*/

MERGE(pos)                      /* function for ascending  item order */
MERGE(neg)                      /* function for descending item order */

/*--------------------------------------------------------------------*/

#define PRUNEX(dir) \
static void prunex_##dir (PFXNODE *node, PFXTREE *pxt)                 \
{                               /* --- prune a prefix tree */          \
  PFXNODE *n, *t;               /* to traverse the node list */        \
  PFXNODE *o, **p = &o;         /* nodes to keep and end pointer */    \
                                                                       \
  assert(node && pxt);          /* check the function arguments */     \
  n = node->children; node->children = NULL;                           \
  while (n) {                   /* traverse the sibling list */        \
    if (n->children)            /* recursively prune the children */   \
      prunex_##dir(n, pxt);     /* before pruning the node itself */   \
    if (n->supp >= pxt->mins[n->item]) {                               \
      *p = n; p = &n->sibling;  /* if item (set) may be frequent, */   \
      n = *p; }                 /* keep the corresponding node */      \
    else {                      /* if item (set) is infrequent */      \
      node->children =          /* (and cannot become frequent) */     \
        merge_##dir(node->children, n->children, pxt->mem);            \
      t = n; n = n->sibling;    /* merge the child nodes */            \
      ms_free(pxt->mem, t);     /* with the pruned subtrees */         \
    }                           /* and delete the processed node */    \
  }                             /* finally terminate the keep list */  \
  *p = NULL;                    /* and merge the two output lists */   \
  node->children = merge_##dir(node->children, o, pxt->mem);           \
}  /* prunex() */

/*--------------------------------------------------------------------*/

PRUNEX(pos)                     /* function for ascending  item order */
PRUNEX(neg)                     /* function for descending item order */

/*--------------------------------------------------------------------*/

int pxt_prunex (PFXTREE *pxt, SUPP supp, const SUPP *frqs)
{                               /* --- prune infrequent item sets */
  ITEM i;                       /* loop variable */

  assert(pxt && (supp > 0) && frqs);  /* check the function arguments */
  for (i = 0; i < pxt->size;i++)/* traverse the items and */
    pxt->mins[i] = supp-frqs[i];/* compute minimum support values */
  if (pxt->dir < 0) prunex_neg(&pxt->root, pxt);
  else              prunex_pos(&pxt->root, pxt);
  return 0;                     /* prune prefix tree and return 'ok' */
}  /* pxt_prunex() */

/*--------------------------------------------------------------------*/

static void prune (PFXNODE **node, SUPP supp, MEMSYS *mem)
{                               /* --- recursively prune tree */
  PFXNODE *t;                   /* temporary buffer for deletion */

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

void pxt_prune (PFXTREE *pxt, SUPP supp)
{                               /* --- prune item set repository */
  assert(pxt && (supp >= 0));   /* check the function arguments */
  prune(&pxt->root.children, supp, pxt->mem);
}  /* pxt_prune() */            /* recursively prune prefix tree */

/*--------------------------------------------------------------------*/

static int closed (PFXTREE *pxt, PFXNODE *node)
{                               /* --- report closed item sets */
  int  r, x = 0;                /* error status, perfect ext. flag */
  SUPP supp;                    /* support of current item set */

  assert(pxt && node);          /* check the function arguments */
  supp = node->supp;            /* get current item set support */
  if (isr_xable(pxt->rep, 1)) { /* if item set may be extended */
    for (node = node->children; node; node = node->sibling) {
      if (node->supp < pxt->min)/* traverse the node nodes, */
        continue;               /* but skip infrequent item sets */
      x |= (node->supp >= supp);/* set perfect extension flag */
      r  = isr_addnc(pxt->rep, node->item, node->supp);
      if (r < 0) return r;      /* add current item to the reporter */
      r  = closed(pxt, node);   /* recursively report item sets, */
      isr_remove(pxt->rep, 1);  /* then remove the current item */
      if (r < 0) return r;      /* from the item set reporter, */
    } }                         /* finally check for an error */
  else {                        /* if item set may not be extended */
    for (node = node->children; node; node = node->sibling)
      if (node->supp >= supp) { x = -1; break; }
  }                             /* check for a perfect extension */
  return (x) ? 0 : isr_report(pxt->rep);
}  /* closed() */               /* report sets without perfect exts. */

/*--------------------------------------------------------------------*/

static int maximal (PFXTREE *pxt, PFXNODE *node)
{                               /* --- report maximal item sets */
  int r, x = 0;                 /* error status, freq. child flag */

  assert(pxt && node);          /* check the function arguments */
  if (isr_xable(pxt->rep, 1)) { /* if item set may be extended */
    for (node = node->children; node; node = node->sibling) {
      if (node->supp < pxt->min)/* traverse the child nodes, */
        continue;               /* but skip infrequent item sets */
      x = -1;                   /* set flag for a frequent child */
      r = isr_addnc(pxt->rep, node->item, node->supp);
      if (r < 0) return r;      /* add current item to the reporter */
      r = maximal(pxt, node);   /* recursively report item sets */
      isr_remove(pxt->rep, 1);  /* and remove the current item */
      if (r < 0) return r;      /* from the item set reporter */
    } }                         /* finally check for an error */
  else {                        /* if item set may not be extended */
    for (node = node->children; node; node = node->sibling)
      if (node->supp >= pxt->min) { x = -1; break; }
  }                             /* check for a frequent superset */
  return (x) ? 0 : isr_report(pxt->rep);
}  /* maximal() */              /* report sets w/o frequent supersets */

/*--------------------------------------------------------------------*/

static int maxonly (PFXTREE *pxt, PFXNODE *node)
{                               /* --- report maximal item sets */
  int        r, x = 0;          /* error status, freq. child flag */
  ITEM       n;                 /* number of items in set */
  const ITEM *s;                /* current item set */
  PFXNODE    *curr = node;      /* buffer for current node */

  assert(pxt && node);          /* check the function arguments */
  if (isr_xable(pxt->rep, 1)) { /* if item set may be extended */
    for (node = node->children; node; node = node->sibling) {
      if (node->supp < pxt->min)/* traverse the child nodes, */
        continue;               /* but skip infrequent item sets */
      x = -1;                   /* set flag for a frequent child */
      r = isr_addnc(pxt->rep, node->item, node->supp);
      if (r < 0) return r;      /* add current item to the reporter */
      r = maxonly(pxt, node);   /* recursively report item sets */
      isr_remove(pxt->rep, 1);  /* and remove the current item */
      if (r < 0) return r;      /* from the item set reporter */
    } }                         /* finally check for an error */
  else {                        /* if item set may not be extended */
    for (node = node->children; node; node = node->sibling)
      if (node->supp >= pxt->min) { x = -1; break; }
  }                             /* check for a frequent superset */
  if (x) return 0;              /* if there is a freq. child, abort */
  curr->supp = -curr->supp;     /* mark the current node as to ignore */
  n = isr_cnt  (pxt->rep);      /* get the current number of items */
  s = isr_items(pxt->rep);      /* and the current item set and */
  x = (pxt->dir < 0)            /* recursively search for a superset */
    ? super_neg(pxt->root.children, s, n, pxt->min)
    : super_pos(pxt->root.children, s, n, pxt->min);
  curr->supp = -curr->supp;     /* unmark the current node */
  return (x) ? 0 : isr_report(pxt->rep);
}  /* maxonly() */              /* report sets w/o frequent supersets */

/*--------------------------------------------------------------------*/

int pxt_report (PFXTREE *pxt, int max, SUPP supp, ISREPORT *rep)
{                               /* --- report (closed) item sets */
  assert(pxt && rep);           /* check the function arguments */
  pxt->min = supp;              /* note the minimum support and */
  pxt->rep = rep;               /* the item set reporter */
  if      (max < 0)             /* if maximal item sets only */
    return maxonly(pxt, &pxt->root);
  else if (max > 0)             /* if maximal item sets (ext. filter) */
    return maximal(pxt, &pxt->root);
  else                          /* if closed  item sets */
    return closed (pxt, &pxt->root);
}  /* pxt_report() */

/*--------------------------------------------------------------------*/
#ifndef NDEBUG

static void indent (int k)
{ while (--k >= 0) printf("   "); }

/*--------------------------------------------------------------------*/

static void show (PFXNODE *node, ITEMBASE *base, int ind)
{                               /* --- recursively show nodes */
  assert(ind >= 0);             /* check the function arguments */
  while (node) {                /* traverse the node list */
    indent(ind);                /* indent the output line */
    if (base) printf("%s/", ib_name(base, node->item));
    printf("%"ITEM_FMT":",  node->item); /* print the node */
    printf("%"SUPP_FMT"\n", node->supp); /* information */
    show(node->children, base, ind+1);
    node = node->sibling;       /* recursively show the child nodes, */
  }                             /* then go to the next node */
}  /* show() */

/*--------------------------------------------------------------------*/

void pxt_show (PFXTREE *pxt, ITEMBASE *base)
{                               /* --- print a prefix tree */
  assert(pxt);                  /* check the function arguments */
  if (!pxt) {                   /* check whether tree exists */
    printf("(null)\n"); return; }
  show(pxt->root.children, base, 0); /* recursively show the nodes */
  printf("supp:  %"SUPP_FMT"\n", pxt->root.supp);
  printf("nodes: %"SIZE_FMT"\n", ms_used(pxt->mem));
}  /* pxt_show() */             /* print global information */

#endif
