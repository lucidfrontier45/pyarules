/*----------------------------------------------------------------------
  File    : pattree.c
  Contents: patricia tree management for item sets
  Author  : Christian Borgelt
  History : 2012.07.06 file created from pfxtree.c
            2012.07.09 first version completed (without pruning)
            2012.07.13 pruning added (merge(), prunex(), pat_prunex())
            2012.07.16 pruning improved (full reduction of items)
            2013.04.01 adapted to type changes in module tract
            2013.10.15 checks of return code of isr_report() added
----------------------------------------------------------------------*/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <limits.h>
#include <assert.h>
#include "arrays.h"
#include "pattree.h"
#ifdef STORAGE
#include "storage.h"
#endif

/*----------------------------------------------------------------------
  Preprocessor Definitions
----------------------------------------------------------------------*/
#define pos(x,y)  ((x) < (y))   /* macros for item comparison */
#define neg(x,y)  ((x) > (y))   /* (ascending and descending) */

#define PATNODESIZE(n)  (sizeof(PATNODE) +(size_t)((n)-1)*sizeof(int))

/*----------------------------------------------------------------------
  Main Functions
----------------------------------------------------------------------*/

PATTREE* pat_create (ITEM size, int dir)
{                               /* --- create a patricia tree */
  PATTREE *pat;                 /* created patricia tree */

  assert(size >= 0);            /* check the function arguments */
  pat = (PATTREE*)malloc(sizeof(PATTREE)+(size_t)(size-1)*sizeof(SUPP)
                                        +(size_t) size   *sizeof(ITEM));
  if (!pat) return NULL;        /* create a patricia tree and */
  pat->size  = size;            /* initialize its fields */
  pat->cnt   = pat->max = 0;
  pat->dir   = (dir < 0) ? -1 : +1;
  pat->step  = 0; pat->last = 0;
  pat->min   = pat->supp = 0;
  pat->items = (ITEM*)(pat->mins +size);
  pat->root.supp     = 0; pat->root.step = 0;
  pat->root.children = pat->root.sibling = NULL;
  pat->root.cnt      = 0;       /* initialize the root node */
  pat->root.items[0] = ITEM_MAX;/* (empty patricia tree) */
  return pat;                   /* return the created patricia tree */
}  /* pat_create() */

/*--------------------------------------------------------------------*/

static void delete (PATNODE *node)
{                               /* --- recursively delete nodes */
  PATNODE *tmp;                 /* buffer for deallocation */

  while (node) {                /* sibling list deletion loop */
    delete(node->children);     /* recursively delete the children */
    tmp = node; node = node->sibling; free(tmp);
  }                             /* finally delete the node itself */
}  /* delete() */

/*--------------------------------------------------------------------*/

void pat_delete (PATTREE *pat)
{                               /* --- delete a patricia tree */
  assert(pat);                  /* check the function argument */
  delete(pat->root.children);   /* delete the tree nodes */
  free(pat);                    /* delete the base structure */
}  /* pat_delete() */

/*--------------------------------------------------------------------*/

static PATNODE* split (PATNODE *node, ITEM n)
{                               /* --- split a node at the n-th item */
  ITEM    k;                    /* number of items in child node */
  PATNODE *child;               /* created child node */

  assert(node                   /* check the function arguments */
  &&    (n > 0) && (n < node->cnt));
  k     = node->cnt -n;         /* get number of items in child */
  child = (PATNODE*)malloc(PATNODESIZE(k));
  if (!child) return NULL;      /* create a child node and */
  child->step     = node->step; /* copy and init. the fields */
  child->supp     = node->supp;
  child->sibling  = NULL;
  child->children = node->children;
  child->cnt      = k;          /* copy the last k items */
  memcpy(child->items, node->items +n, (size_t)k *sizeof(ITEM));
  node->children  = child;
  node->cnt       = n;          /* shrink the original node */
  return (PATNODE*)realloc(node, PATNODESIZE(n));
}  /* split() */

/*--------------------------------------------------------------------*/

static PATNODE* expand (PATNODE *node, const ITEM *items, ITEM n)
{                               /* --- expand a node with add. items */
  assert(node && items && (n > 0));  /* check the function arguments */
  node = (PATNODE*)realloc(node, PATNODESIZE(node->cnt +n));
  if (!node) return NULL;       /* resize the node */
  memcpy(node->items +node->cnt, items, (size_t)n *sizeof(ITEM));
  node->cnt += n;               /* copy the additional items */
  return node;                  /* return the resized node */
}  /* expand() */

/*--------------------------------------------------------------------*/

int pat_add (PATTREE *pat, const ITEM *items, ITEM n, SUPP supp)
{                               /* --- add item set to patricia tree */
  ITEM    i;                    /* item buffer, */
  PATNODE **p;                  /* pointer to insertion position */
  PATNODE *node;                /* to insert new nodes */

  assert(pat                    /* check the function arguments */
  &&    (items || (n <= 0)) && (supp >= 0));
  node = &pat->root;            /* start at the root node */
  while (1) {                   /* traverse the items of the set */
    if (supp > node->supp)      /* adapt the node/item set support */
      node->supp = supp;        /* (note: root represents empty set) */
    if (n <= 0) return 0;       /* if all items are processed, abort */
    i = *items++;               /* get the next item in the set */
    p = &node->children;        /* traverse the list of children and */
    if (pat->dir < 0)           /* find the item/insertion position */
         while (*p && ((*p)->items[0] > i)) p = &(*p)->sibling;
    else while (*p && ((*p)->items[0] < i)) p = &(*p)->sibling;
    node = *p;                  /* if no node with the next item */
    if (!node || (node->items[0] != i))
      break;                    /* abort the search loop */
    for (i = 1; (--n > 0) && (i < node->cnt); items++, i++)
      if (node->items[i] != *items)
        break;                  /* compare other items in node */
    if (i < node->cnt) {        /* if not all items were matched */
      if ((n <= 0) && (node->supp >= supp))
        return 0;               /* if items are covered, abort */
      node = split(node, i);    /* split the current node */
      if (!node) return -1;     /* at the first unmatched item */
      *p = node;                /* and count the added node */
      if (++pat->cnt > pat->max) pat->max = pat->cnt; }
    else if (!node->children    /* if the node may be extended */
    &&      (n > 0) && (supp >= node->supp)) {
      node->supp = supp;        /* update the node support */
      node = expand(node, items, n);
      if (!node) return -1;     /* expand the tree node */
      *p = node; return  0;     /* with the remaining items */
    }                           /* and update the node pointer */
  }                             /* (may be changed due to realloc) */
  node = (PATNODE*)malloc(PATNODESIZE(n));
  if (!node) return -1;         /* create a new patricia tree node */
  node->step     = 0;           /* clear the update step and */
  node->supp     = supp;        /* store the support of the item set */
  node->sibling  = *p; *p = node;
  node->children = NULL;        /* insert node into sibling list */
  node->cnt      = n;           /* copy the unprocessed items */
  memcpy(node->items, items-1, (size_t)n *sizeof(ITEM));
  if (++pat->cnt > pat->max) pat->max = pat->cnt;
  return 0;                     /* count added node and return 'ok' */
}  /* pat_add() */

/*--------------------------------------------------------------------*/

#define INSERT(dir) \
static PATNODE** insert_##dir (PATNODE **ins, const ITEM *items,       \
                               ITEM n, SUPP supp, PATTREE *pat)        \
{                               /* --- insert for intersection */      \
  ITEM    i;                    /* item buffer, loop variable */       \
  PATNODE *node;                /* to traverse/insert nodes */         \
                                                                       \
  assert(ins                    /* check the function arguments */     \
  &&     items && (n > 0) && (supp >= 0) && pat);                      \
  while (1) {                   /* item insertion loop */              \
    i = *items++;               /* get the first item to insert */     \
    while ((node = *ins) && dir(node->items[0], i))                    \
      ins = &node->sibling;     /* find the insertion position */      \
    if (!node || (node->items[0] != i)) { /* if no node exists */      \
      node = (PATNODE*)malloc(PATNODESIZE(n));                         \
      if (!node) return NULL;   /* create a new patricia tree node */  \
      node->step = pat->step;   /* store support of the item set */    \
      node->supp = pat->supp +supp;                                    \
      node->sibling  = *ins; *ins = node;                              \
      node->children = NULL;    /* insert node into sibling list */    \
      node->cnt      = n;       /* copy the items to insert */         \
      memcpy(node->items, items-1, (size_t)n *sizeof(ITEM));           \
      pat->cnt++;               /* count the new node (just added) */  \
      return &node->children;   /* return new insertion position */    \
    }                                                                  \
    for (i = 1; (--n > 0) && (i < node->cnt); items++, i++)            \
      if (node->items[i] != *items)                                    \
        break;                  /* compare other items in the node */  \
    if (i < node->cnt) {        /* if not all items were matched, */   \
      node = split(node, i);    /* split the current node */           \
      if (!node) return NULL;   /* at the first unmatched item */      \
      *ins = node;              /* store the split node (top part) */  \
      pat->cnt++;               /* (may be changed due to realloc) */  \
    }                           /* and count the added node */         \
    if (node->step >= pat->step) node->supp -= pat->supp;              \
    if (node->supp <  supp)      node->supp  = supp;                   \
    node->supp += pat->supp;    /* update the intersection support */  \
    node->step  = pat->step;    /* and set the current update step */  \
    if (n <= 0)                 /* if all items are processed, */      \
      return &node->children;   /* return new insertion position */    \
    if (!node->children         /* if the node may be extended */      \
    && (pat->supp +supp >= node->supp)) {                              \
      node = expand(node, items, n);     /* expand the node */         \
      if (!node)   return NULL; /* with the remaining items */         \
      *ins = node; return &node->children;                             \
    }                           /* return new insertion position */    \
    ins = &node->children;      /* continue with the child */          \
  }                             /* (insert remaining items) */         \
}  /* insert() */

/* Unlike pat_add(), the insertion always has to split a node that  */
/* contains only part of the items to insert, because the insertion */
/* position (children) for future insertions needs to be returned.  */

/*--------------------------------------------------------------------*/

INSERT(pos)                     /* function for ascending  item order */
INSERT(neg)                     /* function for descending item order */

/*--------------------------------------------------------------------*/

#define ISECT(dir) \
static int isect_##dir (PATNODE *node, PATNODE **ins, PATTREE *pat)    \
{                               /* --- intersect with an item set */   \
  ITEM    i, k, n;              /* item buffer, loop variables */      \
  PATNODE *c, *p, **x;          /* to allocate new nodes */            \
                                                                       \
  assert(node && ins && pat);   /* check the function arguments */     \
  for ( ; node; node = node->sibling) {                                \
    if (node->step >= pat->step) { /* if the node has been visited */  \
      if (!dir(node->items[0], pat->last))                             \
        break;                  /* if last item processed, abort */    \
      if (node->children        /* if the node has child nodes */      \
      && (isect_##dir(node->children, &node->children, pat) < 0))      \
        return -1; }            /* intersect subtree with item set */  \
    else {                      /* if node has not been visited */     \
      i = node->items[0];       /* get the first item in the node */   \
      while ((p = *ins) && dir(p->items[0], i))                        \
        ins = &p->sibling;      /* advance the insertion position */   \
      for (k = n = 0; k < node->cnt; k++)                              \
        if (pat->mins[i = node->items[k]])                             \
          pat->items[n++] = i;  /* collect the intersection items */   \
      if (n <= 0) {             /* if the intersection is empty */     \
        if (!dir(node->items[0], pat->last))                           \
          break;                /* if last item processed, abort */    \
        if (node->children      /* if the node has child nodes */      \
        && (isect_##dir(node->children, ins, pat) < 0))                \
          return -1; }          /* intersect subtree with item set */  \
      else if (node->supp < pat->mins[pat->items[0]]) {                \
        if (!dir(node->items[0], pat->last))                           \
          break; }              /* if not enough support, skip item */ \
      else {                    /* if intersection is not empty */     \
        c = node->children;     /* note the child nodes */             \
        /* The child nodes have to be retrieved here, because the  */  \
        /* insertion function may split the current node, changing */  \
        /* the child pointer. This can lead to double processing.  */  \
        x = insert_##dir(ins, pat->items, n, node->supp, pat);         \
        if (!x) return -1;      /* insert the intersection */          \
        if (node == p)          /* if inserting may have updated */    \
          node = *ins;          /* the current node, get it again */   \
        if (!dir(node->items[0], pat->last))                           \
          break;                /* if last item processed, abort */    \
        if (c && (isect_##dir(c, x, pat) < 0))                         \
          return -1;            /* if there are child nodes */         \
      }                         /* recursively intersect subtree */    \
    }                           /* with the rest of the item set */    \
  }                                                                    \
  return 0;                     /* return 'ok' */                      \
}  /* isect() */

/*--------------------------------------------------------------------*/

ISECT(pos)                      /* function for ascending  item order */
ISECT(neg)                      /* function for descending item order */

/*--------------------------------------------------------------------*/

int pat_isect (PATTREE *pat, const ITEM *items, ITEM n, SUPP supp,
               SUPP min, const SUPP *frqs)
{                               /* --- intersect with an item set */
  int  r;                       /* result of isect_dir() */
  ITEM i;                       /* to traverse the items */
  SUPP s;                       /* to compute limiting support */

  assert(pat && (items || (n <= 0)));  /* check function arguments */
  pat->root.supp += supp;       /* update the empty set support */
  if (n <= 0) return 0;         /* handle empty item set specially */
  if (pat_add(pat, items, n, 0) < 0)
    return -1;                  /* add the item set to the tree */
  pat->last = items[n-1];       /* note the last item in the set */
  pat->supp = supp;             /* and the transaction support */
  memset(pat->mins, 0, (size_t)pat->size *sizeof(ITEM));
  if (!frqs) min = 0;           /* clear all item flags/supports */
  for (s = 0; --n >= 0; ) {     /* traverse items in transaction */
    i = items[n];               /* get the next item */
    if (frqs && (frqs[i] > s)) s = frqs[i];
    pat->mins[i] = (min > s) ? min-s : (SUPP)-1;
  }                             /* compute minimum support values */
  /* For the pruning, it is not sufficient to use only the frequency  */
  /* of an item in the future transactions. Rather the frequency of   */
  /* all items with codes following the code of the item have to be   */
  /* taken into account (find maximum frequency). Otherwise items in  */
  /* subtrees rooted at nodes with the item, which could still become */
  /* frequent, may not be processed, and then results may be lost.    */
  pat->step++;                  /* increment the update step */
  r = (pat->dir < 0)            /* intersect transaction with tree */
    ? isect_neg(pat->root.children, &pat->root.children, pat)
    : isect_pos(pat->root.children, &pat->root.children, pat);
  if (pat->cnt > pat->max) pat->max = pat->cnt;
  return r;                     /* return the error status */
}  /* pat_isect() */

/*--------------------------------------------------------------------*/

SUPP pat_get (PATTREE *pat, const ITEM *items, ITEM n)
{                               /* --- get support of an item set */
  ITEM    i;                    /* item buffer, loop variable */
  PATNODE *p;                   /* to traverse the nodes */

  assert(pat && (items || (n <= 0))); /* check function arguments */
  p = &pat->root;               /* start at the root node */
  while (--n >= 0) {            /* while not at the last item */
    p = p->children;            /* continue with the child nodes */
    i = *items++;               /* try to find a corresp. child node */
    if (pat->dir < 0) while (p && (p->items[0] > i)) p = p->sibling;
    else              while (p && (p->items[0] < i)) p = p->sibling;
    if (!p || (p->items[0] != i))  /* if no node with the next item, */
      return -1;                   /* abort the search with failure */
    for (i = 1; i < p->cnt; i++) { /* compare other items in node */
      if (--n < 0)                 return p->supp;
      if (p->items[i] != *items++) return -1;
    }                           /* if an item in the node differs, */
  }                             /* abort the search with failure */
  return p->supp;               /* return support of the item set */
}  /* pat_get() */

/*--------------------------------------------------------------------*/

#define SUPER(dir) \
static int super_##dir (PATNODE *node,                                 \
                        const ITEM *items, ITEM n, SUPP supp)          \
{                               /* --- check for a superset */         \
  ITEM i;                       /* loop variable */                    \
                                                                       \
  assert(items && (n > 0) && (supp > 0));   /* check arguments */      \
  while (node                   /* while there is another node with */ \
  &&    !dir(*items, node->items[0])) {  /* a first item not after */  \
    if ((*items != node->items[0])       /* the next item to match */  \
    &&  super_##dir(node->sibling, items, n, supp))                    \
      return -1;                /* process siblings of current node */ \
    if (node->supp < supp)      /* if the support is too low, */       \
      return 0;                 /* no superset will be found */        \
    for (i = 0; (n > 0) && (i < node->cnt); i++)                       \
      if (node->items[i] == *items) {                                  \
        ++items; --n; }         /* check other items in node */        \
    if (n <= 0) return -1;      /* if all items were matched, abort */ \
    node = node->children;      /* continue with the child nodes */    \
  }                             /* (match the remaining items) */      \
  return 0;                     /* return 'no superset exists' */      \
}  /* super() */

/*--------------------------------------------------------------------*/

SUPER(pos)                      /* function for ascending  item order */
SUPER(neg)                      /* function for descending item order */

/*--------------------------------------------------------------------*/

int pat_super (PATTREE *pat, const ITEM *items, ITEM n, SUPP supp)
{                               /* --- check for a superset */
  assert(pat                    /* check the function arguments */
  &&    (items || (n <= 0)) && (supp > 0));
  if (n <= 0)                   /* if no item, check root support */
    return (pat->root.supp >= supp) ? -1 : 0;
  return (pat->dir < 0)         /* recursively check for a superset */
       ? super_neg(pat->root.children, items, n, supp)
       : super_pos(pat->root.children, items, n, supp);
}  /* pat_super() */

/*--------------------------------------------------------------------*/

#define MERGE(dir) \
static PATNODE* merge_##dir (PATNODE *s1, PATNODE *s2, PATTREE *pat)   \
{                               /* --- merge two node list */          \
  ITEM    i, k;                 /* loop variables */                   \
  PATNODE *out, **end;          /* output node list and end pointer */ \
  PATNODE *node;                /* split and deallocation buffer */    \
                                                                       \
  if (!s1) return s2;           /* if there is only one node list, */  \
  if (!s2) return s1;           /* simply return the other list */     \
  end = &out;                   /* start the output list */            \
  while (1) {                   /* node list merge loop */             \
    if      (dir(s1->items[0], s2->items[0])) {                        \
      *end = s1; end = &s1->sibling; if (!(s1 = *end)) break; }        \
    else if (dir(s2->items[0], s1->items[0])) {                        \
      *end = s2; end = &s2->sibling; if (!(s2 = *end)) break; }        \
    else {                      /* copy nodes with singular items */   \
      k = (s1->cnt < s2->cnt) ? s1->cnt : s2->cnt;                     \
      for (i = 1; i < k; i++)   /* count the equal items */            \
        if (s1->items[i] != s2->items[i]) break;                       \
      if (i < s1->cnt) {        /* if not all items were matched, */   \
        node = split(s1, i);    /* split node from first list */       \
        if (!node) {            /* if the node split failed */         \
          pat->err = -1; *end = s1; end = &s1->sibling;                \
          if (!(s1 = *end)) break; else continue;                      \
        }                       /* keep node and set error flag */     \
        s1 = node; pat->cnt++;  /* get top part of split node */       \
      }                         /* and count the added node */         \
      if (s1->supp < s2->supp)  /* matched items: keep the maximum */  \
        s1->supp = s2->supp;    /* of the support values */            \
      if (i < s2->cnt) {        /* if not all items were matched, */   \
        s2->cnt -= i;           /* move the remaining items */         \
        memmove(s2->items, s2->items+i, (size_t)s2->cnt*sizeof(ITEM)); \
        node = (PATNODE*)realloc(s2, PATNODESIZE(s2->cnt));            \
        s2   = node->sibling; node->sibling = NULL;                    \
        s1->children = merge_##dir(s1->children, node, pat); }         \
      else {                    /* if all items were matched */        \
        s1->children = merge_##dir(s1->children, s2->children, pat);   \
        node = s2; s2 = s2->sibling; free(node); pat->cnt--;           \
      }                         /* merge and delete second node */     \
      *end = s1; end = &s1->sibling; s1 = *end;                        \
      if (!s1 || !s2) break;    /* if an item occurs in both lists, */ \
    }                           /* rec. merge the lists of children */ \
  }                             /* and delete one of the two nodes */  \
  *end = (s1) ? s1 : s2;        /* append the remaining nodes */       \
  return out;                   /* return the merged patricia tree */  \
}  /* merge() */

/*--------------------------------------------------------------------*/

MERGE(pos)                      /* function for ascending  item order */
MERGE(neg)                      /* function for descending item order */

/*--------------------------------------------------------------------*/
#if 1                           /* full reduction of items */

#define SORT(dir) \
static PATNODE* sort_##dir (PATNODE *list, PATTREE *pat)               \
{                               /* --- sort a transaction list */      \
  PATNODE *b, *a = list;        /* to traverse the lists */            \
                                                                       \
  assert(list && pat);          /* check the function arguments */     \
  for (b = list->sibling; b; ){ /* traverse the list to sort */        \
    b = b->sibling;             /* two steps on b, one step on list */ \
    if (b) { b = b->sibling; list = list->sibling; }                   \
  }                             /* (split list into two halves) */     \
  b = list->sibling;            /* get the second list and */          \
  list->sibling = NULL;         /* terminate the first list */         \
  if (a->sibling) a = sort_##dir(a, pat); /* if more than 1 elem., */  \
  if (b->sibling) b = sort_##dir(b, pat); /* sort lists recursively */ \
  return merge_##dir(a,b, pat); /* return the merged lists */          \
}  /* sort() */

/*--------------------------------------------------------------------*/

SORT(pos)                       /* function for ascending  item order */
SORT(neg)                       /* function for descending item order */

/*--------------------------------------------------------------------*/

#define PRUNEX(dir) \
static void prunex_##dir (PATNODE *node, PATTREE *pat)                 \
{                               /* --- prune a patricia tree */        \
  int     s = 0;                /* flag for item sorting */            \
  ITEM    i, k;                 /* loop variables, item counter */     \
  PATNODE *n, *t, *x;           /* to traverse the node list */        \
  PATNODE *o, **p = &o;         /* nodes to keep and end pointer */    \
                                                                       \
  assert(node && pat);          /* check the function arguments */     \
  n = node->children; node->children = NULL;                           \
  while (n) {                   /* traverse the sibling list */        \
    if (n->children)            /* recursively prune the children */   \
      prunex_##dir(n, pat);     /* before pruning the node itself */   \
    s |= (n->supp >= pat->mins[n->items[0]]);                          \
    for (i = k = 0; i < n->cnt; i++)                                   \
      if (n->supp >= pat->mins[n->items[i]])                           \
        n->items[k++] = n->items[i];                                   \
    if (k > 0) {                /* collect items that may be freq. */  \
      while ((t = n->children)  /* while merger with child possible */ \
      &&     !t->sibling && (n->supp <= t->supp)) {                    \
        x = (PATNODE*)realloc(n, PATNODESIZE(k +t->cnt));              \
        if (!x) break;          /* try to enlarge the node */          \
        n = x;                  /* on success get the new node */      \
        memcpy(n->items +k, t->items, (size_t)t->cnt *sizeof(ITEM));   \
        n->cnt = k += t->cnt;   /* append the child node items */      \
        n->children = t->children;                                     \
        free(t); pat->cnt--;    /* unlink and delete the child node */ \
      }                         /* and decrement the node counter */   \
      if (k < n->cnt) {         /* if the node is too large, */        \
        n->cnt = k;             /* set the new number of items */      \
        n = (PATNODE*)realloc(n, PATNODESIZE(k));                      \
      }                         /* try to shrink the node */           \
      *p = n; p = &n->sibling;  /* move the current node */            \
      n = *p; }                 /* to the output node list */          \
    else {                      /* if item (set) is infrequent */      \
      node->children =          /* (and cannot become frequent) */     \
        merge_##dir(node->children, n->children, pat);                 \
      t = n; n = n->sibling;    /* merge the child nodes */            \
      free(t); pat->cnt--;      /* with the pruned subtrees */         \
    }                           /* and delete the processed node */    \
  }                                                                    \
  *p = NULL;                    /* finally terminate the keep list, */ \
  if (s && o && o->sibling)     /* sort the keep list if necessary, */ \
    o = sort_##dir(o,pat);      /* and merge the two output lists */   \
  node->children = merge_##dir(node->children, o, pat);                \
}  /* prunex() */

/*--------------------------------------------------------------------*/
#else                           /* partial reduction of items */

#define PRUNEX(dir) \
static void prunex_##dir (PATNODE *node, PATTREE *pat)                 \
{                               /* --- prune a patricia tree */        \
  ITEM    i, k;                 /* loop variables, item counter */     \
  PATNODE *n, *t, *x;           /* to traverse the node list */        \
  PATNODE *o, **p = &o;         /* nodes to keep and end pointer */    \
                                                                       \
  assert(node && pat);          /* check the function arguments */     \
  n = node->children; node->children = NULL;                           \
  while (n) {                   /* traverse the sibling list */        \
    if (n->children)            /* recursively prune the children */   \
      prunex_##dir(n, pat);     /* before pruning the node itself */   \
    for (i = k = 1; i < n->cnt; i++)                                   \
      if (n->supp >= pat->mins[n->items[i]])                           \
        n->items[k++] = n->items[i];                                   \
    if (k <= 1)                 /* collect the qualifying items */     \
      k = (n->supp >= pat->mins[n->items[0]]) ? 1 : 0;                 \
    /* First item can be removed only if the node becomes empty. */    \
    /* Otherwise the output list needs sorting and node merging. */    \
    if (k > 0) {                /* collect items that may be freq. */  \
      while ((t = n->children)  /* while merger with child possible */ \
      &&     !t->sibling && (n->supp <= t->supp)) {                    \
        x = (PATNODE*)realloc(n, PATNODESIZE(k +t->cnt));              \
        if (!x) break;          /* try to enlarge the node */          \
        n = x;                  /* on success get the new node */      \
        memcpy(n->items +k, t->items, (size_t)t->cnt *sizeof(ITEM));   \
        n->cnt = k += t->cnt;   /* append the child node items */      \
        n->children = t->children;                                     \
        free(t); pat->cnt--;    /* unlink and delete the child node */ \
      }                         /* and decrement the node counter */   \
      if (k < n->cnt) {         /* if the node is too large, */        \
        n->cnt = k;             /* set the new number of items */      \
        n = (PATNODE*)realloc(n, PATNODESIZE(k));                      \
      }                         /* try to shrink the node */           \
      *p = n; p = &n->sibling;  /* move the current node */            \
      n = *p; }                 /* to the output node list */          \
    else {                      /* if item (set) is infrequent */      \
      node->children =          /* (and cannot become frequent) */     \
        merge_##dir(node->children, n->children, pat);                 \
      t = n; n = n->sibling;    /* merge the child nodes */            \
      free(t); pat->cnt--;      /* with the pruned subtrees */         \
    }                           /* and delete the processed node */    \
  }                             /* finally terminate the keep list */  \
  *p = NULL;                    /* and merge the two output lists */   \
  node->children = merge_##dir(node->children, o, pat);                \
}  /* prunex() */

#endif
/*--------------------------------------------------------------------*/

PRUNEX(pos)                     /* function for ascending  item order */
PRUNEX(neg)                     /* function for descending item order */

/*--------------------------------------------------------------------*/

int pat_prunex (PATTREE *pat, SUPP supp, const SUPP *frqs)
{                               /* --- prune infrequent item sets */
  ITEM i;                       /* loop variable */

  assert(pat && (supp > 0) && frqs);  /* check the function arguments */
  for (i = 0; i < pat->size;i++)/* traverse the items and */
    pat->mins[i] = supp-frqs[i];/* compute minimum support values */
  pat->err = 0;                 /* recursively prune the tree */
  if (pat->dir < 0) prunex_neg(&pat->root, pat);
  else              prunex_pos(&pat->root, pat);
  if (pat->cnt > pat->max)      /* update the maximum number of nodes */
    pat->max = pat->cnt;        /* (just to be sure that max >= cnt) */
  return pat->err;              /* return the error status */
}  /* pat_prunex() */

/*--------------------------------------------------------------------*/

static void prune (PATNODE **node, SUPP supp, PATTREE *pat)
{                               /* --- recursively prune tree */
  PATNODE *t;                   /* temporary buffer for deletion */

  while (*node) {               /* traverse the sibling list */
    if ((*node)->children)      /* prune children recursively */
      prune(&(*node)->children, supp, pat);
    if ((*node)->supp >= supp){ /* keep nodes with sufficient supp. */
      node = &(*node)->sibling; continue; }
    t = *node; *node = (*node)->sibling;
    free(t); pat->cnt--;        /* nodes with insufficient support */
  }                             /* are removed from sibling list */
}  /* prune() */

/*--------------------------------------------------------------------*/

void pat_prune (PATTREE *pat, SUPP supp)
{                               /* --- prune item set repository */
  assert(pat && (supp >= 0));   /* check the function arguments */
  prune(&pat->root.children, supp, pat);
}  /* pat_prune() */            /* recursively prune prefix tree */

/*--------------------------------------------------------------------*/

static int closed (PATTREE *pat, PATNODE *node)
{                               /* --- report closed item sets */
  int  r, x = 0;                /* error status, perfect ext. flag */
  ITEM i;                       /* loop variable */
  SUPP supp;                    /* support of current item set */

  assert(pat && node);          /* check the function arguments */
  supp = isr_supp(pat->rep);    /* get current item set support */
  if (isr_xable(pat->rep, 1)) { /* if item set may be extended */
    for (node = node->children; node; node = node->sibling) {
      if (node->supp < pat->min)/* traverse the node nodes, */
        continue;               /* but skip infrequent item sets */
      x |= (node->supp >= supp);/* set perfect extension flag */
      for (i = 0; i < node->cnt; i++) {
        if (!isr_xable(pat->rep, 1)) break;
        r = isr_addnc(pat->rep, node->items[i], node->supp);
        if (r < 0) return r;    /* add all items in the node */
      }                         /* to the item set reporter */
      r = closed(pat, node);    /* recursively report item sets */
      isr_remove(pat->rep, i);  /* remove the items of the node */
      if (r < 0) return r;      /* from the item set reporter, */
    } }                         /* finally check for an error */
  else {                        /* if item set may not be extended */
    for (node = node->children; node; node = node->sibling)
      if (node->supp >= supp) { x = -1; break; }
  }                             /* check for a perfect extension */
  return (x) ? 0 : isr_report(pat->rep);
}  /* closed() */               /* report sets without perfect exts. */

/*--------------------------------------------------------------------*/

static int maximal (PATTREE *pat, PATNODE *node)
{                               /* --- report maximal item sets */
  int  r, x = 0;                /* error status, freq. child flag */
  ITEM i;                       /* loop variable */

  assert(pat && node);          /* check the function arguments */
  if (isr_xable(pat->rep, 1)) { /* if item set may be extended */
    for (node = node->children; node; node = node->sibling) {
      if (node->supp < pat->min)/* traverse the child nodes, */
        continue;               /* but skip infrequent item sets */
      x = -1;                   /* set flag for a frequent child */
      for (i = 0; i < node->cnt; i++) {
        if (!isr_xable(pat->rep, 1)) break;
        r = isr_addnc(pat->rep, node->items[i], node->supp);
        if (r < 0) return r;    /* add all items in the node */
      }                         /* to the item set reporter */
      r = maximal(pat, node);   /* recursively report item sets */
      isr_remove(pat->rep, i);  /* and remove the current item */
      if (r < 0) return r;      /* from the item set reporter */
    } }                         /* finally check for an error */
  else {                        /* if item set may not be extended */
    for (node = node->children; node; node = node->sibling)
      if (node->supp >= pat->min) { x = -1; break; }
  }                             /* check for a frequent superset */
  return (x) ? 0 : isr_report(pat->rep);
}  /* maximal() */              /* report sets w/o frequent supersets */

/*--------------------------------------------------------------------*/

static int maxonly (PATTREE *pat, PATNODE *node)
{                               /* --- report maximal item sets */
  int        r, x = 0;          /* error status, freq. child flag */
  ITEM       i;                 /* loop variable, number of items */
  const ITEM *s;                /* current item set */
  PATNODE    *curr = node;      /* buffer for current node */

  assert(pat && node);          /* check the function arguments */
  if (isr_xable(pat->rep, 1)) { /* if item set may be extended */
    for (node = node->children; node; node = node->sibling) {
      if (node->supp < pat->min)/* traverse the child nodes, */
        continue;               /* but skip infrequent item sets */
      x = -1;                   /* set flag for a frequent child */
      for (i = 0; i < node->cnt; i++) {
        if (!isr_xable(pat->rep, 1)) break;
        r = isr_addnc(pat->rep, node->items[i], node->supp);
        if (r < 0) return r;    /* add all items in the node */
      }                         /* to the item set reporter */
      r = maxonly(pat, node);   /* recursively report item sets */
      isr_remove(pat->rep, i);  /* and remove the current item */
      if (r < 0) return r;      /* from the item set reporter */
    } }                         /* finally check for an error */
  else {                        /* if item set may not be extended */
    for (node = node->children; node; node = node->sibling)
      if (node->supp >= pat->min) { x = -1; break; }
  }                             /* check for a frequent superset */
  if (x) return 0;              /* if there is a freq. child, abort */
  curr->supp = -curr->supp;     /* mark the current node as to ignore */
  i = isr_cnt  (pat->rep);      /* get the current number of items */
  s = isr_items(pat->rep);      /* and the current item set and */
  x = (pat->dir < 0)            /* recursively search for a superset */
    ? super_neg(pat->root.children, s, i, pat->min)
    : super_pos(pat->root.children, s, i, pat->min);
  curr->supp = -curr->supp;     /* unmark the current node */
  return (x) ? 0 : isr_report(pat->rep);
}  /* maxonly() */              /* report sets w/o frequent supersets */

/*--------------------------------------------------------------------*/

int pat_report (PATTREE *pat, int max, SUPP supp, ISREPORT *rep)
{                               /* --- report (closed) item sets */
  assert(pat && rep);           /* check the function arguments */
  pat->min = supp;              /* note the minimum support and */
  pat->rep = rep;               /* the item set reporter */
  if      (max < 0)             /* if maximal item sets only */
    return maxonly(pat, &pat->root);
  else if (max > 0)             /* if maximal item sets (ext. filter) */
    return maximal(pat, &pat->root);
  else                          /* if closed  item sets */
    return closed (pat, &pat->root);
}  /* pat_report() */

/*--------------------------------------------------------------------*/
#ifndef NDEBUG

static void indent (int k)
{ while (--k >= 0) printf("   "); }

/*--------------------------------------------------------------------*/

static void show (PATNODE *node, ITEMBASE *base, int ind)
{                               /* --- recursively show nodes */
  ITEM i;                       /* loop variable */

  assert(ind >= 0);             /* check the function arguments */
  while (node) {                /* traverse the node list */
    indent(ind);                /* indent the output line */
    for (i = 0; i < node->cnt; i++) { /* traverse the items */
      if (base) printf("%s/", ib_name(base, node->items[i]));
      printf("%"ITEM_FMT" ", node->items[i]);
    }                           /* print the node information */
    printf(": %"SUPP_FMT"\n", node->supp);
    show(node->children, base, ind+1);
    node = node->sibling;       /* recursively show the child nodes, */
  }                             /* then go to the next node */
}  /* show() */

/*--------------------------------------------------------------------*/

void pat_show (PATTREE *pat, ITEMBASE *base)
{                               /* --- print a patricia tree */
  assert(pat);                  /* check the function arguments */
  if (!pat) {                   /* check whether tree exists */
    printf("(null)\n"); return; }
  show(pat->root.children, base, 0); /* recursively show the nodes */
  printf("supp:  %"SUPP_FMT"\n", pat->root.supp);
  printf("nodes: %"SIZE_FMT"\n", pat->cnt);
}  /* pat_show() */             /* print global information */

#endif
