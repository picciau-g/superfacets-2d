#ifndef SORTING_H
#define SORTING_H

#include <math.h>

typedef struct {
    int v1,v2;
    int t;
} aux;

static int cmp_aux( const void *p, const void *q )
{
  aux *a = (aux *) p ;
  aux *b = (aux *) q ;
  if ( a->v1 > b->v1 ) return 1 ;
  if ( a->v1 < b->v1 ) return -1 ;
  if ( a->v2 > b->v2 ) return 1;
  if ( a->v2 < b->v2 ) return -1;
  return 0;
}

static bool CompareAux(aux p, aux q)
{
    if(p.v1 < q.v1)
        return true;
    if(p.v1 == q.v1 && p.v2 < q.v2)
        return true;
    if(p.v1 == q.v1 && p.v2 == q.v2 && p.t < q.t)
        return true;

    return false;

}

#endif // SORTING_H
