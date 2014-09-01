#ifndef SORTING_H
#define SORTING_H

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

#endif // SORTING_H
