
namespace LaPackCPP {

static const double log10e = 0.43429448190325182765;

// todo: templatize all these functions

template<class T>
T d_lg10(T *x)
{
  return( log10e * log(*x) );
}

template<class T>
T d_sign(T *a, T *b)
{
  T x;
  x = (*a >= 0 ? *a : - *a);
  return( *b >= 0 ? x : -x);
}

template<class T>
integer i_dnnt(T *x)
{
  return (integer)(*x >= 0. ? floor(*x + .5) : -floor(.5 - *x));
}

//template<class T>
//integer i_nint(T *x)
//{
//  return (integer)(*x >= 0 ? floor(*x + .5) : -floor(.5 - *x));
//}


template<class T>
T pow_di(T *ap, integer *bp)
{
  T pow, x;
  integer n;
  unsigned long u;

  pow = 1;
  x = *ap;
  n = *bp;

  if(n != 0)
  {
    if(n < 0)
    {
      n = -n;
      x = 1/x;
    }
    for(u = n; ; )
    {
      if(u & 01)
        pow *= x;
      if(u >>= 1)
        x *= x;
      else
        break;
    }
  }
  return(pow);
}

}