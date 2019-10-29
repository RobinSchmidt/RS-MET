template<class T>
void rsLine2D<T>::twoPointToImplicit(T x0, T y0, T x1, T y1, T* A, T* B, T* C, bool normalize)
{
  *A = y0 - y1; // == -dy
  *B = x1 - x0; // ==  dx
  *C = -(x0 * *A + y0 * *B);
  if(normalize)
  {
    T s = 1 / sqrt(*A * *A + *B * *B);
    *A *= s;
    *B *= s;
    *C *= s;
  }
}

template<class T>
void rsLine2D<T>::reflectPointInLine(T x, T y, T A, T B, T C, T *xr, T *yr)
{
  T d = 2*(A*x + B*y + C) / (A*A + B*B);
  *xr = x - A*d;
  *yr = y - B*d;

  // formula taken from:
  // https://math.stackexchange.com/questions/1013230/how-to-find-coordinates-of-reflected-point
  // its also in Salomon's "Computer Graphics...", page 71, Eq. 3.12

  // maybe make a version that assumes A^2 + B^2 = 1 so we can get rid of the division
}
