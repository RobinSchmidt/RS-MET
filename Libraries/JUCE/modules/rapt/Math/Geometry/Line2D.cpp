
template<class T>
void rsLine2D<T>::twoPointToImplicit(T x0, T y0, T x1, T y1, T& A, T& B, T& C)
{
  A = y0 - y1; // == -dy
  B = x1 - x0; // ==  dx
  C = -(x0*A + y0*B);
  T s = 1 / sqrt(A*A + B*B);
  A *= s;
  B *= s;
  C *= s;
}
