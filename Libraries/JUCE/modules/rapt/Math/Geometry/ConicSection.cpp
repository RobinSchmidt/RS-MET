template<class T>
rsConicSection<T>::rsConicSection(T _A, T _B, T _C, T _D, T _E, T _F)
  : A(_A), B(_B), C(_C), D(_D), E(_E), F(_F)
{

}

template<class T>
void rsConicSection<T>::lineIntersectionParameter(T x, T dx, T y, T dy, T* t1, T* t2)
{
  // Coeffs of quadratic equation: a*t^2 + b*t + c = 0:
  T a = A*dx*dx + B*dx*dy + C*dy*dy;
  T b = 2*A*x*dx + B*(x*dy+y*dx) + 2*C*y*dy + D*dx + E*dy;
  T c = A*x*x + B*x*y + C*y*y + D*x + E*y + F;
  
  // something to do...solve quadratic...
}
