template<class T>
rsConicSection<T>::rsConicSection(T _A, T _B, T _C, T _D, T _E, T _F)
  : A(_A), B(_B), C(_C), D(_D), E(_E), F(_F)
{

}

template<class T>
void rsConicSection<T>::lineIntersectionParameter(T x, T dx, T y, T dy, T* t1, T* t2) const
{
  // Coeffs of quadratic equation: a*t^2 + b*t + c = 0:
  T a = A*dx*dx + B*dx*dy + C*dy*dy;
  T b = 2*A*x*dx + B*(x*dy+y*dx) + 2*C*y*dy + D*dx + E*dy;
  T c = A*x*x + B*x*y + C*y*y + D*x + E*y + F;
  rsPolynomial<T>::rootsQuadraticReal(c, b, a, t1, t2);
}

template<class T>
T rsConicSection<T>::evaluate(T x, T y) const
{
  return A*x*x + B*x*y + C*y*y + D*x + E*y + F;
}

template<class T>
void rsConicSection<T>::getTangentCoeffs(T x, T y, T* a, T* b, T* c) const
{
  *a = 2*A*x + B*y + D;
  *b = 2*C*y + B*x + E;
  *c = -(*a * x + *b * y);

  // formula derived from:
  // https://en.wikipedia.org/wiki/Implicit_curve#Tangent_and_normal_vector
}

/*

see this video for conics that go through given points or touch given lines and more:
https://www.youtube.com/watch?v=X83vac2uTUs
https://github.com/HackerPoet/Conics   code to accompany the video

*/