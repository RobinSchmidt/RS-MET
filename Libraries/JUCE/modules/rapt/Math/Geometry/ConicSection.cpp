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
  
  // ...solve quadratic...(factor out and optimize):
  // Solutions: t_1,2 = (-b +- sqrt(b^2-4*a*c)) / (2*a):
  T d = b*b - 4*a*c;         // discriminant
  d   = sqrt(d);
  *t1 = (-b+d) / (2*a);        
  *t2 = (-b-d) / (2*a);
}

template<class T>
T rsConicSection<T>::evaluate(T x, T y)
{
  return A*x*x + B*x*y + C*y*y + D*x + E*y + F;
}

template<class T>
void rsConicSection<T>::getTangentCoeffs(T x, T y, T* a, T* b, T* c)
{
  *a = 2*A*x + B*y + D;
  *b = 2*C*y + B*x + E;
  *c = -(*a * x + *b * y);

  // formula derived from:
  // https://en.wikipedia.org/wiki/Implicit_curve#Tangent_and_normal_vector
}
