template<class T>
T rsEllipticSectorArea(T a, T b, T p)
{
  T A = 0.0;
  while(p >= 2*PI) // handle full turnarounds
  {
    A += PI*a*b;
    p -= 2*PI;
  }
  if(p >= PI)      // p is now in [0...2*pi[, handle lower half-ellipse
  {
    A += PI*a*b/2;
    p -= PI;
  }
  if(p >= PI/2)    // handle cases where angle is above 90 degrees
  {
    p -= PI/2;
    A += (a*b/2)*atan((b*tan(p))/a) + PI*a*b/4;
  }
  else               // angle is below 90 degrees
    A += (a*b/2)*atan((a*tan(p))/b);
  return A;
}
