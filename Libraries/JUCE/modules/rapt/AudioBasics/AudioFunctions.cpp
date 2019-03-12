

// useful for freq de-biasing and for phase interpolation in synthesis
template<class T>
T rsFindCosistentPhase(T storedPhase, T computedPhase) 
{
  T p = storedPhase;     // 0..2pi
  T q = computedPhase;   // 0..inf
  T tmp = p;
  T k   = 0;
  T d;
  while(true) {
    tmp = p + 2*k*PI;
    d   = q - tmp;
    if(rsAbs(d) <= PI + 1.e-13)
      break;
    k  += T(1);
  }
  return tmp;
  // this may lead to an endless loop when 
  // storedPhase = 6.2804880025844536, computedPhase = 0.0044408067227543532
  // -> make unit test - i think, the condition rsAbs(d) < PI is not enough ..maybe
  // it needs a tolerance? nope: using rsAbs(d) <= PI + 1.e-13 also doesn't work
}

// whoa - this is very tricky - isn't there a simpler way for this?
template<class T>
T rsPhaseError(T p1, T p2)
{
  p1  = RAPT::rsWrapToInterval(p1, 0, 2*PI);  // 0..2pi  // rename to rsWrapToRange...or just rsWrap
  p2  = RAPT::rsWrapToInterval(p2, 0, 2*PI);  // 0..2pi
  T d = p2-p1;                                // -2pi..2pi
  d   = RAPT::rsWrapToInterval(d,  0, 2*PI);  // 0..2pi
  return d;
}

template<class T>
bool rsArePhasesConsistent(T p1, T p2, T tol)
{
  T d = rsPhaseError(p1, p2);
  if(d < tol)
    return true;
  if(d-2*PI < tol)
    return true;
  return false;
}