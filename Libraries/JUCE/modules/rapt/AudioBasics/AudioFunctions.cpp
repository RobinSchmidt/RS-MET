
template<class T>
T rsFindCosistentPhase(T target, T value) 
{
  T p = target;  // 0..2pi ..rename to t for target
  T q = value;   // 0..inf ..rename to v for value
  T tmp = p;
  T k   = 0;
  T d;
  while(true) {
    tmp = p + 2*k*PI;
    d   = q - tmp;
    if(rsAbs(d) <= PI)
      break;
    k  += T(1);
  }
  return tmp;
  // this may lead to an endless loop when 
  // storedPhase = 6.2804880025844536 (target), computedPhase = 0.0044408067227543532 (value)
  // -> make unit test - i think, the condition rsAbs(d) < PI is not enough ..maybe
  // it needs a tolerance? nope: using rsAbs(d) <= PI + 1.e-13 also doesn't work
  // maybe the check has to be done *before* adding 2*k*PI?

  // maybe generalize to an arbitrary range: 0..2pi, -pi..pi, 0..1, -1..1, etc. - the result 
  // should be such that abs(result - target) <= 0.5*range, or maybe <= (0.5*range)*(1+tol)
  // ..how large should tol be? is one machine epsilon enough? make unit tests that really
  // push the boundary ..but maybe, if the float resolution get too small at very high values,
  // that may fail, too - maybe wee need to additionally check that result <= value + margin 
  // for some margin based on the range? ...but what if the range is so small than 
  // value + margin = value after rounding?

  // rename to findConsistentUnwrappedValue(T preliminaryUnwrappedValue, T targetWrappedValue, 
  //  T range
}

template<class T>
T rsFindCosistentUnwrappedValue(T value, T target, T rangeMin, T rangeMax)
{
  T rangeSize = rangeMax - rangeMin;
  T maxDelta  = T(0.5) * rangeSize;
  T result    = target;
  T sign      = value > result ? T(-1) : T(+1);
  T k         = 0;
  while(true) {
    result  = target + sign*k*rangeSize;
    T delta = result - value;
    if( rsAbs(delta) <= maxDelta )
      break;
    k += T(1);
  }
  return result;
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