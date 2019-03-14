
/*
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
*/

template<class T>
T rsConsistentUnwrappedValue(T value, T target, T rangeMin, T rangeMax)
{
  T rangeSize = rangeMax - rangeMin;
  T maxDelta  = T(0.5) * rangeSize;  // maybe we should have a safety margin based on eps?

  //T margin    = RS_EPS(T);         // nope! hangs! ...
  T margin    = 1.e-13;              // ...but that's rather arbitrary!
  maxDelta   *= 1 + margin;          // should that depend on value, too? ..or maybe we should use 
                                     // as 2nd break condition that the sign of delta changes?

  T result    = target;
  T sign      = value > result ? T(+1) : T(-1);
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
// can we do something based on fmod? this function has linear complexitity in the size of the 
// passed value - the larger the value, the longer it takes - bad! ...but is fmod actually of O(1)
// complexity anyway?

template<class T>
T rsConsistentUnwrappedValue0(T value, T target, T rangeSize)
{
  T over     = fmod(value, rangeSize);             // overshoot over an integer number of cycles
  T cycles   = round((value - over) / rangeSize);  // number of full cycles through the range
  T result   = cycles * rangeSize + target;

  // re-adjustement:
  T delta    = value - result;
  T maxDelta = T(0.5)*rangeSize;
  if(delta > maxDelta)
    result += rangeSize;
  else if(delta < -maxDelta)
    result -= rangeSize;
  rsAssert(rsAbs(value - result) <= T(0.5)*rangeSize);  // self-check

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