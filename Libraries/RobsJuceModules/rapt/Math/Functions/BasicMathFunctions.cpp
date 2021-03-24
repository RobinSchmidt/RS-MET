//template <class T>
//T rsAbs(T x)
//{
//  if(x < rsZeroValue(x))
//    return -x;
//  else
//    return  x;
//}
//inline double  rsAbs(double  x) { return fabs(x); }
//inline float   rsAbs(float   x) { return fabs(x); }
//inline rsInt8  rsAbs(rsInt8  x) { return  abs(x); }
//inline rsInt16 rsAbs(rsInt16 x) { return  abs(x); }
//inline rsInt32 rsAbs(rsInt32 x) { return  abs(x); }
//inline rsInt64 rsAbs(rsInt64 x) { return  abs(x); } // doesn't work with MinGW gcc 4.7

// for fast sign function, see here:
//http://stackoverflow.com/questions/1903954/is-there-a-standard-sign-function-signum-sgn-in-c-c

/*
template<class T>
T rsClip(T x, T min, T max)
{
  rsAssert(min <= max);
  if(x < min)
    return min;
  if(x > max)
    return max;
  return x;
}
*/

//inline rsInt32 rsClippedSum(rsInt32 a, rsInt32 b, rsInt32 min, rsInt32 max)
//{
//  rsInt64 tmp = (rsInt64)a + (rsInt64)b;
//  if(tmp < min)
//    return (rsInt32)min;
//  else if(tmp > max)
//    return (rsInt32)max;
//  else return (rsInt32)tmp;
//}
//
//inline rsInt32 rsClippedDifference(rsInt32 a, rsInt32 b, rsInt32 min, rsInt32 max)
//{
//  rsInt64 tmp = (rsInt64)a - (rsInt64)b;
//  if(tmp < min)
//    return (rsInt32)min;
//  else if(tmp > max)
//    return (rsInt32)max;
//  else return (rsInt32)tmp;
//}

inline bool rsExOr(bool a, bool b)
{
  return (a && !b) || (!a && b);
}

/*
inline double rsExpToLin(double in, double inMin, double inMax, double outMin, double outMax)
{
  double tmp = log(in / inMin) / log(inMax / inMin);
  return outMin + tmp * (outMax - outMin);
}

inline double rsExpToLinWithOffset(double in, double inMin, double inMax, double outMin,
  double outMax, double offset)
{
  double tmp = in * (inMax + offset) / inMax;
  tmp -= offset;
  return rsExpToLin(tmp, inMin, inMax, outMin, outMax);

  //double tmp = linToExp(in, inMin, inMax, outMin, outMax);
  //tmp += offset;
  //tmp *= outMax/(outMax+offset);
  //return tmp;
}
*/

inline bool rsIsCloseTo(double x, double targetValue, double tolerance)
{
  return rsAbs(x - targetValue) <= tolerance;
  /*
  if(fabs(x - targetValue) <= tolerance)
    return true;
  else
    return false;
    */
}
/*
template<class T>
inline bool rsIsCloseTo(T x, T targetValue, double tolerance)
{
  if(rsAbs(x - targetValue) <= tolerance)
    return true;
  else
    return false;
}
*/
/*
template<class T>
inline bool rsIsEven(T x)
{
  return x % 2 == 0;
}

template<class T>
inline bool rsIsOdd(T x)
{
  return x % 2 != 0;
}
*/

inline bool rsIsPowerOfTwo(unsigned int x)
{
  unsigned int currentPower = 1;
  while(currentPower <= x)
  {
    if(currentPower == x)
      return true;
    currentPower *= 2;
  }
  return false;
}

/*
template<class T>
inline bool rsIsInRange(T x, T min, T max)
{
  if(x >= min && x <= max)
    return true;
  else
    return false;
}
*/

inline int rsIntAndFracPart(double x, double &frac)
{
  int i = (int)x;
  frac  = x - (double)i;
  return i;
  // this version takes roughly 40 cycles - use reference for i, too - avoid local variable

  /*
  double f = floor(x);
  frac     = x - f;
  return rsRoundToInt(f);
  // this version takes roughly 50 cycles
  */
}

 // \todo supeseded by rsPow - get rid of it:
inline double rsIntegerPower(double x, int exponent)
{
  double accu = 1.0;
  for(int i=0; i<exponent; i++)
    accu *= x;
  return accu;
}

/*
// redundant with rsClip
template <class T>
inline T rsLimit(T x, T lowerBound, T upperBound)
{
  if(x < lowerBound)
    return lowerBound;
  if(x > upperBound)
    return upperBound;
  return x;
}
*/

/*
inline double rsLinToLin(double in, double inMin, double inMax, double outMin, double outMax)
{
  // map input to the range 0.0...1.0:
  double tmp = (in - inMin) / (inMax - inMin);

  // map the tmp-value to the range outMin...outMax:
  tmp *= (outMax - outMin);
  tmp += outMin;

  return tmp;
}
*/
/*
inline double rsLinToExp(double in, double inMin, double inMax, double outMin, double outMax)
{
  // map input to the range 0.0...1.0:
  double tmp = (in - inMin) / (inMax - inMin);

  // map the tmp-value exponentially to the range outMin...outMax:
  return outMin * std::exp(tmp * (log(outMax / outMin)));
}

inline double rsLinToExpWithOffset(double in, double inMin, double inMax, double outMin,
  double outMax, double offset)
{
  double tmp = rsLinToExp(in, inMin, inMax, outMin, outMax);
  tmp += offset;
  tmp *= outMax / (outMax + offset);
  return tmp;
}
*/

template <class T>
inline T rsNextEvenNumber(T x)
{
  if(rsIsEven(x))
    return x;
  return x+1;
}

/*
template <class T>
inline T rsNextPowerOfTwo(T x)
{
  T accu = 1;
  while(accu < x)
    accu *= 2;
  return accu;
}
*/

template <class T>
T rsPow(const T& base, int exponent)
{
  T result = rsUnityValue(base);
  T square(base);
  while(true)
  {
    if(exponent & 1)   // \todo use if( rsIsOdd(exponent) ), such that we don't rely on
      result *= square;  // bit-operations and then let exponent by of type T, too
    exponent /= 2;
    if(exponent == 0)
      break;
    square *= square;
  }
  return result;
}

inline double rsRandomUniform(double min, double max, int seed)
{
  static rsUint32 state = 0;
  if(seed >= 0)
    state = seed;
  state = 1664525 * state + 1013904223;
  return min + (max - min) * ((1.0 / 4294967296.0) * state);
}

template<class T>
inline void rsRangeConversionCoefficients(T inMin, T inMax, T outMin, T outMax, T *scale, T *shift)
{
  *scale = (outMax-outMin) / (inMax-inMin);
  *shift = outMin - (*scale * inMin);
}

//template <class T>
//inline int rsRoundToInt(T x)
//{
//  T xFloor = floor(x);
//  T xFrac  = x-xFloor;
//  if( xFrac >= 0.5 )
//    return (int) xFloor + 1;
//  else
//    return (int) xFloor;
//}

template <class T>
inline T rsSquare(T x)
{
  return x*x;
}

//template <class T>
//inline void rsSwap(T &in1, T &in2)
//{
//  T tmp = in1;
//  in1   = in2;
//  in2   = tmp;
//}

//template<class T>
//inline T rsUnityValue(T value)
//{
//  return T(1);
//}

template <class T>
inline void rsSortAscending(T &in1, T &in2)
{
  if(in1 > in2)
    rsSwap(in1, in2);
}

inline double rsWrapAround(double numberToWrap, double length)
{
  while(numberToWrap < 0.0)
    numberToWrap += length;
  return fmod(numberToWrap, length);
}

inline double rsWrapToInterval(double x, double min, double max)
{
  double r   = max-min;   // range
  double tmp = x-min;
  if(tmp >= 0.0)
    tmp = fmod(tmp, r);
  else
  {
    tmp = fmod(fabs(tmp), r);
    tmp = r-tmp;
  }
  return tmp + min;
}
// why no template? because of fmod?

inline double rsZeroFunction(double x)
{
  return x = 0.0;
}

//template<class T>
//inline T rsZeroValue(T /*value*/)
//{
//  return T(0);
//}

