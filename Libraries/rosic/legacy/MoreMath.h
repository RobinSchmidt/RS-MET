// maybe the command-line option /FORCE:MULTIPLE has to be added
// to the linker options

#ifndef MoreMath_h
#define MoreMath_h

#include <math.h>
#include <stdlib.h>
#include <float.h>
#include <emmintrin.h>
#include "Definitions.h"
//#include "Cephes/adapted/jn.c"  // bessel function from the cephes math library

namespace MoreMath
{

 INLINE double amp2dB(double amp);
 ///< Converts a raw amplitude value/factor to a value in decibels.

 INLINE int arrayMaxIndex(double* doubleArray, int numValues);
 /**< Returns the index of the maximum value in an array of doubles where the
      array should be of length numValues. */

 INLINE int arrayMinIndex(double* doubleArray, int numValues);
 /**< Returns the index of the minimum value in an array of doubles where the
      array should be of length numValues. */

 INLINE double belowOrAbove(double x, double low, double high);
 /**< Returns -1.0 if x is below low, 0.0 if x is between low and high and
      1.0 if x is above high. */   

 // bessel function n */
 INLINE double besselJ(double x, int n);

 INLINE double cheby(double x, double n);
 ///< Returns the value of the n-th chebychev polynomial.

 INLINE double ceil(double x);
 ///< Returns the next greater integer number.

 INLINE double centerClip(double x, double t);
 ///< Center clips the x at threshold t.

 INLINE double clip(double x, double min, double max);
 ///< Clips x at the threshold t.

 INLINE double dB2amp(double x);
 ///< Converts a value in decibels to a raw amplitude value/factor.

 INLINE double foldOver(double x, double min, double max);
 ///< foldover at the specified value */

 INLINE double frac(double x);
 ///< Returns the fractional part of x.

 INLINE double freqToPitch(double freq);
 ///< Converts a frequency in Hz into a MIDI-note value.

 // gamma function
 //double gamma(double x);

 INLINE double gauss(double x, double mu, double sigma);
 ///< Gaussian function with mean mu and standard deviation sigma.

 INLINE bool isEven(int x);
 INLINE bool isOdd(int x);
 INLINE bool isPowerOfTwo(unsigned int x);

 INLINE double log2(double x);
 ///< Calculates the logarithm to base 2.

 INLINE double logB(double x, double b);
 ///< Calculates logarithm to an arbitrary base b.

 INLINE double logistic(double x, double b);
 ///< Calculates the logistic function with slope parameter b.

 INLINE double mapLinearToLinear(double in, double inMin, double inMax,
                                 double outMin, double outMax);
 /**< Converts a value between inMin and inMax into a value between outMin
      and outMax where the mapping is linear for the input and the output.
      an example: y = mapLinearToLinear(x, 0.0, 1.0, -96.0, 24.0) will map
      the input x assumed to lie inside 0.0...1.0 to the range between
      -96.0...24.0. This function is useful to convert between parameter
      representations between 0.0...1.0 and the clear-text parameters. */

 INLINE double mapLinearToExponential(double in, double inMin, double inMax,
                                      double outMin, double outMax);
 /**< Converts a value between inMin and inMax into a value between outMin
      and outMax where the mapping of the output is exponential.
      an example: y = mapLinearToLinear(x, 0.0, 1.0, 20.0, 20000.0) will map
      the input x assumed to lie inside 0.0...1.0 to the range between
      20.0...20000.0 where equal differences in the input lead to equal 
      factors in the output. Make sure that the outMin value is greater than 
      zero! */

 INLINE double mapExponentialToLinear(double in, double inMin, double inMax,
                                      double outMin, double outMax);
 /**< The Inverse of "mapLinearToExponential" */


 //INLINE double max(double in1, double in2);
 /**< The maximum of two double numbers. */

 //INLINE int max(int in1, int in2);
 /**< The maximum of two int numbers. */

 INLINE double minkowski(double x, double mu, double sigma, double k);
 /**< Minkowski distribution - a generalization of the gaussian distribution,
      where the exponent k inside the exponential is an addtional parameter.
      For k=2, it reduces to the well known gaussian distribution. */

 INLINE double pitchOffsetToFreqFactor(double pitchOffset);
 /**< Converts a picth-offset in semitones value into a frequency 
      multiplication factor. */

 INLINE double pitchToFreq(double pitch);
 ///< Converts a MIDI-note value into a frequency in Hz.

 INLINE double quant(double x, double i);
 ///< Quantizes x to the interval i.

 INLINE double quantToBits(double x, double nBits);
 ///< Quantizes x (assumed to be in -1.0 - 1.0) to a grid which represents a resolution of nBits bits.

 INLINE double randDbl();
 ///< Generates a random number between -1.0...+1.0 in double format.

 INLINE double round(double x);
 ///< Returns the nearest integer (as double).

 INLINE int roundToInt(double x);
 ///< Returns the nearest integer.

 INLINE int roundToIntSSE2(double x);
 ///< Returns the nearest integer.

 INLINE double sawWave(double x);
 ///< Generates a 2*pi periodic saw wave.

 INLINE double sign(double x);
 ///< Returns the sign of x as double.

 INLINE void sinCos(double x, double* sinResult, double* cosResult);
 /**< Calculates sine and cosine of x more efficiently than calling sin(x) 
      and cos(x) seperately. */

 INLINE void sinCosApprox(double x, double* sinResult, double* cosResult);
 /**< Calculates a parabolic approximation of the sine and cosine of x. */

 INLINE double sqrWave(double x);
 ///< Generates a 2*pi periodic square wave.

 INLINE double step(double x);
 ///< The unit step function - returns 1 for x>=0, 0 else.

 INLINE double triWave(double x);
 ///< Generates a 2*pi periodic triangle wave.


 // wrapAround 



//---------------------------------------------------------------------------
 //implementation:

 INLINE double MoreMath::amp2dB(double amp)
 {
  return 20*log10(amp);
 }

 INLINE int arrayMaxIndex(double* doubleArray, int numValues)
 {
  int    maxIndex = 0;
  double maxValue = doubleArray[0];

  for(int i=0; i<numValues; i++)
  {
   if( doubleArray[i] > maxValue )
   {
    maxValue = doubleArray[i];
    maxIndex = i;
   }
  }

  return maxIndex;
 }

 INLINE int arrayMinIndex(double* doubleArray, int numValues)
 {
  int    minIndex = 0;
  double minValue = doubleArray[0];

  for(int i=0; i<numValues; i++)
  {
   if( doubleArray[i] < minValue )
   {
    minValue = doubleArray[i];
    minIndex = i;
   }
  }

  return minIndex;
 }


 INLINE double MoreMath::belowOrAbove(double x, double low, double high)
 {
  if( x < low )
   return -1.0;
  else if ( x > high )
   return 1.0;
  else
   return 0.0;
 }

 INLINE double cheby(double x, double n)
 {
  return cos( n * acos(x) );
 }

 INLINE double MoreMath::ceil(double x)
 {
  //double frac = x - floor(x);
  return x + (1-frac(x));  
 }

 INLINE double MoreMath::clip(double x, double min, double max)
 {
  if( x > max )
   return max;
  else if ( x < min )
   return min;
  else return x;
 }

 INLINE double dB2amp(double dB)
 {
  return pow(10.0, (0.05*dB));
 }

 INLINE double MoreMath::foldOver(double x, double min, double max)
 {
  if( x > max )
   return max - (x-max);
  else if( x < min )
   return min - (x-min);
  else return x;
 }

 INLINE double MoreMath::frac(double x)
 {
  return x - floor(x);
 }

 INLINE double MoreMath::freqToPitch(double freq)
 {
  return 12.0 * log2(freq/440.0) + 69.0;
 }

 INLINE double MoreMath::gauss(double x, double mu, double sigma)
 {
  return (1.0 / (sqrt(2*PI)*sigma)) * exp( - ((x-mu)*(x-mu)) / (2*sigma*sigma) );
 }

 INLINE bool MoreMath::isEven(int x)
 {
  if( x%2 == 0 )
   return true;
  else
   return false;
 }

 INLINE bool MoreMath::isOdd(int x)
 {
  if( x%2 != 0 )
   return true;
  else
   return false;
 }

 INLINE bool MoreMath::isPowerOfTwo(unsigned int x)
 {
  unsigned int currentPower    = 1;

  while( currentPower <= x )
  {
   if( currentPower == x )
    return true;

   currentPower *= 2;
  }

  return false;
 }

 INLINE double MoreMath::log2(double x)
 {
  return 1.4426950408889634073599246810019*log(x); // 1.44... = 1/log(2) 
 }

 INLINE double MoreMath::logB(double x, double b)
 {
  return log(x)/log(b);
 }

 INLINE double MoreMath::logistic(double x, double b)
 {
  return 1.0 / ( 1 + exp(-b*x) );
 }

 INLINE double mapLinearToLinear(double in, double inMin, double inMax,
                                 double outMin, double outMax)
 {
  double tmp;

  // map input to the range 0.0...1.0:
  tmp = (in-inMin) / (inMax-inMin);

  // map the tmp-value to the range outMin...outMax:
  tmp *= (outMax-outMin);
  tmp += outMin;

  return tmp;
 }

 INLINE double mapLinearToExponential(double in, double inMin, double inMax,
                                      double outMin, double outMax)
 {
  double tmp;

  // map input to the range 0.0...1.0:
  tmp = (in-inMin) / (inMax-inMin);

  // map the tmp-value exponentially to the range outMin...outMax:
  //tmp = outMin * exp( tmp*(log(outMax)-log(outMin)) );
  tmp = outMin * exp( tmp*(log(outMax/outMin)) );

  return tmp;
 }

 INLINE double mapExponentialToLinear(double in, double inMin, double inMax,
                                      double outMin, double outMax)
 {
  double tmp;

  tmp = log(in/inMin) / log(inMax/inMin);

  tmp = outMin + tmp * (outMax-outMin);

  return tmp;
 }

 /*
 INLINE double MoreMath::max(double in1, double in2)
 {
  if( in1 >= in2 )
   return in1;
  else
   return in2;
 }

 INLINE int MoreMath::max(int in1, int in2)
 {
  if( in1 >= in2 )
   return in1;
  else
   return in2;
 }
 */

 INLINE double MoreMath::minkowski(double x, double mu, double sigma, double k)
 {
   return (1 /(sqrt(2*PI)*sigma)) * exp( - ( pow(fabs(x-mu),k) ) / (2*pow(sigma,k)) );
   //return (1 /(sqrt(2*PI)*sigma)) * exp( - ( pow(abs(x-mu),k) ) / (2*pow(sigma,k)) );
   //return (1 /(sqrt(2*PI*pow(sigma,k)))) * exp( - ( pow(abs(x-mu),k) ) / (2*pow(sigma,k)) );
 }

 INLINE double pitchOffsetToFreqFactor(double pitchOffset)
 {
  return pow(2.0, pitchOffset/12.0);
 }

 INLINE double pitchToFreq(double pitch)
 {
  return 440.0*( pow(2.0, (pitch-69.0)/12.0) );
 }

 INLINE double MoreMath::quant(double x, double i)
 {
  double absX  = fabs(x);
  double signX = sign(x);
  double tmp   = fmod(x, i);
  if(x>=0)
  {
   if(tmp < 0.5*i)
    return x - tmp;
   else
    return x - tmp + i;
  }
  else
  {
   if(tmp < -0.5*i)
    return x - tmp - i;
   else
    return x - tmp;
  }
 }

 INLINE double MoreMath::quantToBits(double x, double nBits)
 {
  double interval = 2 / pow(2.0, nBits);
  return MoreMath::quant(x, interval);
 }

 INLINE double randDbl()
 {
  return (  ( (double)rand() - 16383.5 ) / 16383.5 );
 }

 INLINE double round(double x)
 {
  if( x-floor(x) >= 0.5 )
   return ceil(x);
  else
   return floor(x);
 }

 INLINE int roundToInt(double const x) 
 {
  int n;
  #if defined(__unix__) || defined(__GNUC__)
   // 32-bit Linux, Gnu/AT&T syntax:
   __asm ("fldl %1 \n fistpl %0 " : "=m"(n) : "m"(x) : "memory" );
  #else
   // 32-bit Windows, Intel/MASM syntax:
   __asm fld qword ptr x;
   __asm fistp dword ptr n;
  #endif
  return n;
 }

 INLINE int roundToIntSSE2(double const x) 
 {
  return _mm_cvtsd_si32(_mm_load_sd(&x));
 }

 INLINE double MoreMath::sign(double x)
 {
  if(x<0)
   return -1.0;
  else if(x>0)
   return 1.0;
  else
   return 0;
 }

 INLINE void sinCos(double x, double* sinResult, double* cosResult)
 {

  // this is the assembly-version:
  double s, c;
		__asm fld x
		__asm fsincos
		__asm fstp c
		__asm fstp s
  *sinResult = s;
  *cosResult = c;

  /*
  // this is the version, which fuzzes around with the square-root:
  *sinResult = sin(x);
  if( x < 0.5*PI )
   *cosResult = sqrt(1 - (*sinResult) * (*sinResult));
  else if ( x < (3.0/2.0)*PI )
   *cosResult = -sqrt(1 - (*sinResult) * (*sinResult));
  else
   *cosResult = sqrt(1 - (*sinResult) * (*sinResult));
   */

  /*
  // this is the unoptimized standard version:
  *sinResult = sin(x);
  *cosResult = cos(x);
  */
 }

 INLINE void sinCosApprox(double x, double* sinResult, double* cosResult)
 {
  static const double c = 0.70710678118654752440;

  // restrict input x to the range 0.0...2*PI:
  while( x > 2.0*PI )
   x -= 2*PI;
  while( x < 0.0 )
   x += 2*PI;

  if( x < PI/2 )
  {
   double tmp1 = x;
   double tmp2 = (2/PI) * tmp1 - 0.5;
   double tmp3 = (2-4*c)*tmp2*tmp2 + c;
   *sinResult  = tmp3 + tmp2;
   *cosResult  = tmp3 - tmp2;
  }
  else if( x < PI )
  {
   double tmp1 = (x-PI/2);
   double tmp2 = 0.5 - (2/PI) * tmp1;
   double tmp3 = (2-4*c)*tmp2*tmp2 + c;
   *sinResult  = tmp2 + tmp3;
   *cosResult  = tmp2 - tmp3;
  }
  else if( x < 1.5*PI )
  {
   double tmp1 = (x-PI);
   double tmp2 = (2/PI) * tmp1 - 0.5;
   double tmp3 = (4*c-2)*tmp2*tmp2 - c;
   *sinResult  = tmp3 - tmp2;
   *cosResult  = tmp3 + tmp2; 
  }
  else
  {
   double tmp1 = (x-1.5*PI);
   double tmp2 = (2/PI) * tmp1 - 0.5;
   double tmp3 = (2-4*c)*tmp2*tmp2 + c;
   *sinResult  = tmp2 - tmp3;
   *cosResult  = tmp2 + tmp3;  
  }

 }

 INLINE double MoreMath::sawWave(double x)
 {
  double tmp = fmod(x, 2*PI);
  if( tmp < PI )
   return tmp/PI;
  else
   return (tmp/PI)-2.0;
 }

 INLINE double MoreMath::sqrWave(double x)
 {
  double tmp = fmod(x, 2*PI);
  if( tmp < PI )
   return 1.0;
  else
   return -1.0;
 }

 INLINE double MoreMath::triWave(double x)
 {
  double tmp = fmod(x, 2*PI);
  if( tmp < 0.5*PI )
   return tmp/(0.5*PI);
  else if( tmp < 1.5*PI )
   return 1.0 - ((tmp-0.5*PI)/(0.5*PI));
  else
   return -1.0 + ((tmp-1.5*PI)/(0.5*PI));
 }


}

#endif // #ifndef MoreMath_h