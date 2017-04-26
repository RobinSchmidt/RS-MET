#ifndef rosic_ElementaryFunctionsReal_h
#define rosic_ElementaryFunctionsReal_h

// standard library includes:
#include <math.h>
#include <stdlib.h>

// rosic includes:
#include "../basics/GlobalFunctions.h"
#include "../basics/rosic_Constants.h"
#include "rosic_CephesDeclarations.h"  // move this over to rosic_SpecialFunctionsReal
#include "rosic_RealFunctionEvaluationAlgorithms.h"
#include "rosic_IntegerFunctions.h"
#include "../basics/rosic_FunctionTemplates.h"
#include "../basics/rosic_NumberManipulations.h"

namespace rosic
{

  /** This file contains functions for real arguments that are either themselves elementary functions (such as sin, exp, log, pow, tanh, 
  etc. - although, we don't duplicate these c-standard functions here) or can be expressed as combinations thereof (involving composition 
  and arithmetic operations), like gauss, log2, cheby, etc...  */

  // \todo maybe have special polynomials like cheby in an extra file - but maybe not
  // \todo move some of them into the cpp file

  /** A table of reciprocal values of the factorial of some integer number n, that is: 1/n!, where n should be between 0...31 (inclusive) 
  - the values are tabulated and any value outside the range 0...31 will be an access violation. These numbers are useful for making 
  Taylor series approximations. */
  static const int numInverseFactorials = 32;
  static long double inverseFactorials[numInverseFactorials] =
  {
    1.0,
    1.0,
    1.0/2.0,
    1.0/6.0,
    1.0/24.0,
    1.0/120.0,
    1.0/720.0,
    1.0/5040.0,
    1.0/40320.0,
    1.0/362880.0,
    1.0/3628800.0,
    1.0/39916800.0,
    1.0/479001600.0,
    1.0/6227020800.0,
    1.0/87178291200.0,
    1.0/1307674368000.0,
    1.0/20922789888000.0,
    1.0/355687428096000.0,
    1.0/6402373705728000.0,
    1.0/121645100408832000.0,
    1.0/2432902008176640000.0,
    1.0/51090942171709440000.0,
    1.0/1124000727777607680000.0,
    1.0/25852016738884976640000.0,
    1.0/620448401733239439360000.0,
    1.0/15511210043330985984000000.0,
    1.0/403291461126605635584000000.0,
    1.0/10888869450418352160768000000.0,
    1.0/304888344611713860501504000000.0,
    1.0/8841761993739701954543616000000.0,
    1.0/265252859812191058636308480000000.0,
    1.0/8222838654177922817725562880000000.0  // == 1.2161250415535179 * 10^-34
  };

  /** Inverse hyperbolic cosine. The argument must be >= 1, for the general case use the complex version acosh_c defined in 
  rosic_ComplexFunctions.h */
  INLINE double acosh(double x);

  /** Symmetrized inverse hyperbolic cosine. */
  INLINE double acoshs(double x);

  ///** Inverse hyperbolic sine. */
  //INLINE double asinh(double x);

  /** Returns -1.0 if x is below low, 0.0 if x is between low and high and 1.0 if x is above high. */
  INLINE double belowOrAbove(double x, double low, double high);

  /** Binomial coefficient for real n. */
  INLINE double binomialCoefficientReal(double n, int k);

  // Bessel function J_n(x) */
  //INLINE double besselj(double x, double order);

  /** Returns the value of the n-th chebychev polynomial. */
  INLINE double cheby(double x, double n);

  /** Returns the value of the n-th chebychev polynomial. */
  INLINE double cheby(double x, int n);

  /** Returns the next greater integer number. */
  //INLINE double ceil(double x);

  /** Center clips the x at threshold t. */
  INLINE double centerClip(double x, double t);

  /** Clips x into the range min...max. */
  template <class T>
  INLINE T clip(T x, T min, T max);

  /** Approximates cos^2(0.5*pi*x) over the domain -1.0...+1.0 via y = a*(1-x^2)^2 + (1-a)*(1-x^2)^3 with a == 0.565... (found visually - 
  may be refined via optimization). */
  INLINE double cosSquaredApprox(double x);

  /** Computes a sin/cos based equal power crossfade according to the value of x. The gain for the first signal will be stored in gain1, 
  the gain for the second signal in gain2. If the range of x is not 0.0...1.0, use optional the xMin and xMax arguments to specify the 
  range. For example, if x is a panning value between -100...+100, a constant power pan can be achieved via 
  equalPowerGainFactors(x, gainL, gainR, -100.0, 100.0), or if x is a dry/wet control between 0.0...100.0: 
  equalPowerGainFactors(x, gainDry, gainWet, 0.0, 100.0) */
  INLINE void equalPowerGainFactors(double x, double *gain1, double *gain2,
    double xMin = 0.0, double xMax = 1.0);

  /** Evaluates the quartic polynomial y = a4*x^4 + a3*x^3 + a2*x^2 + a1*x + a0 at x. */
  INLINE double evaluateQuartic(double x, double a0, double a1, double a2, double a3, double a4);

  /** foldover at the specified value */
  INLINE double foldOver(double x, double min, double max);

  /** Returns the fractional part of x. */
  INLINE double frac(double x);

  // gamma function
  //double gamma(double x);

  /** Gaussian function (aka normal distribution) with mean mu and standard deviation sigma. */
  INLINE double gauss(double x, double mu, double sigma);

  /** Splits a double into integer and fractional part, the integer part is returned in the
  return value, the fractional part in the second argument. */
  INLINE int intAndFracPart(double x, double &frac);

  /** Computes an integer power of x by successively multiplying x with itself. */
  INLINE double integerPower(double x, int exponent);

  /** Limits the value of an object.   ...redundant with clip */
  template <class T>
  INLINE T limit(T x, T lowerBound, T upperBound);

  /** Calculates the logistic function with slope parameter b. */
  INLINE double logistic(double x, double b);

  /** Fills the window-array with a Hamming-window. */
  INLINE void makeBlackmanWindow(double *window, int length);

  /** Fills the window-array with a cosine power window. */
  INLINE void makeCosinePowerWindow(double *window, int length, double power = 2.0);

  /** Fills the window-array with a Hamming-window. */
  INLINE void makeHammingWindow(double *window, int length);

  /** Fills the window-array with a Hanning-window. */
  INLINE void makeHanningWindow(double *window, int length);

  /** Fills the window-array with a rectangular window. ...ahem trivial */
  INLINE void makeRectangularWindow(double *window, int length);

  /** Minkowski distribution - a generalization of the gaussian distribution, where the exponent k inside the exponential is an addtional 
  parameter. For k=2, it reduces to the Gaussian distribution. k<2 makes the distribution more spikey and the tails fatter, values k>2 
  make the distribution more rectangular. */
  INLINE double minkowski(double x, double mu, double sigma, double k);

  /** Calculates the power of the absolute-value of some number and re-applies its original sign afterwards, if base==0, the result will be 
  0 also. */
  INLINE double powBipolar(double base, double exponent);

  /** Power function for intergers. Multiplies the base exponent-1 time with itself and returns the result. */
  //INLINE int powInt(int base, int exponent);

  /** Generates a 2*pi periodic pulse wave with a duty-cycle (high level fraction) of d. The DC amount of -2*d is compensated such that the 
  resulting waveform is DC-free. */
  INLINE double pulseWave(double x, double d);

  /** Quantizes x to the interval i. */
  INLINE double quant(double x, double i);

  /** Quantizes x (assumed to be in -1.0 - 1.0) to a grid which represents a resolution of nBits bits. */
  INLINE double quantToBits(double x, double nBits);

  /** Returns the absolute value of x. */
  INLINE double rabs(double x);

  /** Generates a pseudo-random number between min and max. */
  INLINE double random(double min=0.0, double max=1.0);

  /** Generates a random number between -1.0...+1.0 in double format. */
  INLINE double randDbl();

  /** Calculates the sine of x assuming that x is in the range 0...pi/2 - this is not checked for. */
  INLINE double rsSin(double x, int numTerms = 16);

  /** Generates a 2*pi periodic saw wave. */
  INLINE double sawWave(double x);

  /** Calculates sine and cosine of x - this is more efficient than calling sin(x) and cos(x) seperately. */
  INLINE void sinCos(double x, double* sinResult, double* cosResult);

  /** Calculates hyperbloic sine and cosine of x -this is more efficient than calling sinh(x) and cosh(x) seperately. */
  INLINE void sinhCosh(double x, double* sinhResult, double* coshResult);

  /** Calculates a parabolic approximation of the sine and cosine of x. */
  INLINE void sinCosApprox(double x, double* sinResult, double* cosResult);

  /** Applies a softclipping with a tanh shaped transtion region between the linear part and the clipping value. */
  INLINE double softClip(double x, double lowClamp   = -1.0, double highClamp  =  1.0,
                                   double lowKnee    =  0.1, double highKnee   =  0.1 );

  /** Generates a 2*pi periodic square wave. */
  INLINE double sqrWave(double x);

  /** Same syntax as fmod but brings the result into the symmetric interval -y/2...y/2 */
  INLINE double srem(double x, double y);

  /** The unit step function - returns 1 for x>=0, 0 else. */
  INLINE double step(double x);

  /** Rational approximation of the hyperbolic tangent. */
  INLINE double tanhApprox(double x);

  /** Generates a 2*pi periodic triangle wave. */
  INLINE double triWave(double x);

  /** Periodically wraps a number into an interval between min and max, for example, an arbitrary angle x may be wrapped into the interval 
  0...2*PI via xWrapped = wrapToInterval(x, 0.0, 2*PI) or into the interval -PI...PI via xWrapped = wrapToInterval(x, -PI, PI). */
  INLINE double wrapToInterval(double x, double min, double max);

  /** Just outputs the constant value 0.0 for all inputs - used as default function pointer when client code selects an invalid 
  function-index, for example in the waveform-renderers. */
  INLINE double zeroFunction(double x);

  // wrapAround

  //===============================================================================================
  //implementation:

  INLINE double acosh(double x)
  {
    if( x >= 1 )
      return log(x + sqrt(x*x-1) );
    else
    {
      DEBUG_BREAK; // argument must be >= 1, for the general case use the complex version acosh_c defined in rosic_ComplexFunctions.h
      return 0.0;
    }
  }

  INLINE double acoshs(double x)
  {
    if( fabs(x) >= 1 )
      return log(fabs(x + sqrt(x*x-1)) );
    else
    {
      DEBUG_BREAK; // fabs(x) must be >= 1, for the general case use the complex version acosh_c defined in rosic_ComplexFunctions.h
      return 0.0;
    }
  }

  // not needed anymore - provided in math.h
  //INLINE double asinh(double x)
  //{
  //  return log(x + sqrt(x*x+1) );
  //}

  INLINE double belowOrAbove(double x, double low, double high)
  {
    if( x < low )
      return -1.0;
    else if ( x > high )
      return 1.0;
    else
      return 0.0;
  }

  //INLINE double besselj(double x, double order) { return jv(order, x); }

  INLINE double binomialCoefficientReal(double n, int k)
  {
    if( k == 0 )
      return 1.0;
    else
    {
      double num = n;
      for(int i=1; i<k; i++)
        num *= (n - (double)i);
      return num / (double) factorial(k);
    }
  }

  INLINE double cheby(double x, double n)
  {
    if( fabs(x) <= 1.0 )
      return cos( n * acos(x) );
    else if( x > 0.0 )
      return cosh( n * acosh(x) );
    else
      return pow(-1.0, roundToInt(n) ) * cosh( n * acosh(-x) );
  }

  INLINE double cheby(double x, int n)
  {
    double t0 = 1.0;
    double t1 = x;
    double tn = 1.0;
    for(int i=0; i<n; i++)
    {
      tn = 2*x*t1 - t0;
      t0 = t1;
      t1 = tn;
    }
    return t0;
  }

  /*
  INLINE double ceil(double x)
  {
    //double frac = x - floor(x);
    return x + (1-frac(x));
  }
  */

  template <class T>
  INLINE T clip(T x, T min, T max)
  {
    if( x > max )
      return max;
    else if ( x < min )
      return min;
    else return x;
  }

  INLINE double cosSquaredApprox(double x)
  {
    // approximant: y = a*(1-x^2)^2 + (1-a)*(1-x^2)^3;

    double a  = 0.564765625;     // ...minimax-optimized
    double b  = 1.0 - a;
    double y1 = 1.0 - x*x;  // y1 = (1-x^2)
    double y2 = y1*y1;      // y2 = (1-x^2)^2
           y1 = y1*y2;      // y1 = (1-x^2)^3
    return a*y2 + b*y1;
  }

  INLINE void equalPowerGainFactors(double x, double *gain1, double *gain2,
    double xMin, double xMax)
  {
    double tmp = linToLin(x, xMin, xMax, 0.0, PI/2.0);
    sinCos(tmp, gain2, gain1);
  }

  INLINE double evaluateQuartic(double x, double a0, double a1, double a2, double a3, double a4)
  {
    double x2 = x*x;
    return x*(a3*x2+a1) + x2*(a4*x2+a2) + a0;
  }

  INLINE double foldOver(double x, double min, double max)
  {
    if( x > max )
      return max - (x-max);
    else if( x < min )
      return min - (x-min);
    else return x;
  }

  INLINE double frac(double x)
  {
    return x - floor(x);
  }

  INLINE double gauss(double x, double mu, double sigma)
  {
    return (1.0 / (sqrt(2*PI)*sigma)) * exp( - ((x-mu)*(x-mu)) / (2*sigma*sigma) );
  }

  INLINE int intAndFracPart(double x, double &frac)
  {
    int i = (int) x;
    frac  = x - (double) i;
    return i;
    // this version takes roughly 40 cycles

    /*
    double f = floor(x);
    frac     = x - f;
    return roundToInt(f);
    // this version takes roughly 50 cycles
    */
  }

  INLINE double integerPower(double x, int exponent)
  {
    double accu = 1.0;
    for(int i=0; i<exponent; i++)
      accu *= x;
    return accu;
  }

  template <class T>
  INLINE T limit(T x, T lowerBound, T upperBound)
  {
    if( x < lowerBound )
      return lowerBound;
    if( x > upperBound )
      return upperBound;
    return x;
  }

  /*
  INLINE double rosic::log10(double x)
  {
    return ONE_OVER_LN10*log(x);
  }
  */

  INLINE double logistic(double x, double b)
  {
    return 1.0 / ( 1 + exp(-b*x) );
  }

  INLINE void makeBlackmanWindow(double *window, int length)
  {
    for(int n=0; n<length; n++)
      window[n] = 0.42 - 0.5*cos( 2.0*PI*n / (double) (length-1))
                       + 0.08*cos(4.0*PI*n / (double) (length-1)) ;
  }

  INLINE void makeCosinePowerWindow(double *window, int length, double power)
  {
    for(int n=0; n<length; n++)
      window[n] = pow( sin(PI * (double) n / (double) length), power );
  }

  INLINE void makeHammingWindow(double *window, int length)
  {
    for(int n=0; n<length; n++)
      window[n] = 0.54 -  0.46 * cos(2.0*PI*n / (double) (length-1)) ;
  }

  INLINE void makeHanningWindow(double *window, int length)
  {
    for(int n=0; n<length; n++)
      window[n] = 0.5 * ( 1.0 - cos(2.0*PI*n / (double) (length-1)) );
  }

  INLINE void makeRectangularWindow(double *window, int length)
  {
    for(int n=0; n<length; n++)
      window[n] = 1.0;
  }

  INLINE double minkowski(double x, double mu, double sigma, double k)
  {
    return (1 /(sqrt(2*PI)*sigma)) * exp( - ( pow(fabs(x-mu),k) ) / (2*pow(sigma,k)) );
    //return (1 /(sqrt(2*PI*pow(sigma,k)))) * exp( - ( pow(abs(x-mu),k) ) / (2*pow(sigma,k)) );
  }

  INLINE double powBipolar(double base, double exponent)
  {
    if( base > 0.0 )
      return pow(base, exponent);
    else if( base < 0.0 )
      return sign(base) * pow(fabs(base), exponent);
    else
      return 0.0;
  }

  INLINE double pulseWave(double x, double d)
  {
    double tmp = fmod(x, 2*PI);
    if( tmp < PI*d || tmp >= 2*PI-PI*d )
      return 1.0 + 2*d;
    else
      return -1.0 + 2*d;
  }

  /*
  INLINE int powInt(int base, int exponent)
  {
    int result = base;
    for(int p=1; p<exponent; p++)
      result *= base;
    return result;
  }
  */

  INLINE double quant(double x, double i)
  {
    if( i <= DBL_MIN )
    {
      DEBUG_BREAK;
      return x;
    }

    double tmp = fmod(x, i);

    if( x >= 0 ) 
    {
      if( tmp < 0.5*i )
        return x - tmp;
      else
        return x - tmp + i;
    }
    else
    {
      if( tmp < -0.5*i )
        return x - tmp - i;
      else
        return x - tmp;
    }
  }

  INLINE double quantToBits(double x, double nBits)
  {
    double interval = 2 / pow(2.0, nBits);
    return quant(x, interval);
  }

  INLINE double rabs(double x)  // actually slower than fabs - grrrr
  {
    unsigned long long *pInt = reinterpret_cast<unsigned long long *> (&x);
    *pInt &= 0x7FFFFFFFFFFFFFFFULL;
    return *reinterpret_cast<double *> (pInt);
  }

  INLINE double randDbl()
  {
    return (  ( (double)rand() - 16383.5 ) / 16383.5 );
  }

  INLINE double random(double min, double max)
  {
    double tmp = (1.0/RAND_MAX) * rand() ;  // between 0...1
    return linToLin(tmp, 0.0, 1.0, min, max);
  }

  /*
  INLINE int rosic::roundToInt(double const x)
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
  */

  /*
  INLINE int rosic::roundToIntSSE2(double const x)
  {
    return _mm_cvtsd_si32(_mm_load_sd(&x));
  }
  */

  INLINE void sinCos(double x, double* sinResult, double* cosResult)
  {
    double s, c;     // do we need these intermediate variables?
    #ifdef LINUX     // assembly-version causes compiler errors on linux (gcc-thing?)
      *sinResult = sin(x);
      *cosResult = cos(x);
    #else
      __asm fld x
      __asm fsincos
      __asm fstp c
      __asm fstp s
      *sinResult = s;
      *cosResult = c;
    #endif
  }

  INLINE void sinhCosh(double x, double* sinhResult, double* coshResult)
  {
    double ep = exp( x);
    double em = exp(-x);
    *sinhResult = 0.5*(ep-em);
    *coshResult = 0.5*(ep+em);
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

  INLINE double rsSin(double x, int numTerms)
  {
    /*
    double accu = 1.0;
    double n    = 1.0;
    double nSq  = 1.0;
    double xSq  = x*x;
    double piSq = PI*PI;
    for(int i=1; i<=numTerms; i++)
    {
      accu *= (1.0 - xSq / (piSq*nSq));
      n    += 1.0;
      nSq   = n*n;
    }
    return x*accu;
    */

    long double accu = 0.0;
    long double num  = x;    
    long double xSq  = x*x;  
    for(int k=0; k<numTerms; k+=2)
    {
      accu += inverseFactorials[2*k+1] * num;
      num  *= xSq;
      accu -= inverseFactorials[2*k+3] * num;
      num  *= xSq;
    }
    return (double) accu;

    /*
    double xSq   = x*x;
    double num1  = x;      // first term is x^1
    double accu1 = 0.0;    // accumulates the positive terms
    double num2  = x*x*x;  // first term is x^3
    double accu2 = 0.0;    // accumulates the negative terms
    double x4    = num2*x; // x^4
    for(int k=0; k<numTerms; k+=2)
    {
      accu1 += inverseFactorials[2*k+1] * num1;
      num1  *= xSq;
      accu2 -= inverseFactorials[2*k+3] * num2;
      num2  *= xSq;
    }
    return accu1 - accu2;
    */
  }

  INLINE double sawWave(double x)
  {
    double tmp = fmod(x, 2*PI);
    if( tmp < PI )
      return tmp/PI;
    else
      return (tmp/PI)-2.0;
  }

  INLINE double sqrWave(double x)
  {
    double tmp = fmod(x, 2*PI);
    if( tmp < PI )
      return 1.0;
    else
      return -1.0;
  }

  INLINE double softClip(double x,  double lowClamp, double highClamp,
                                    double lowKnee,  double highKnee)
  {
    double mid        = 0.5*(highClamp+lowClamp);
    double high       = highClamp-mid;
    double low        = lowClamp-mid;
    double highThresh = mid + (1.0-highKnee)*high;
    double lowThresh  = mid + (1.0-lowKnee)*low;

    // we want to deal with the absolute values of the thresholds and clamp values:
    lowThresh = mid-lowThresh;
    lowClamp  = mid-lowClamp;

    double highAlpha = highClamp - highThresh;
    if(highAlpha == 0.0)
      highAlpha += DBL_MIN;
    double highBeta = 1.0/highAlpha;

    double lowAlpha = lowClamp - lowThresh;
    if(lowAlpha == 0.0)
      lowAlpha += DBL_MIN;
    double lowBeta = 1.0/lowAlpha;

    double out;
    if(x >= mid)     // positive half-wave
    {
      if(x <= highThresh)
        out = x;
      else
        out = highThresh + highAlpha*tanh(highBeta*(x-highThresh));
    }
    else              // negative half-wave
    {
      x = mid-x;
      if(x <= lowThresh)
        out = x;
      else
        out = lowThresh + lowAlpha*tanh(lowBeta*(x-lowThresh));
      out = mid-out;
    }

    return out;
  }

  INLINE double srem(double x, double y)
  {
    double z = fmod(x, y);

    if( fabs(z) > 0.5*y )
      z = z - y*sign(z);

    return z;
  }

  INLINE double step(double x)
  {
    if( x >= 0.0 )
      return 1.0;
    else
      return 0.0;
  }

  INLINE double tanhApprox(double x)
  {
    double a = fabs(2*x);
    double b = 24+a*(12+a*(6+a));
    return 2*(x*b)/(a*b+48);
  }

  INLINE double triWave(double x)
  {
    double tmp = fmod(x, 2*PI);
    if( tmp < 0.5*PI )
      return tmp/(0.5*PI);
    else if( tmp < 1.5*PI )
      return 1.0 - ((tmp-0.5*PI)/(0.5*PI));
    else
      return -1.0 + ((tmp-1.5*PI)/(0.5*PI));
  }

  INLINE double wrapToInterval(double x, double min, double max)
  {
    double r   = max-min;   // range
    double tmp = x-min;
    if( tmp >= 0.0 )
      tmp = fmod(tmp, r);
    else
    {
      tmp = fmod(fabs(tmp), r);
      tmp = r-tmp;
    }
    return tmp + min;
  }

  INLINE double zeroFunction(double x)
  {
    return x = 0.0;
  }

} // end namespace rosic

#endif // #ifndef rosic_ElementaryFunctionsReal_h
