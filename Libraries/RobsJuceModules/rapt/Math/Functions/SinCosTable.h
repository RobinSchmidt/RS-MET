#ifndef RAPT_SINCOSTABLE_H_INCLUDED
#define RAPT_SINCOSTABLE_H_INCLUDED

/** A class for getting approximate values of the sine and cosine of a given angle using lookup 
tables. It's meant for optimizing sin/cos calculations. ..but you should really measure, whether
or not it gives any better performance than using the standard sin/cos functions. Often, it 
doesn't seem to be the case. 

maybe we should also try this:
http://lab.polygonal.de/2007/07/18/fast-and-accurate-sinecosine-approximation/

*/

template<class T>
class rsSinCosTable
{

public:

  /** Constructor. The tableSize should be a power of two. */
  rsSinCosTable(int tableSize = 1024);

  /** Produces a pair of sin/cos values using truncation of the continuous table index. It produces
  equal errors for the positive and negative range. */
  inline void getValuesTruncated(T x, T* sinValue, T* cosValue)
  {    
    int i = ((int)(scaler*x)) & mask; 
    *sinValue = sinTbl[i];
    *cosValue = cosTbl[i];
  }

  /** Produces a pair of sin/cos values using rounding of the continuous table index. It produces
  less error in the positive and more error in negative range compared to 
  getValuesTruncated. */
  inline void getValuesRounded(T x, T* sinValue, T* cosValue)
  {
    int i = ((int)(scaler*x+0.5)) & mask;
    *sinValue = sinTbl[i];
    *cosValue = cosTbl[i];
  }

  /** Produces a pair of sin/cos values using linear interpolation. Recommended only for 
  nonnegative inputs (for negative inputs, the error is large). */
  inline void getValuesLinear(T x, T* sinValue, T* cosValue)
  {
    T pos  = scaler * x;   // continuous readout index     
    int i  = (int)pos;
    //int i  = (int)floor(pos); // this would work equally well for negative inputs - weird: when
                                // using this, the performance test for the getValuesRounded is 
                                // affected
    T frac = pos-i;
    T wi   = 1-frac;
    i      =  i    & mask;
    int i1 = (i+1) & mask;
    *sinValue = wi * sinTbl[i] + frac * sinTbl[i1];
    *cosValue = wi * cosTbl[i] + frac * cosTbl[i1];
  }

  /** Produces a pair of sin/cos values using cubic interpolation. Recommended only for 
  nonnegative inputs (for negative inputs, the error is large). */
  inline void getValuesCubic(T x, T* sinValue, T* cosValue)
  {
    T pos  = scaler * x;     // continuous readout index     
    int i  = (int)pos;
    x      = pos-i;          // fractional part (input to polynomial)
    T x2   = x*x;
    T x3   = x*x2;
    i      =  i    & mask;
    int i1 = (i+1) & mask;

    // use the cosine-table as derivative table for the sine and the negative sine-table as 
    // derivative table for the cosine:
    T k0, k1;
    T y0, yp0, y1, yp1;

    y0  = sinTbl[i];
    y1  = sinTbl[i1];
    yp0 = cosTbl[i]  * scalerInv;
    yp1 = cosTbl[i1] * scalerInv;
    k0  = y1  - yp0 - y0;
    k1  = yp1 - yp0;
    *sinValue = y0 + yp0*x + (3*k0-k1)*x2 + (k1-2*k0)*x3;

    y0  =  cosTbl[i];
    y1  =  cosTbl[i1];
    yp0 = -sinTbl[i]  * scalerInv;
    yp1 = -sinTbl[i1] * scalerInv;
    k0  = y1  - yp0 - y0;
    k1  = yp1 - yp0;;
    *cosValue = y0 + yp0*x + (3*k0-k1)*x2 + (k1-2*k0)*x3;
  }



protected:

  T scaler, scalerInv;  
  int mask;
  std::vector<T> sinTbl, cosTbl;

};

#endif
