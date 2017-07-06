#ifndef RAPT_SINCOSTABLE_H_INCLUDED
#define RAPT_SINCOSTABLE_H_INCLUDED

/** A class for getting approximate values of the sine and cosine of a given angle using lookup 
tables. It's meant for optimizing sin/cos calculations. */

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
    //int i  = (int)floor(pos); // this would work equally well for negative inputs
    T frac = pos-i;
    T wi   = 1-frac;
    i      =  i    & mask;
    int i1 = (i+1) & mask;
    *sinValue = wi * sinTbl[i] + frac * sinTbl[i1];
    *cosValue = wi * cosTbl[i] + frac * cosTbl[i1];
  }




protected:

  T scaler;  
  int mask;
  std::vector<T> sinTbl, cosTbl;

};

#endif
