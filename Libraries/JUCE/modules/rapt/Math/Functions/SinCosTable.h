#ifndef RAPT_SINCOSTABLE_H_INCLUDED
#define RAPT_SINCOSTABLE_H_INCLUDED

/** A class for getting approximate values of the sine and cosine of a given angle using linearly 
interpolating lookup tables. It's meant for optimizing sin/cos calculations. */

template<class T>
class rsSinCosTable
{

public:

  /** Constructor. The tableSize should be a power of two. */
  rsSinCosTable(int tableSize = 1024);


  /** Produces a value of the sine and cosine for the given input value x. */
  inline void getSineAndCosine(T x, T* sinValue, T* cosValue)
  {
    //// preliminary - use tables later:
    //*sinValue = sin(x);
    //*cosValue = cos(x);

    T pos = scaler * x;  // continuous readout index
    //int i = ((int)(pos)) & mask; // truncated - equally good/bad in positiv/negative range
    int i = ((int)(pos+0.5)) & mask; // rounded - better for positive, worse for negative x

    *sinValue = sinTbl[i];
    *cosValue = cosTbl[i];
  }

protected:

  T scaler;  
  int mask;
  std::vector<T> sinTbl, cosTbl;

};

#endif
