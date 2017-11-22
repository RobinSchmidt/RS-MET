#ifndef RAPT_ELLIPTICSUBBANDFILTERDIRECTFORM_H_INCLUDED
#define RAPT_ELLIPTICSUBBANDFILTERDIRECTFORM_H_INCLUDED

/** This is an elliptic subband filter of 12th order using a Direct Form II implementation 
structure. */

template<class TSig, class TPar>
class rsEllipticSubBandFilterDirectForm
{

public:

  /** Constructor. Initializes coefficients for an elliptic halfband filter. */
  rsEllipticSubBandFilterDirectForm();

  /** Sets the subdivision factor, for example 2 for a halfband filter (which passes everything
  below half the Nyquist frequency and stops everything above) or 4 for a quarterband filter. */
  void setSubDivision(TPar newSubDivision);

  /** Resets the filter state. */
  void reset();

  /** Calculates a single filtered output-sample. */
  inline TSig getSample(TSig in)
  {
    /*
    // this is the straightforward, non-optimized version:
    double tmp =   in + TINY
                 - a[1]*w[0] - a[2]*w[1] - a[3]*w[2] - a[4]*w[3]  - a[5]*w[4]   - a[6]*w[5]
                 - a[7]*w[6] - a[8]*w[7] - a[9]*w[8] - a[10]*w[9] - a[11]*w[10] - a[12]*w[11];

    double y =     b[0]*tmp
                 + b[1]*w[0] + b[2]*w[1] + b[3]*w[2] + b[4]*w[3]  + b[5]*w[4]   + b[6]*w[5]
                 + b[7]*w[6] + b[8]*w[7] + b[9]*w[8] + b[10]*w[9] + b[11]*w[10] + b[12]*w[11];
    */

    // calculate intermediate and output sample via direct form II - the parentheses facilitate 
    // out-of-order execution of the independent additions (for performance optimization):
    TSig tmp = in
      - ((a[1]*w[0] + a[2]*w[1]) + (a[3]*w[2]   + a[4]*w[3]))
      - ((a[5]*w[4] + a[6]*w[5]) + (a[7]*w[6]   + a[8]*w[7]))
      - ((a[9]*w[8] + a[10]*w[9]) + (a[11]*w[10] + a[12]*w[11]));

    TSig y = b[0]*tmp
      + ((b[1]*w[0] + b[2]*w[1])  +  (b[3]*w[2]   + b[4]*w[3]))
      + ((b[5]*w[4] + b[6]*w[5])  +  (b[7]*w[6]   + b[8]*w[7]))
      + ((b[9]*w[8] + b[10]*w[9]) +  (b[11]*w[10] + b[12]*w[11]));

    // update state variables:
    memmove(&w[1], &w[0], 11*sizeof(TSig));
    w[0] = tmp;

    return y;
  }

protected:

  TSig w[12];   // state buffer
  TPar a[13];   // feedback coefficients
  TPar b[13];   // feedforward coefficients

};

#endif
