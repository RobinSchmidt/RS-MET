#ifndef RAPT_BANDSPLITTER_H_INCLUDED
#define RAPT_BANDSPLITTER_H_INCLUDED

/** A filter pair to split an incoming signal into lowpass- and a highpass part. */

template<class TSig, class TPar>
class rsTwoBandSplitter
{

public:

  /** Sets the normalized radian frequency at which the split occurs. */
  void setOmega(TPar newOmega);

  /** Resets state buffer variables */
  void reset() { x1 = y1 = 0; }

  inline void getSamplePair(TSig in, TSig* lo, TSig* hi)
  {
    // state update:
    y1 = b0*in + b1*x1 - a1*y1;
    x1 = in;

    // assign outputs:
    *lo = y1;
    *hi = in - y1;

    // maybe factor out a function getLowpassSample
  }

protected:

  TSig x1 = 0, y1 = 0; // buffers
  TPar b0 = 0, b1 = 0; // feedforward coeffs
  TPar a1 = 0;         // feedback coeff

};

//=================================================================================================

/** Uses an arbitrary number of two-band splitters to split a signal into an arbitrary number of 
bands. */

template<class TSig, class TPar>
class rsMultiBandSplitter
{

public:

protected:

  std::vector<rsTwoBandSplitter<TSig, TPar>*> splitters;

};

/*
-split the signal in a hierarchical way
 -always split the low band further
 -..or maybe always split the high band further - this makes the final high pass steepest and i 
  think it's more important to keep the high band free from low-frequency leakage than vice versa
 -or maybe split in a binary tree like way (distributes the slopes evenly between high and low 
  bands)
 -maybe make this a user choice
*/


#endif