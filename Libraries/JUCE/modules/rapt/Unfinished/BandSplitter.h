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


  inline TSig getLowpassSample(TSig in)
  {
    y1 = b0*in + b1*x1 - a1*y1;
    x1 = in;
    return y1;
  }

  inline void getSamplePair(TSig in, TSig* lo, TSig* hi)
  {
    *lo = getLowpassSample(in);
    *hi = in - *lo;
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



  /** Sets up the sample rate. */
  void setSampleRate(TPar newSampleRate);

  /** Sets up all the splitting frequencies at once. */
  void setSplitFrequencies(const std::vector<TPar>& newFrequencies);

  void setNumBands(int newNumBands);

  void setSplitFrequency(int bandIndex, TPar newFrequency);

  /** Produces one output sample frame. The frequency bands are in ascending order and the called 
  must make sure that the output array is at least as long as the number of bands. */
  void processSampleFrame(TSig in, TSig* outs)
  {
    TSig lo, hi;  // temporaries
    size_t N = splitters.size();

    // slope accumulates into lowpass band:
    lo = in;
    for(size_t k = 0; k < N; k++) {
      splitters[N-1-k]->getSamplePair(lo, &lo, &hi);
      outs[N-1-k] = hi; }

    // slope accumualtes into highpass band:





  }

protected:

  /** Updates our array of two-way splitters. */
  //void updateSplitters();

  std::vector<rsTwoBandSplitter<TSig, TPar>*> splitters;

  int mode = 0;

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