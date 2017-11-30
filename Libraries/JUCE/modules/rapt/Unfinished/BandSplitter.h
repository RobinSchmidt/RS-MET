#ifndef RAPT_BANDSPLITTER_H_INCLUDED
#define RAPT_BANDSPLITTER_H_INCLUDED

/** A filter pair to split an incoming signal into lowpass- and a highpass part. The highpass part
is obtained by subtracting the lowpass part from the original signal. This means that adding the
highpass and lowpass outputs together gives the original signal back - the splitter provides
perfect reconstruction. That's a feature that not all band-splitters have. For example, 
Linkwitz/Riley splitters provide only allpass reconstruction (adding high and low bands gives an
allpassed version of the input). */

template<class TSig, class TPar>
class rsTwoBandSplitter
{

public:

  /** Sets the normalized radian frequency at which the split occurs. */
  void setOmega(TPar newOmega);

  //void setOrder(int newOrder); // to be added later

  /** Returns the normalized radian frequency at which the split occurs. */
  TPar getOmega() const { return w; }

  /** Resets state buffer variables */
  void reset() { x1 = y1 = 0; }


  inline TSig getLowpassSample(TSig in)
  {
    y1 = b0*in + b1*x1 - a1*y1;
    x1 = in;
    return y1;
  }

  /** Produces a pair of a lowpass and a highpass output sample from an input sample. */
  inline void getSamplePair(TSig in, TSig* lo, TSig* hi)
  {
    *lo = getLowpassSample(in);
    *hi = in - *lo; // ensures tha in = lo + hi (perfect reconstruction)
  }

protected:

  TSig x1 = 0, y1 = 0; // buffers
  TPar w  = 0;         // normalized radian frequency ("omega")
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

  enum slopeAccumulationModes
  {
    ACCUMULATE_INTO_HIGHPASS = 0,
    ACCUMULATE_INTO_LOWPASS
  };

  /** Destructor. Clears array of bandsplitter objects. */
  ~rsMultiBandSplitter();

  /** Sets up the sample rate. */
  void setSampleRate(TPar newSampleRate);

  /** Sets the splitting frequency for the band with given index. */
  void setSplitFrequency(int bandIndex, TPar newFrequency);

  /** Sets up all the splitting frequencies at once. */
  void setSplitFrequencies(const std::vector<TPar>& newFrequencies);

  //void setNumBands(int newNumBands);

  /** Adds a new band with the given splitting frequency. */
  void addBand(TPar splitFrequency);


  void setSlopeAccumulationMode(int newMode) { mode = newMode; }






  int getNumBands() { return (int)splitters.size() + 1; }


  /** Produces one output sample frame. The frequency bands are in ascending order (from lowpass 
  through the variosu bandpasses up to highpass). The caller must make sure that the output array 
  is at least as long as the number of bands. */
  void processSampleFrame(TSig in, TSig* outs)
  {
    TSig lo, hi;  // temporaries
    size_t N = splitters.size();
    switch(mode)
    {
    case ACCUMULATE_INTO_HIGHPASS: {   // slope accumulates into highpass band
      hi = in;
      for(size_t k = 0; k < N; k++) {
        splitters[k]->getSamplePair(hi, &lo, &hi);
        outs[k] = lo; }
      outs[N] = hi; // is that right?
    } break;
    case ACCUMULATE_INTO_LOWPASS: {   // slope accumulates into lowpass band
      lo = in;
      for(size_t k = 0; k < N; k++) {
        splitters[N-1-k]->getSamplePair(lo, &lo, &hi);
        outs[N-1-k] = hi; }
      outs[0] = lo; // ?
    } break;
    }
  }

protected:

  /** Clears or arrays of band-splitter objects and frequencies. */
  void clearArrays();

  /** Updates our array of two-way splitters. */
  void updateSplitters();

  std::vector<TPar> splitFreqs;  // splitting frequencies
  std::vector<rsTwoBandSplitter<TSig, TPar>*> splitters;

  TPar sampleRate = 44100;

  int mode = ACCUMULATE_INTO_HIGHPASS;

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