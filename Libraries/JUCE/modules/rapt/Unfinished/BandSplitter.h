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

  enum slopeAccumulationModes // rename to splitmodes
  {
    ACCUMULATE_INTO_HIGHPASS = 0,  // always splits the high band further
    ACCUMULATE_INTO_LOWPASS,       // always splits the low band further
    BINARY_TREE                    // builds a binary tree of splitters (assumes numBands = 2^k)
  };

  /** Destructor. Clears array of bandsplitter objects. */
  ~rsMultiBandSplitter();

  /** Sets up the sample rate. */
  void setSampleRate(TPar newSampleRate);

  /** Sets the splitting frequency for the band with given index. */
  void setSplitFrequency(int bandIndex, TPar newFrequency);

  /** Sets up all the splitting frequencies at once. */
  void setSplitFrequencies(const std::vector<TPar>& newFrequencies);

  /** Sets the number of bands into which the signal will be split. Bands are ordered from the 
  lowest to the highest, so when the new number is less than the old number, bands will be taken 
  away from the high side. When the new number is greater than the old ones, the newly added bands
  will their split frequencies evenly distributed on a log-frequency scale between the old topmost
  split-frequency and the Nyquist frequency. */
  void setNumberOfBands(int newNumber);
    // rename to setNumberOfAvailableBands

  /** Sometimes it may be convenient (or more efficient) to not change the arrays of splitters 
  (causing deletions or allocations of new objects) but just use a lesser number of splitters than 
  available. Use this function to set the number of af active bands to less than or equal the 
  number of available bands. */
  void setNumberOfActiveBands(int newNumber);

  /** Adds a new band with the given splitting frequency. */
  void addBand(TPar splitFrequency);


  void setSlopeAccumulationMode(int newMode) { mode = newMode; }

  int getNumBands() { return (int)splitters.size() + 1; }

  /** Resets the states of the band-splitter filters. */
  void reset();


  /** Produces one output sample frame. The frequency bands are in ascending order (from lowpass 
  through the various bandpasses up to highpass). The caller must make sure that the output array 
  is at least as long as the number of bands. */
  void processSampleFrame(TSig in, TSig* outs)
  {
    //int numSplitters = (int) splitters.size(); // use getNumActiveBands
    int numSplitters = numActiveBands - 1;
    switch(mode)
    {
    case ACCUMULATE_INTO_HIGHPASS: {   // slope accumulates into highpass band
      for(int k = 0; k < numSplitters; k++)
        splitters[k]->getSamplePair(in, &outs[k], &in);
      outs[numSplitters] = in; // "outs" has to be one slot longer than splitters array
    } break;
    case ACCUMULATE_INTO_LOWPASS: {   // slope accumulates into lowpass band
      for(int k = numSplitters-1; k >= 0; k--) // here, k indeed needs to be a signed integer
        splitters[k]->getSamplePair(in, &in, &outs[k+1]);
      outs[0] = in;
    } break;
    case BINARY_TREE:
    {
      int numBands = getNumBands(); // currently assumes numBands to be a power of 2
      int inc = numBands;
      outs[0] = in;
      while(inc > 1){
        int pos = 0;
        while(pos < numBands){
          //int splitterIndex = pos+inc/2-1; // for debug
          splitters[pos+inc/2-1]->getSamplePair(outs[pos], &outs[pos], &outs[pos+inc/2]);
          pos += inc; }
        inc /= 2; }

      // maybe here we should collect/recombine the topmost bands, if the number of active bands is
      // less than the number of available bands - available bands should always be a power of 2

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
  int numActiveBands = 1;

};

#endif