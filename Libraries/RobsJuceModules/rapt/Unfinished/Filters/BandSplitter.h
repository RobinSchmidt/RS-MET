#ifndef RAPT_BANDSPLITTER_H_INCLUDED
#define RAPT_BANDSPLITTER_H_INCLUDED

/** A filter pair to split an incoming signal into lowpass- and a highpass part. The highpass part
is obtained by subtracting the lowpass part from the original signal. This means that adding the
highpass and lowpass outputs together gives the original signal back - the splitter provides
perfect reconstruction. That's a feature that not all band-splitters have. For example, 
Linkwitz/Riley splitters provide only allpass reconstruction (adding high and low bands gives an
allpassed version of the input). 

ToDo:
-use an rsFirstOrderFilter object for the lowpass
-use a design formula with prescribed Nyquist gain (currently we use BLT) */

template<class TSig, class TPar>
class rsTwoBandSplitter
{
  typedef const TSig& CRSig;  // const reference to a signal value
  typedef const TPar& CRPar;  // const reference to a parameter value

public:

  /** Sets the normalized radian frequency at which the split occurs. */
  void setOmega(CRPar newOmega);

  //void setOrder(int newOrder); // to be added later

  /** Copies the settings (split frequency and coeffs) from another splitter. */
  void copySettingsFrom(const rsTwoBandSplitter& splitter);

  /** Copies the state (past in/out samples) from another splitter. */
  void copyStateFrom(const rsTwoBandSplitter& splitter);

  /** Returns the normalized radian frequency at which the split occurs. */
  TPar getOmega() const { return w; }


  /** Returns the z-domain transfer function value H(z) at the given z. The transfer function value
  of the highpass is then just 1-H(z). */
  std::complex<TPar> getLowpassTransferFunctionAt(const std::complex<TPar>& z) const;

  std::complex<TPar> getHighpassTransferFunctionAt(const std::complex<TPar>& z) const
  {
    return TPar(1) - getLowpassTransferFunctionAt(z);
  }





  /** Resets state buffer variables */
  void reset() { x1 = y1 = 0; }





  /** Produces one output sample at a time. */
  inline TSig getLowpassSample(CRSig in)
  {
    y1 = b0*in + b1*x1 - a1*y1;
    x1 = in;
    return y1;
  }


  /** Produces a pair of a lowpass and a highpass output sample from an input sample. */
  inline void getSamplePair(CRSig in, TSig* lo, TSig* hi)
  {
    TSig tmp = in;
    *lo = getLowpassSample(tmp);
    *hi = tmp - *lo; // ensures tha in = lo + hi (perfect reconstruction)

    // the test function bandSplittingMultiWay() fails, if we don't make a copy of "in" in "tmp"
    // ...i don't know why - perhaps something about overwriting something that shouldn't be?
    // well - the multiband splitter uses the same variable for input and one of the outputs...
    // -> figure out -> write better comment...yeah - i see - it's because we take the input signal
    // by (const) reference
  }

  // old:
  //inline void getSamplePair(TSig in, TSig* lo, TSig* hi)
  //{
  //  *lo = getLowpassSample(in);
  //  *hi = in - *lo; // ensures tha in = lo + hi (perfect reconstruction)
  //}

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
  typedef const TSig& CRSig;  // const reference to a signal value
  typedef const TPar& CRPar;  // const reference to a parameter value

public:

  enum slopeAccumulationModes // rename to splitmodes
  {
    ACCUMULATE_INTO_HIGHPASS = 0,  // always splits the high band further
    ACCUMULATE_INTO_LOWPASS,       // always splits the low band further
    BINARY_TREE                    // builds a binary tree of splitters (assumes numBands = 2^k)
  };

  /** Destructor. Clears array of bandsplitter objects. */
  ~rsMultiBandSplitter();

  //-----------------------------------------------------------------------------------------------
  // \name Setup

  /** Sets up the sample rate. */
  void setSampleRate(CRPar newSampleRate);

  /** Sets the splitting frequency for the band with given index. */
  void setSplitFrequency(int bandIndex, CRPar newFrequency);

  /** Sets up all the splitting frequencies at once. */
  void setSplitFrequencies(const std::vector<TPar>& newFrequencies);

  /** Sets the number of bands into which the signal will be split. Bands are ordered from the 
  lowest to the highest, so when the new number is less than the old number, bands will be taken 
  away from the high side. When the new number is greater than the old ones, the newly added bands
  will their split frequencies evenly distributed on a log-frequency scale between the old topmost
  split-frequency and the Nyquist frequency. */
  void setNumberOfBands(int newNumber);
  //void setNumberOfBands(int newNumber);
    // rename to setNumberOfAvailableBands

  /** Sometimes it may be convenient (or more efficient) to not change the arrays of splitters 
  (causing deletions or allocations of new objects) but just use a lesser number of splitters than 
  available. Use this function to set the number of af active bands to less than or equal the 
  number of available bands. */
  void setNumberOfActiveBands(int newNumber);

  /** Adds a new band with the given splitting frequency. */
  void addBand(CRPar splitFrequency);


  void insertBand(int index, CRPar splitFrequency); // not yet tested
  void removeBand(int index, bool mergeWithRightNeighbour = false); // not yet tested

  /** Copies the settings of the band with index src into the band with index dst. When copyState is 
  true, the filter's state variables will also be copied. This function is for facilitating the 
  insertion and removal of bands by an outside class. */
  void copyBandSettings(int src, int dst, bool copyState = true);



  //void setSlopeAccumulationMode(int newMode) { mode = newMode; }

  void setSplitMode(int newMode) { mode = newMode; }

  //-----------------------------------------------------------------------------------------------
  // \name Inquiry

  /** Returns the number of available bands, i.e. the maximum number of bands that is supported. */
  int getNumAvailableBands() const { return (int)splitters.size() + 1; }

  /** Returns the number of bands that are currently in use. */
  int getNumActiveBands() const { return numActiveBands; }


  /** Returns the upper cutoff frequency for the band with given index. */
  TPar getSplitFrequency(int index) const
  { 
    if(index < (int)splitFreqs.size())
      return splitFreqs[index];
    else
      return 0;
  }

  /** Returns the complex frequency response of the band with given index at the given 
  frequency. */
  std::complex<TPar> getBandFrequencyResponseAt(int bandIndex, CRPar frequency) const;

  /** Returns the magnitude response of the band with given index at the given frequency. */
  TPar getBandMagnitudeAt(int bandIndex, CRPar frequency) const
  {
    return abs(getBandFrequencyResponseAt(bandIndex, frequency));
  }

  /** Returns the smaple rate at which the splitter is currently running. */
  inline TPar getSampleRate() const { return sampleRate; }

  //-----------------------------------------------------------------------------------------------
  // \name Processing

  /** Produces one output sample frame. The frequency bands are in ascending order (from lowpass 
  through the various bandpasses up to highpass). The caller must make sure that the output array 
  is at least as long as the number of bands. */
  void processSampleFrame(CRSig input, TSig* outs)
  {
    //if(numActiveBands == 0)
    //  return;

    TSig in = input;


    //int numSplitters = (int) splitters.size(); // use getNumActiveBands
    int numSplitters = numActiveBands - 1; // assume numActiveBands >= 1
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


    // currently assumes numBands to be a power of 2 - not yet usable - maybe, we can keep using 
    // this code (without access violations) even, when splitters.size() is less than numBands as 
    // long as splitters.capacity() is larger?
    // maybe handle this by ensuring that the number of available bands is always a  power of two
    // the number of active bands can be less in which case the excess (available) bands are summed
    // into the topmost active band
    case BINARY_TREE:
    {
      int numBands = getNumAvailableBands();
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

  /** Resets the states of the band-splitter filters. */
  void reset();

protected:

  /** Clears or arrays of band-splitter objects and frequencies. */
  void clearArrays();

  /** Updates our array of two-way splitters. */
  void updateSplitters();

  std::vector<TPar> splitFreqs;  // splitting frequencies
  std::vector<rsTwoBandSplitter<TSig, TPar>*> splitters;

  TPar sampleRate = 44100;

  int mode = ACCUMULATE_INTO_HIGHPASS;
  int numActiveBands = 1; // make sure, that this is always right - or better: get rid of it

};

//=================================================================================================
/* Ideas

To make it more general in order to allow each 2-way splitter to be any kind of splitter (FIR, 
Linkwitz/Riley, perfect-reconstruction-IIR and maybe more) and also to allow each splitter to be
of a different kind than the others, maybe we could devise a suitable datastructure to encapsulate
the nesting of splitters together with compensation allpasses and delays.

maybe like this:

template<class TSig, class TPar>
class rsBandSplitterBase
{
public:

  virtual int getNumBands() const = 0;
  virtual void split(TSig in, TSig* bandOuts) const = 0;
  virtual void recombine(TSig* bandIns, TSig* sum) const = 0;

  // used for compensation delays and allpasses :

  // Returns the amount of delay that this splitter subclass produces.
  virtual int getDelay() const { return 0; }

  // Returns the order of the allpass that will be inflicted on the recombined signal.
  virtual int getAllpassOrder() const { return 0; }

  // Fills the given array with the allpass poles. Should be of length equal to the return value of 
  // getAllpassOrder()
  virtual void getAllpassPoles(std::complex<TPar>* poles) const {}

  // maybe it's too restrictive to require the compensation filter to be an allpass - allow any 
  // kind of compensation filter
};

class rsBandSplitterTwoWay : public rsBandSplitterBase
{
public:
  virtual int getNumBands() const override { return 2; }
};

class rsBandSplitCompensationFilter
{
// must implement an integer delayline and a chain of allpass filters
}

class rsBandSplitter : public rsBandSplitterBase
{
public:

  virtual int getNumBands() const override; // recursively calls getNumBands on nested splitters
  virtual void split(TSig in, TSig* bandOuts) const override;
  virtual void recombine(TSig* bandIns, TSig* sum) const override;

  virtual int getDelay() const override;
  virtual int getAllpassOrder() const override;
  virtual void getAllpassPoles(std::complex<TPar>* poles) const override;


  virtual void split(TSig in, TSig* bandOuts) const override
  {
    // apply twoWaySplitter to produce low and high band
    // pass these signals to splitterLow and splitterHigh
  }

  virtual void recombine(TSig* bandIns, TSig* sum) const override
  {
    // obtain recombinded low and high band from splitterLow, splitterHigh
    // apply compensation filters
    // add results and store in sum
  }

protected:

  rsBandSplitterTwoWay* twoWaySplitter;
  rsBandSplitter *splitterLow, *splitterHigh;
  rsBandSplitCompensationFilter compensationFilterLow, compensationFilterHigh;

};


// the actual splitting filter classes:
class rsBansSplitterTwoWayPerfect : public rsBandSplitterTwoWay
{
  // implements perfect reconstruction IIR band splitter
};

class rsBandSplitterTwoWayLinear : public rsBandSplitterTwoWay
{
  // implements linear phase FIR band splitter
public:
  virtual int getDelay() const override;
};

class rsBandSplitterTwoWayLinkwitzRiley : public rsBandSplitterTwoWay
{
  // implements Linkwitz/Riley band splitter
public:
  virtual int getAllpassOrder() const override;
  virtual void getAllpassPoles(std::complex<TPar>* poles) const override;
};

*/

#endif