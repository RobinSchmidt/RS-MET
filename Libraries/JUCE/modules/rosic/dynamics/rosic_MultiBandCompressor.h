#ifndef rosic_MultiBandCompressor_h
#define rosic_MultiBandCompressor_h

namespace rosic
{

/** Baseclass for multiband effects such as multiband-compressors, -distortion. etc. */

class rsMultiBandEffect
{

public:

  /** Constructor. */
  rsMultiBandEffect();

  //---------------------------------------------------------------------------------------------
  /** \name Setup */

  void setSampleRate(double newSampleRate);

  void setNumberOfBands(int newNumber);

  void setSplitMode(int newMode);

  void setSplitFrequency(int bandIndex, double newFrequency);

  //---------------------------------------------------------------------------------------------
  /** \name Inquiry */

  /** Retruns the number of bands that are currently active. */
  int getNumberOfBands() const { return numBands; }

  /** Returns the maximum number of bands that is supported. */
  int getMaxNumberOfBands() const { return maxNumBands; }

  /** Returns the upper cutoff frequency for the band with given index. */
  double getSplitFrequency(int index) { return splitterL.getSplitFrequency(index); }

  //---------------------------------------------------------------------------------------------
  /** \name Processing */

  /** Splits an incoming stereo signal into "numBands" bands and fills our member arrays 
  tmpL, tmpR with the individual band outputs. Subclasses are supposed to call this before 
  processing each band. */
  INLINE void split(double *inOutL, double *inOutR)
  {
    splitterL.processSampleFrame(*inOutL, &tmpL[0]);
    splitterR.processSampleFrame(*inOutR, &tmpR[0]);
  }

  /** Recombines the bands from our member arrays tmpL, tmpR into the outputs. Subclasses are 
  supposed to call this after processing each indiviudal band to produce the final recombined
  output. */
  INLINE void recombine(double *inOutL, double *inOutR)
  {
    *inOutL = *inOutR = 0;
    for(int k = 0; k < numBands; k++)
    {
      *inOutL += tmpL[k];
      *inOutR += tmpR[k];
    }
  }

  /** Resets the states of the band-splitting filters. */
  void reset();

protected:

  RAPT::rsMultiBandSplitter<double, double> splitterL, splitterR;
  std::vector<double> tmpL, tmpR; // temporary buffers
  int numBands = 1;
  int maxNumBands = 16; // preliminary - make indefinite in the future

};

//=================================================================================================

/** A multiband compressor with an arbitrary number of bands and perfect reconstruction, i.e. with 
neutral settings, it will be totally transparent (up to roundoff error) and not color the signal in
any way. */

class rsMultiBandCompressor : public rsMultiBandEffect
{

public:

  //---------------------------------------------------------------------------------------------
  /** \name Construction/Destruction */

  /** Constructor. */
  rsMultiBandCompressor();

  /** Destructor */
  ~rsMultiBandCompressor();

  //---------------------------------------------------------------------------------------------
  /** \name Setup */

  void setSampleRate(double newSampleRate);

  void setThreshold(int bandIndex, double newThreshold);

  void setRatio(int bandIndex, double newRatio);

  void setAttackTime(int bandIndex, double newAttackTime);

  void setReleaseTime(int bandIndex, double newReleaseTime);

  //-----------------------------------------------------------------------------------------------
  /** \name Processing */

  /** Calculates one output stereo sample-frame at a time. */
  INLINE void getSampleFrameStereo(double *inOutL, double *inOutR);

  void reset();

  //===============================================================================================

protected:

  std::vector<Compressor*> compressors;

};

//-----------------------------------------------------------------------------------------------
// inlined functions:

INLINE void rsMultiBandCompressor::getSampleFrameStereo(double *inOutL, double *inOutR)
{
  split(inOutL, inOutR);
  for(int k = 0; k < numBands; k++)  // compress individual bands
    compressors[k]->getSampleFrameStereo(&tmpL[k], &tmpR[k]);
  recombine(inOutL, inOutR);
}

}

#endif 
