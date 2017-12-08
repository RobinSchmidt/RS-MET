#ifndef rosic_MultiBandCompressor_h
#define rosic_MultiBandCompressor_h

namespace rosic
{

/** A multiband compressor with an arbitrary number of bands and perfect reconstruction, i.e. with 
neutral settings, it will be totally transparent (up to roundoff error) and not color the signal in
any way. */

class rsMultiBandCompressor
{

public:

  //---------------------------------------------------------------------------------------------
/** \name Construction/Destruction */

  /** Constructor - constructs a dynamics processor with a given maximum number of samples
  lookahead. */
  rsMultiBandCompressor();

  /** Destructor */
  ~rsMultiBandCompressor();

  //---------------------------------------------------------------------------------------------
  /** \name Setup */

  void setSampleRate(double newSampleRate);

  void setNumberOfBands(int newNumber);

  void setSplitMode(int newMode);

  void setSplitFrequency(int bandIndex, double newFrequency);

  void setThreshold(int bandIndex, double newThreshold);

  void setRatio(int bandIndex, double newRatio);

  void setAttackTime(int bandIndex, double newAttackTime);

  void setReleaseTime(int bandIndex, double newReleaseTime);

  //---------------------------------------------------------------------------------------------
  /** \name Inquiry */

  /** Returns the maximum number of bands that is supported. */
  int getMaxNumberOfBands() { return maxNumBands; }

  //-----------------------------------------------------------------------------------------------
  /** \name Processing */

  /** Calculates one output stereo sample-frame at a time. */
  INLINE void getSampleFrameStereo(double *inOutL, double *inOutR);

  void reset();

  //===============================================================================================

protected:

  RAPT::rsMultiBandSplitter<double, double> splitterL, splitterR;
  std::vector<Compressor*> compressors;
  std::vector<double> tmpL, tmpR; // temporary buffers
  int numBands = 1;
  int maxNumBands = 16; // preliminary - make indefinite in the future
};

//-----------------------------------------------------------------------------------------------
// inlined functions:

INLINE void rsMultiBandCompressor::getSampleFrameStereo(double *inOutL, double *inOutR)
{
  // split inputs into frequency bands:
  splitterL.processSampleFrame(*inOutL, &tmpL[0]);
  splitterR.processSampleFrame(*inOutR, &tmpR[0]);

  // compress individual bands:
  for(int k = 0; k < numBands; k++)
    compressors[k]->getSampleFrameStereo(&tmpL[k], &tmpR[k]);

  // recombine compressed bands into output:
  *inOutL = *inOutR = 0;
  for(int k = 0; k < numBands; k++)
  {
    *inOutL += tmpL[k];
    *inOutR += tmpR[k];
  }
}

}

#endif 
