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

  void setNumberOfBands(int newNumber);

  void setSplitFrequency(int bandIndex, double newFrequency);

  void setThreshold(int bandIndex, double newThreshold);

  void setRatio(int bandIndex, double newRatio);

  void setAttackTime(int bandIndex, double newAttackTime);

  void setReleaseTime(int bandIndex, double newReleaseTime);

  //-----------------------------------------------------------------------------------------------
  /** \name Processing */

  /** Calculates one output stereo sample-frame at a time. */
  INLINE void getSampleFrameStereo(double *inOutL, double *inOutR);

  //===============================================================================================

protected:

  rapt::rsMultiBandSplitter<double, double> splitterL, splitterR;
  std::vector<Compressor*> compressors;
  std::vector<double> tmpL, tmpR; // temporary buffers
  int numBands = 1;
};

//-----------------------------------------------------------------------------------------------
// inlined functions:

INLINE void rsMultiBandCompressor::getSampleFrameStereo(double *inOutL, double *inOutR)
{
  splitterL.processSampleFrame(*inOutL, &tmpL[0]);
  splitterR.processSampleFrame(*inOutR, &tmpR[0]);


}

}

#endif 
