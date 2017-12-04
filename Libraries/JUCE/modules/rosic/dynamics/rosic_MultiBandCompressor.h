#ifndef rosic_MultiBandCompressor_h
#define rosic_MultiBandCompressor_h

namespace rosic
{

/**  */

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

  void setNumerOfBands(int newNumber);

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

};

//-----------------------------------------------------------------------------------------------
// inlined functions:



INLINE void rsMultiBandCompressor::getSampleFrameStereo(double *inOutL, double *inOutR)
{

}

}

#endif 
