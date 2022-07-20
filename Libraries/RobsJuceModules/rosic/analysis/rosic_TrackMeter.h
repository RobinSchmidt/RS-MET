#ifndef rosic_TrackMeter_h
#define rosic_TrackMeter_h

namespace rosic
{

/** This is a meaurement device which measures certain parameters of an incoming audio signal. It 
is mainly intended for visualization purposes.

ToDo: maybe switch all datatypes from double to float. Or move it to RAPT and templatize it.  */

class TrackMeter
{

public:

  //---------------------------------------------------------------------------------------------
  // construction/destruction:

  /** Constructor. */
  TrackMeter();

  /** Destructor. */
  ~TrackMeter();

  //---------------------------------------------------------------------------------------------
  // parameter settings:

  /** Sets the sample-rate. */
  void setSampleRate(double newSampleRate);

  /** Sets the attack time for the envelope followers */
  void setAttackTimeInMilliseconds(double newAttackTime);

  /** Sets the release time for the envelope followers */
  void setReleaseTimeInMilliseconds(double newReleaseTime);

  //---------------------------------------------------------------------------------------------
  // inquiry:

  /** Returns the current measurement value which is taken to be the maximum over all
  measurement values since the last call to this function. */
  SignalMeasures getCurrentMeasurement(bool reset = true);

  //---------------------------------------------------------------------------------------------
  // audio processing:

  /** Measures a stereo-ouput frame. */
  INLINE void measureSampleFrameStereo(double* inL, double* inR);

  //=============================================================================================

protected:

  double maxLeftLevel, maxRightLevel, maxMidLevel, maxSideLevel;
  double sumOfSquaresLeft, sumOfSquaresRight, sumOfProducts;
  int    sampleCounter;

  EnvelopeFollower leftLevelExtractor, rightLevelExtractor, midLevelExtractor,
    sideLevelExtractor;
  rsOnePoleFilterDD meanSquareExtractorLeft, meanSquareExtractorRight, productLevelExtractor;

  SignalMeasures currentMeasures;

};

//-----------------------------------------------------------------------------------------------
// from here: definitions of the functions to be inlined, i.e. all functions which are supposed 
// to be called at audio-rate (they can't be put into the .cpp file):

INLINE void TrackMeter::measureSampleFrameStereo(double* inL, double* inR)
{
  double tmp;

  tmp = leftLevelExtractor.getSample(*inL);
  if(tmp > maxLeftLevel)
    maxLeftLevel = tmp;

  tmp = rightLevelExtractor.getSample(*inR);
  if(tmp > maxRightLevel)
    maxRightLevel = tmp;

  tmp = midLevelExtractor.getSample(SQRT2_INV * (*inL + *inR));
  if(tmp > maxMidLevel)
    maxMidLevel = tmp;

  tmp = sideLevelExtractor.getSample(SQRT2_INV * (*inL - *inR));
  if(tmp > maxSideLevel)
    maxSideLevel = tmp;

  tmp                = productLevelExtractor.getSample((*inL) * (*inR));
  sumOfProducts     += tmp;
  tmp                = meanSquareExtractorLeft.getSample((*inL) * (*inL));
  sumOfSquaresLeft  += tmp;
  tmp                = meanSquareExtractorRight.getSample((*inR) * (*inR));
  sumOfSquaresRight += tmp;
  /*
  sumOfSquaresLeft  += (*inL) * (*inL);
  sumOfSquaresRight += (*inR) * (*inR);
  sumOfProduct      += (*inL) * (*inR);
  */

  sampleCounter++;
}

} // end namespace rosic

#endif // rosic_TrackMeter_h
