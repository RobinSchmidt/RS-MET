#ifndef rosic_NoiseGeneratorOld_h
#define rosic_NoiseGeneratorOld_h

//// standard-library indcludes:
//#include <stdlib.h>
//
//// rosic-indcludes:
////#include "../filters/rosic_IirFilter.h"
//#include "../filters/rosic_FourPoleFilter.h"


namespace rosic
{

 /**

 This is a noise generator capapble of producing various kinds of noise.

 */

 class NoiseGeneratorOld 
 {

 public:

  //---------------------------------------------------------------------------
  // construction/destruction:

  NoiseGeneratorOld();  ///< Constructor.
  ~NoiseGeneratorOld();  ///< Destructor.

  //---------------------------------------------------------------------------
  // parameter settings:

  void setSampleRate(double newSampleRate);
  ///< Overrides the setSampleRate() method of the AudioModule base class.

  //void setLowpassOrder(int newLowpassOrder);
  /**< sets the order of the lowpass-filter. */

  void setLowpassCutoff(double newLowpassCutoff);
  /**< sets the cutoff-frequency of the lowpass-filter. */

  //void setHighpassOrder(int newHighpassOrder);
  /**< sets the order of the highpass-filter. */

  void setHighpassCutoff(double newHighpassCutoff);
  /**< sets the cutoff-frequency of the highpass-filter. */

  //---------------------------------------------------------------------------
  // audio processing:

  INLINE double getSampleWhite();
  /**< Calculates a sample of white noise in the range -1.0...+1.0. */

  INLINE double getSampleLowpassWhite();
  /**< Calculates a sample of lowpass-filtered white noise. A gain factor is 
       applied to compensate for the loss of energy in the filter. */

  INLINE double getSampleHighpassWhite();
  /**< Calculates a sample of highpass-filtered white noise. A gain factor is 
       applied to compensate for the loss of energy in the filter. */

  INLINE double getSampleBandpassWhite();
  /**< Calculates a sample of bandpass-filtered white noise. A gain factor is 
       applied to compensate for the loss of energy in the filter. */

  INLINE double getSamplePink();
  /**< Calculates a sample of pink noise in the range, A gain factor is 
       applied to compensate for the loss of energy in the pinking-filter. */

  INLINE double getSampleLowpassPink();
  /**< Calculates a sample of lowpass-filtered pink noise. A gain factor is 
       applied to compensate for the loss of energy in the filter. */

  INLINE double getSampleHighpassPink();
  /**< Calculates a sample of highpass-filtered pink noise. A gain factor is 
       applied to compensate for the loss of energy in the filter. */

  INLINE double getSampleBandpassPink();
  /**< Calculates a sample of bandpass-filtered pink noise. A gain factor is 
       applied to compensate for the loss of energy in the filter. */

 protected:

  // embedded audio-modules:
  FourPoleFilter lowpassFilter, highpassFilter;

  double scale, offset, sampleRate;

 };

 //-----------------------------------------------------------------------------
 // from here: definitions of the functions to be inlined, i.e. all functions
 // which are supposed to be called at audio-rate (they can't be put into
 // the .cpp file):

 INLINE double NoiseGeneratorOld::getSampleWhite()
 {
  return scale * ((double) rand()) + offset;
 }

 INLINE double NoiseGeneratorOld::getSampleLowpassWhite()
 {
  return lowpassFilter.getSample(getSampleWhite());
 }

 INLINE double NoiseGeneratorOld::getSampleHighpassWhite()
 {
  return highpassFilter.getSample(getSampleWhite());
 }

 INLINE double NoiseGeneratorOld::getSampleBandpassWhite()
 {
  return highpassFilter.getSample(getSampleLowpassWhite());
 }

 INLINE double NoiseGeneratorOld::getSamplePink()
 {
  return scale * ((double) rand()) + offset;
 }

 INLINE double NoiseGeneratorOld::getSampleLowpassPink()
 {
  return lowpassFilter.getSample(getSamplePink());
 }

 INLINE double NoiseGeneratorOld::getSampleHighpassPink()
 {
  return highpassFilter.getSample(getSamplePink());
 }

 INLINE double NoiseGeneratorOld::getSampleBandpassPink()
 {
  return highpassFilter.getSample(getSampleLowpassPink());
 }

} // end namespace rosic

#endif // rosic_NoiseGeneratorOld_h
