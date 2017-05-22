#ifndef rosic_ToneControl_h
#define rosic_ToneControl_h

//// rosic-indcludes:
//#include "../math/rosic_ElementaryFunctionsReal.h"

namespace rosic
{

 /**

 This class combines a first order low-shelving and a first order high-shelving
 filter into one biquad-section to be used as bass/treble control.

 */

 class ToneControl
 {

 public:

  //---------------------------------------------------------------------------
  // construction/destruction:

  ToneControl();   ///< Constructor.
  ~ToneControl();  ///< Destructor.

  //---------------------------------------------------------------------------
  // parameter settings:

  void setSampleRate(double newSampleRate);
  ///< Sets the sample-rate for this filter.

  void setLowCornerFreq(double newLowCornerFreq);
  /**< Sets the corner frequency of the low-shelving filter - the corner 
       frequency is defined to be the frequency whe the gain in dB is at it's
       half-point. */

  void setLowGain(double newLowGain);
  /**< Sets the gain of the low-shelving filter in dB */

  void setHighCornerFreq(double newHighCornerFreq);
  /**< Sets the corner frequency of the high-shelving filter. */

  void setHighGain(double newHighGain);
  /**< Sets the gain of the high-shelving filter in dB. */


  //---------------------------------------------------------------------------
  // audio processing:

  INLINE double getSample(double in);
  ///< Calculates one output sample at a time.

  //---------------------------------------------------------------------------
  //others:

  void resetBuffers();
  ///< Resets the internal buffers of the filter to zero.

  //===========================================================================

 protected:

  // buffering:
  doubleA out;
  doubleA x1, x2;  // past input values
  doubleA y1, y2;  // past output values

  // overall 2-pole filter-coefficients:
  doubleA b0, b1, b2; // feedforward coeffs
  doubleA     a1, a2; // feedback coeffs

  // 1-pole filter-coefficients:
  doubleA b0Ls, b1Ls, a1Ls; // low-shelf coeffs
  doubleA b0Hs, b1Hs, a1Hs; // high-shelf coeffs

  // filter parameters:
  doubleA lowCornerFreq;
  doubleA lowGain;
  doubleA highCornerFreq;
  doubleA highGain;

  doubleA sampleRate;
  doubleA sampleRateRec;  // reciprocal of the sampleRate
  
  // private member functions:
  void calcLowShelfCoeffs();  // calculates the low-shelving coeffs
  void calcHighShelfCoeffs(); // calculates the low-shelving coeffs
  void calcBiquadCoeffs();    // combines the two first order filters into one
                              // biquad filter
 };

 //-----------------------------------------------------------------------------
 // from here: definitions of the functions to be inlined, i.e. all functions
 // which are supposed to be called at audio-rate (they can't be put into
 // the .cpp file):

 INLINE double ToneControl::getSample(double in)
 {
  // calculate output-sample:
  out = b0*in + b1*x1 + b2*x2 + a1*y1 + a2*y2;

  // update buffer-variables:
  x2 = x1;
  x1 = in;
  y2 = y1;
  y1 = out;

  return out;
 }

} // end namespace rosic

#endif // rosic_ToneControl_h
