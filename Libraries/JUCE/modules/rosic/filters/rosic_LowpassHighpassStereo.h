#ifndef rosic_LowpassHighpassStereo_h
#define rosic_LowpassHighpassStereo_h

// standard-library indcludes:
#include <math.h>

// rosic-indcludes:
#include "../basics/GlobalDefinitions.h"

namespace rosic
{

 /**

 This class combines a first order lowpass-, and a first order highpass-filter
 into a single object. The 2 filters, which are originally in a series 
 connection, are implemented as a single biquad-stage.

 */

 class LowpassHighpassStereo
 {

 public:

  //---------------------------------------------------------------------------
  // construction/destruction:

  LowpassHighpassStereo();   ///< Constructor.
  ~LowpassHighpassStereo();  ///< Destructor.

  //---------------------------------------------------------------------------
  // parameter settings:

  void setSampleRate(double newSampleRate);
  ///< Sets the sample-rate for this filter

  void setLpfCutoff (double newLpfCutoff);
  ///< Sets the cutoff frequency of the lowpass-filter.

  void setHpfCutoff (double newHpfCutoff);
  ///< Sets the cutoff frequency of the highpass-filter.

  //---------------------------------------------------------------------------
  // audio processing:

  INLINE void getSampleFrameStereo(double *inL, double *inR, 
                                   double *outL, double *outR);
  ///< Calculates one output sample-frame at a time.

  //---------------------------------------------------------------------------
  // others:

  void resetBuffers();
  ///< Resets the internal buffers of the filter to zero.

  //===========================================================================

 protected:

  // buffering:
  double x1L,  x1R;           // past input values delayed by one sample
  double x2L,  x2R;           // past input values delayed by two samples
  double y1L,  y1R;           // past output values delayed by one sample
  double y2L,  y2R;           // past output values delayed by two samples

  // biquad-coefficients:
  double b0, b1, b2;       // feedforward coeffs
  double     a1, a2;       // feedback coeffs

  // filter parameters:
  double lpfCutoff;
  double hpfCutoff;

  double sampleRate;
  double sampleRateRec;    // reciprocal of the sampleRate
  
  // private member functions:
  void calcCoeffs();     
   // calculates the filter coefficients from filter parameters. The 
   // design-equations for the first order filters come from "The Scientist and
   // Engineers Guide to Digital Signal Processing" (www.dspguide.com)

 };

 //-----------------------------------------------------------------------------
 // from here: definitions of the functions to be inlined, i.e. all functions
 // which are supposed to be called at audio-rate (they can't be put into
 // the .cpp file):

 INLINE void LowpassHighpassStereo::getSampleFrameStereo(double *inL, 
                                                         double *inR, 
                                                         double *outL, 
                                                         double *outR)
 {
  double tmpL, tmpR;

  // calculate left output-sample:
  tmpL = b0 * (*inL) + b1*x1L + b2*x2L
                     + a1*y1L + a2*y2L 
                     + TINY;          // to avoid denorm problems

  // update left buffer-variables:
  x2L = x1L;
  x1L = *inL;
  y2L = y1L;
  y1L = tmpL;  

  // calculate right output-sample:
  tmpR = b0 * (*inR) + b1*x1R + b2*x2R
                     + a1*y1R + a2*y2R 
                     + TINY;          // to avoid denorm problems

  // update right buffer-variables:
  x2R = x1R;
  x1R = *inR;
  y2R = y1R;
  y1R = tmpR;  

  // store the output-values in their slots:
  *outL = tmpL;
  *outR = tmpR;
 }

} // end namespace rosic

#endif // LowpassHighpassStereo_h
