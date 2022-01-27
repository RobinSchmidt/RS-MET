#ifndef rosic_SlopeFilter_h
#define rosic_SlopeFilter_h

namespace rosic
{

/** This class implements a filter that approximates a (user-adjustable) constant slope in (dB/oct) 
across the whole audible range using two second order shelving filters. 

ToDo: 
-templatize and drag over to RAPT 
-use an SVF-based implementation (but maybe in a subclass - it's more expensive and may not always 
 be needed - it's mostly for better time-variant behavior - maybe a state-vector based version 
 could be useful, too).  */

class SlopeFilter
{

public:

  //-----------------------------------------------------------------------------------------------
  // construction/destruction:

  /** Constructor. */
  SlopeFilter();

  /** Destructor. */
  ~SlopeFilter();

  //-----------------------------------------------------------------------------------------------
  // setup:

  /** Sets the sample-rate on which the filter should operate. */
  void setSampleRate(double newSampleRate);

  /** Selects the slope (in dB/oct) for the filter. */
  void setSlope(double newSlope);

  //-----------------------------------------------------------------------------------------------
  // audio-processing:

  /** Computes one output-sample. */
  INLINE double getSample(double in);

  /** Computes one stereo output sample frame. */
  INLINE void getSampleFrameStereo(double *inOutL, double *inOutR);
  // rename to processFrameStereo or processFrame

  //-----------------------------------------------------------------------------------------------
  // others:

  /** Resets the internal buffers. */
  void reset();

  //===============================================================================================

protected:

  /** Triggers a re-calculation of the filter coefficients. */
  void updateCoefficients();

  double s1b0, s1b1, s1b2, s1a1, s1a2; // filter coefficients for biquad stage 1
  double s2b0, s2b1, s2b2, s2a1, s2a2; // filter coefficients for biquad stage 2

  double s1x1L, s1x2L, s1y1L, s1y2L;   // buffer-variables stage 1, left channel
  double s2x1L, s2x2L, s2y1L, s2y2L;   // buffer-variables stage 2, left channel

  double s1x1R, s1x2R, s1y1R, s1y2R;   // buffer-variables stage 1, right channel
  double s2x1R, s2x2R, s2y1R, s2y2R;   // buffer-variables stage 2, right channel

  // user parameters:
  double sampleRate;
  double slope;

};

//-------------------------------------------------------------------------------------------------
// inlined functions:

INLINE double SlopeFilter::getSample(double in)
{
  double tmp, tmp2;

  // stage 1:
  tmp   = (s1b0*in+TINY) + (s1b1*s1x1L + s1b2*s1x2L) + (s1a1*s1y1L + s1a2*s1y2L); // parentheses facilitate out-of-order execution
  s1x2L = s1x1L;
  s1x1L = in;
  s1y2L = s1y1L;
  s1y1L = tmp;

  // stage 2:
  tmp2  = (s2b0*tmp) + (s2b1*s2x1L + s2b2*s2x2L) + (s2a1*s2y1L + s2a2*s2y2L);
  s2x2L = s2x1L;
  s2x1L = tmp;
  s2y2L = s2y1L;
  s2y1L = tmp2;

  return tmp2;
}

INLINE void SlopeFilter::getSampleFrameStereo(double *inOutL, double *inOutR)
{
  double tmp, tmp2;

  // left channel:
  tmp   = (s1b0 * *inOutL + TINY) + (s1b1*s1x1L + s1b2*s1x2L) + (s1a1*s1y1L + s1a2*s1y2L);
  s1x2L = s1x1L;
  s1x1L = *inOutL;
  s1y2L = s1y1L;
  s1y1L = tmp;

  tmp2  = (s2b0*tmp) + (s2b1*s2x1L + s2b2*s2x2L) + (s2a1*s2y1L + s2a2*s2y2L);
  s2x2L = s2x1L;
  s2x1L = tmp;
  s2y2L = s2y1L;
  s2y1L = tmp2;

  *inOutL = tmp2;

  // right channel:
  tmp   = (s1b0 * *inOutR + TINY) + (s1b1*s1x1R + s1b2*s1x2R) + (s1a1*s1y1R + s1a2*s1y2R);
  s1x2R = s1x1R;
  s1x1R = *inOutR;
  s1y2R = s1y1R;
  s1y1R = tmp;

  tmp2  = (s2b0*tmp) + (s2b1*s2x1R + s2b2*s2x2R) + (s2a1*s2y1R + s2a2*s2y2R);
  s2x2R = s2x1R;
  s2x1R = tmp;
  s2y2R = s2y1R;
  s2y1R = tmp2;

  *inOutR = tmp2;
}

}

#endif 
