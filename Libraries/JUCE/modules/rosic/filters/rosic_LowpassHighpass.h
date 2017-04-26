#ifndef rosic_LowpassHighpass_h
#define rosic_LowpassHighpass_h

// rosic-indcludes:
#include "../math/rosic_ElementaryFunctionsReal.h"
#include "rosic_BiquadMonoDF1.h"

namespace rosic
{

 /**

 This class combines a first order lowpass-, and a first order highpass-filter into a single 
 object. The 2 filters, which are originally in a series connection, are implemented as a 
 single biquad-stage.

 */

  class LowpassHighpass : public BiquadMonoDF1
 {

 public:

  //-----------------------------------------------------------------------------------------------
  // construction/destruction:

  /** Constructor. */
  LowpassHighpass();   

  /** Destructor. */
  ~LowpassHighpass();  

  //-----------------------------------------------------------------------------------------------
  // parameter settings:

  /** Sets the sample-rate for this filter */
  void setSampleRate(double newSampleRate);

  /**Sets the cutoff frequency of the lowpass-filter. */
  void setLowpassCutoff(double newCutoff);

  /** Sets the cutoff frequency of the highpass-filter. */
  void setHighpassCutoff(double newCutoff);

  //-----------------------------------------------------------------------------------------------
  // inquiry:

  /** Returns the cutoff frequency of the lowpass-filter. */
  double getLowpassCutoff() const { return lpfCutoff; }

  /** Returns the cutoff frequency of the highpass-filter. */
  double getHighpassCutoff() const { return hpfCutoff; }

  //===============================================================================================

 protected:

   /** Calculates the filter coefficients from filter parameters. The design-equations for the 
   first order filters come from "The Scientist and Engineers Guide to Digital Signal Processing" 
   (www.dspguide.com). */
  void calcCoeffs();  

  // filter parameters:
  double lpfCutoff;
  double hpfCutoff;
  double sampleRateRec;    // reciprocal of the sampleRate

 };

} // end namespace rosic

#endif // LowpassHighpass_h
