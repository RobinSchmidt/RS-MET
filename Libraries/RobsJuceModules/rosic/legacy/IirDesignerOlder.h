#ifndef IirDesigner_h
#define IirDesigner_h

#include "math.h"
#include "Complex.h"
#include <stdio.h>
#include "AudioModule.h"

/**

This is a class which designs (calculates the coefficients for) scientific
IIR filters. The user can set different parameters such as the sampleRate,
filter-mode, approximation method, cutoff frequency (or frequencies for
bandpass and bandreject mode), etc. When the IirDesigner object is set up,
the user calls the getDirectFormCoeffs()-function with the pointers to the
arrays where the filter coefficients should be written into.

ToDo: 
the frequency-transformation assume the old ordering of the poles and zeros - adapt
them to the new ordering (which has been introduced since the prototype poles and zeros
were hardcoded).


*/

class IirDesigner  
{

public:

 /** This is an enumeration of the 4 available filter modes. */
 enum modes
 {
  BYPASS = 0,   ///< bypass
  LOWPASS,      ///< lowpass filter mode
  HIGHPASS,     ///< highpass filter mode
  BANDPASS,     ///< bandpass filter mode
  BANDREJECT    ///< bandreject filter mode
 };

 /** This is an enumeration of the available approximation methods - at the
     moment only Butterworth is supported. */
 enum approximationMethods
 {
  BUTTERWORTH = 1,  ///< Butterworth approximation method
  CHEBYCHEV,        ///< NOT implemented. Chebychev approximation method
  INV_CHEBYCHEV,    ///< NOT implemented. inverse Chebychev approximation method
  ELLIPTIC,         ///< NOT implemented.
  BESSEL,           ///< NOT implemented.
  PAPOULIS          ///< NOT implemented.
 };

 //---------------------------------------------------------------------------
 //construction/destruction:

	         IirDesigner();  ///< Constructor.
	virtual ~IirDesigner();  ///< Destructor.

 //---------------------------------------------------------------------------
 //parameter settings:

 virtual void setSampleRate(flt64 newSampleRate);
 ///< Sets the sample-rate of the filter.

 virtual void setMode(int64 newMode);       
 ///< Chooses one of the modes (LP, HP, BP, BR) as enumerated above.

 virtual void setSlope(int64 newSlope);      
 /** Determines the order of the filter - 1 results in 6 dB/oct slope, 
     2->12dB/oct 3->18, 4->24 and so on. The actual order of the filter
     depends on its characteristic: in LPF's and HPF's the order is
     identical to the value of "slope", in BPFs and BRFs the order is
     twice as high. */

 //virtual void setMethod(int64 newMethod);     
 /** Chooses one of the approximation methods as enumerated above. At 
     the moment only Butterworth is supported. */

 virtual void setFreq1(flt64 newFreq1);      
 /**< Sets the cutoff frequency for lowpass- and highpass-filters, or the 
      lower cutoff frequency for bandpass- and bandreject-filters. */

 virtual void setFreq2(flt64 newFreq2);
 /**< Sets the upper cutoff frequency for bandpass- and bandreject-filters. */

 //virtual void setPassbandRipple(flt64 newPassbandRipple);   
 /**< NOT yet implemented. Sets the allowed Passband ripple in dB. */

 //virtual void setStopbandRipple(flt64 newStopbandRipple);
 /**< NOT yet implemented. Sets the allowed Passband ripple in dB. */


 //---------------------------------------------------------------------------
 // coefficient calculation:

 virtual void getDirectFormCoeffs(flt64 *FeedforwardCoeffs, 
                                  flt64 *FeedbackCoeffs);
 /**< Calculates filter coefficients and stores them in the passed arrays.
      The feedforward coefficients (b0,...,bN) will be written into the array 
      "FeedforwardCoeffs", and the feedback coefficients (a0,...,aN) where
      a0=1 will be written into the array "FeedbackCoeffs". The calling
      function has to make sure, that these arrays are large enough to hold 
      "slope+1" values for LPFs and HPFs and 2*slope+1 for BPFs and BRFs 
      where "slope" is the parameter passed to the setSlope() function when
      the filter-designer was set up. These coefficients can be used to
      implement the direct form difference equation:

      \f$
      y[n] = b0*x[n] + b1*x[n-1] + b2*x[n-2] + ... + bN*x[n-N]
                     - a1*y[n-1] - a2*y[n-2] - ... - aN*y[n-N]  \f$

      ATTENTION: in some dsp-textbooks the signs of the feedback-coeffs are
      reversed and some dsp-textbooks call the feedforward coeffs a-coeffs
      and the feedback-coeffs b-coeffs, this class calculates the coeffients
      for use with the notation stated above!
 */

 // to be implemented:
 virtual void getBiquadCascadeCoeffs(flt64 *b0, flt64 *b1, flt64 *b2, 
                                                flt64 *a1, flt64 *a2);
 /**< NOT implemented yet!
      Calculates the coefficients for a cascade of biquad filter units. The
      arrays are for the b0, b1, b2, a0, a1, a2 coefficients for the
      individual biquad stages. They can then be used like this:

      y1[n] = b0[0]*x[n] + b1[0]*x[n-1]  + b2[0]*x[n-2]
                         - a1[0]*y1[n-1] - a2[0]*y1[n-2] //this was the 
                                                         //first biquad stage

      y2[n] = b1[1]*y1[n] + b1[1]*y1[n-1] + b2[1]*y1[n-2]
                          - a1[1]*y2[n-1] + a2[1]*y2[n-2] //this was the 
                                                          //second biquad stage
      and so on...

      The individual biquad stages are each normalized such that their ouput 
      level is 1 at z=1 (for LPFs and BRFs), z=-1 (for HPFs) or 
      z=exp(j*omega_c) (for BPFs with omega_c being the normalized radian 
      center-frequency. This makes it meaningful to tap outputs between the
      individual stages. The biquads are ordered in such a way, that....
 */

 virtual void calculateAndPrint();
 /**< This function is for debugging purposes. It goes through all the 
      coefficient calculation steps and prints the results to screen. */


 //===========================================================================

protected:

 //internal parameter variables:

 static const int64 maxOrder = 24; // 24 is the maximum number of poles - the
                                   // slope-parameter will be restricted to 
                                   // the range 1 - 12

 flt64A sampleRate, 
        freq1,
        freq2, 
        passbandRipple, 
        stopbandRipple;

 int64A mode, 
        slope, 
        method;

 int64A numPolesLpfAna; // number of poles of the analog unit lpf

 int64A numZerosLpfAna; // number of zeros of the analog unit lpf

 int64A numPolesAna;    // number of poles of the (frequency transformed)
                        // analog filter

 int64A numZerosAna;    // number of zeros of the (frequency transformed)
                        // analog filter 
 
 // when BPF- or BRF- filters are used, then there are twice as much poles 
 // and zeros as in the LPF/HPF-case

 int64A numPolesDig; // number of poles in the digital filter   
 int64A numZerosDig; // number of zeros in the digital filter

 // Actually a filter has always the same number of poles and zeros and the
 // digital version has the same number of poles and zeros as the analog
 // prototype filter. However, in the H(s) transfer function the numerator
 // may not have explicit zeros at all. In this case all zeros of H(s) are
 // at infinity - these should be mapped to z = -1. So the variable 
 // numZerosAna actually represents the number of explicit zeros in the
 // numerator of H(s) and may be different from numPolesAna and numZerosDig. 
 // That's why we use 4 variables here.

 //internal variables (calculated from the parameters):
 flt64A omegaDig1, omegaDig2,
        omegaAna1, omegaAna2, // normalized (pre-warped) radian frequency(ies)
        omegaAnaC, bwAna;     // center freq and bandwidth for BPF and BRF

 flt64A gainFactor;

 Complex unitLpfPolesAna[maxOrder]; // array of the s-plane poles of the unit
                                    // cutoff analog LPF

 Complex unitLpfZerosAna[maxOrder]; // array of the s-plane zeros of the unit
                                    // cutoff analog LPF

 Complex polesAna[maxOrder]; // array of the s-plane poles of the 
                             // frequency-tranformed analog filter  

 Complex zerosAna[maxOrder]; // array of the s-plane zeros of the 
                             // frequency-tranformed analog filter  

 Complex polesDig[maxOrder]; // array of the z-plane poles of the 
                             // digital filter 

 Complex zerosDig[maxOrder]; // array of the z-plane zeros of the
                             // digital filter 

 flt64A ffCoeffs[maxOrder+1]; // array of the feedforward 
                              // coefficients (b0,...,bN)

 flt64A fbCoeffs[maxOrder+1]; // array of the feedback
                              // coefficients (a0,...,aN) where a0=1


 //---------------------------------------------------------------------------
 // functions for internal calculations:

 void prewarpFreqs();         
 /**< Prewarps the desired cutoff frequencies of the digital filter to the
      corresponding analog frequencies. */

 void calcButterUnitLpfAna(); 
 /**< Calculates the poles and zeros of an analog prototype lowpass-filter
      with unit cutoff the Butterworth approximation method. */

 void assignPrototypePolesAndZeros();
 /**< Assigns the pole- and zero-positions of the analog prototype lowpass 
      filter with unit cutoff-frequency. The positions were calculated with 
      the buttap, cheb1ap, cheb2ap, ellipap and besselap function from MatLabs
      Signal Processing Toolbox and hardcoded into the function. */

 void frequencyTransform();   
 /**< Performs the frequency transform in the s-plane from the unit
      lowpass-filter to the desired characteristic and cutoff frequency(ies).
      That is: lowpass->lowpass, lowpass->highpass, lowpass-bandpass or
      lowpass-bandreject transform. */

 void bilinearTransform();    
 /**< Transforms the poles and zeros from the s-plane to the z-plane by means
      of the bilinear transform. That transforms the analog prototype filter
      into a digital filter. */

 void calcGainFactor();       
 /**< Calculates an overall gain factor for the filter which should be applied
      to the feedforward coefficients. The factor is calculated such that the
      magnitude response at some representative frequency in the passband is
      unity (this frequency is chosen to be dc for LPF's, Nyquist for HPF's,
      the center-frequency for BPF's and dc for BRF's) */

 void calcDirectFormCoeffs(); 
 /**< Calculates the direct form filter coefficients from the poles and zeros
      in the z-plane via a recursively applied convolution of the linear
      factors (polyomials are multiplied by convolving their coefficients). */

 void calcBiquadCascadeCoeffs();   //...work in progress

 void convolve(Complex *Seq1, int64 Length1,
               Complex *Seq2, int64 Length2,
               Complex *Result); 
 /**< Convolves two sequences of complex numbers of Length1 and Length2. The
      result is stored in the array *Result which should be at least
      Length1+Length2-1 long. */

};

#endif // IirDesigner_h
