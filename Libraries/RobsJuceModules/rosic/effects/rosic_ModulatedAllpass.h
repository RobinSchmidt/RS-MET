#ifndef rosic_ModulatedAllpass_h
#define rosic_ModulatedAllpass_h

//// rosic-indcludes:
//#include "../filters/rosic_EllipticSubBandFilterDirectForm.h"
//#include "../filters/rosic_LowpassHighpass.h"

namespace rosic
{

  /**

  This class applies distortion to a signal by means of a coefficient-modulated first order allpass 
  filter. The allpass coefficient at time instant n is given by a modulation signal m[n] which 
  itself is derived from the input signal by filtering it, mutliplying it by a factor and adding
  an offset and finally applying a tanh-shaped saturation function such that: 
  
  m[n] = tanh( offset + factor * filter(x[n]) )

  the output is then obtained via the pair of difference equations:

  y[n] = -m[n] * x[n] + w[n-1] 
  w[n] = x[n] + m[n] * y[n]

  References:
  -Jussi Pekonen: Coefficient-Modulated First-Order Allpass Filter As Distortion Effect (DAFX 2008)


  \todo: oversample, input filter, drive, dc
   
  */

  class ModulatedAllpass
  {

  public:

    //---------------------------------------------------------------------------------------------
    // construction/destruction:

    /** Constructor. */
    ModulatedAllpass();

    /** Destructor */
    ~ModulatedAllpass();

    //---------------------------------------------------------------------------------------------
    // parameter settings (set-functions):

    /** Sets the samplerate. */
    void setSampleRate(double newSampleRate);

    /** Sets the factor by which the filtered input signal is multiplied in the modulation 
    signal. */
    void setFactor(double newFactor) { factor = newFactor; }

    /** Sets the offset for the modulation signal. */
    void setOffset(double newOffset) { offset = newOffset; }



    //---------------------------------------------------------------------------------------------
    // audio processing:

    /** Calculates one output stereo sample-frame at a time. */
    INLINE void getSampleFrameStereo(double *inOutL, double *inOutR);

    //---------------------------------------------------------------------------------------------
    // others:

    /** Resets the internal states of the subband filters. */
    void reset();

    //=============================================================================================

  protected:

    double wL, wR;                // filter state
    double factor, offset;

    /*
    double dry, wet;
    double  driveFactor;              // global gain factor at the input stage
    double  gainFactors[numShapers];  // gain factors for the individual harmonics
    bool   gainIsZero[numShapers];   // array of flags to indicate zero gain
    bool   chebychevMode;            // flag to indicate that chebychev polynomials should be used

    EllipticSubBandFilterDirectForm subbandFiltersL[numShapers], subbandFiltersR[numShapers];
    LowpassHighpass inFilterL, inFilterR, outFilterL, outFilterR;
    */

  };

  //-----------------------------------------------------------------------------------------------
  // from here: definitions of the functions to be inlined, i.e. all functions which are supposed
  // to be called at audio-rate (they can't be put into the .cpp file):

  INLINE void ModulatedAllpass::getSampleFrameStereo(double *inOutL, double *inOutR)
  {
    double x = SQRT2_INV * (*inOutL + *inOutR);  
      // todo: apply drive and filter to x
    double m = RAPT::rsTanhApprox( factor*x + offset );

    // test:
    //double x = *inOutL;  
    //double m = factor*x + offset;

    double yL = -m*(*inOutL) + wL;
    double yR = -m*(*inOutR) + wR;

    wL = *inOutL + m*yL;
    wR = *inOutR + m*yR;

    *inOutL = yL;
    *inOutR = yR;
  }

} // end namespace rosic

#endif // #ifndef rosic_ModulatedAllpass_h
