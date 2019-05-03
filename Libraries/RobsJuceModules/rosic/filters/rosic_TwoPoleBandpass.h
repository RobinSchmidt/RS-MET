#ifndef rosic_TwoPoleBandpass_h
#define rosic_TwoPoleBandpass_h

//// rosic-indcludes:
//#include "../math/rosic_ElementaryFunctionsReal.h"

namespace rosic
{

  /**

  This is a 2-pole bandpass-filter which mimics a chain of a simple RC-lowpass filter and a 
  RC-highpass filter. One could use two objects of class OnePole as well, but it's more
  convenient to have both filters in one single object.

  */

  class TwoPoleBandpass
  {

  public:

    //---------------------------------------------------------------------------------------------
    // construction/destruction:

    TwoPoleBandpass();   ///< Constructor.
    ~TwoPoleBandpass();  ///< Destructor.

    //---------------------------------------------------------------------------------------------
    // parameter settings:

    void setSampleRate(double newSampleRate);
    ///< Sets the setSampleRate() for thsi filter.

    void setLpfCutoff(double newLpfCutoff);
    ///< Sets the cutoff-frequency of the lowpass-filter.

    void setHpfCutoff(double newHpfCutoff);
    ///< Sets the cutoff-frequency of the highpass-filter.

    //---------------------------------------------------------------------------------------------
    // inquiry:

    double getLpfCutoff();
    /**< Returns the cutoff frequency of the lowpass filter. */

    double getHpfCutoff();
    /**< Returns the cutoff frequency of the highpass filter. */

    //---------------------------------------------------------------------------------------------
    // audio processing:

    INLINE double getSample(double in);
    ///< Calculates one output sample at a time.

    //---------------------------------------------------------------------------------------------
    // others:

    void resetBuffers();

    //=============================================================================================

  protected:

    doubleA sampleRate;
    doubleA sampleRateRec;  // reciprocal of the sampleRate

    // filter parameters:
    doubleA lpfCutoff;
    doubleA hpfCutoff;

    // lowpass-filter coefficients:
    doubleA b0Lpf; // feedforward coeffs
    doubleA b1Lpf;
    doubleA a1Lpf; // feedback coeff

    // highpass-filter coefficients:
    doubleA b0Hpf; // feedforward coeffs
    doubleA b1Hpf;
    doubleA a1Hpf; // feedback coeff

    // buffering:
    doubleA x_1Lpf;  // past input sample of LPF (x[n-1])
    doubleA y_1Lpf;  // past output sample of LPF (y[n-1])
    doubleA x_1Hpf;  // past input sample of HPF (x[n-1])
    doubleA y_1Hpf;  // past output sample of HPF (y[n-1])

    // internal functions:
    void calcCoeffs();     // calculates filter coefficients from
    // filter parameters

  };

  //---------------------------------------------------------------------------------------------
  // from here: definitions of the functions to be inlined, i.e. all functions
  // which are supposed to be called at audio-rate (they can't be put into
  // the .cpp file):

  INLINE double TwoPoleBandpass::getSample(double in)
  {
    static doubleA lpfOut, hpfOut;

    // lowpass-filtering:
    lpfOut = b0Lpf*in + b1Lpf*x_1Lpf + a1Lpf*y_1Lpf;

    // update the lpf buffer variables:
    x_1Lpf = in;
    y_1Lpf = lpfOut;

    // highpass-filtering:
    hpfOut = b0Hpf*lpfOut + b1Hpf*x_1Hpf + a1Hpf*y_1Hpf;

    // update the hpf buffer variables:
    x_1Hpf = lpfOut;
    y_1Hpf = hpfOut;

    return hpfOut;
  }

} // end namespace rosic

#endif // rosic_TwoPoleBandpass_h
