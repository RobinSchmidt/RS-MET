#ifndef rosic_FiniteImpulseResponseFilter_h
#define rosic_FiniteImpulseResponseFilter_h

namespace rosic
{

/** This class implements a filter with finite impulse response.It uses a partitioned fast 
convolution algorithm to perform the filtering. */

class FiniteImpulseResponseFilter
{

public:

  //-----------------------------------------------------------------------------------------------
  // construction/destruction:

  /** Constructor. */
  FiniteImpulseResponseFilter();

  /** Destructor. */
  ~FiniteImpulseResponseFilter();

  //-----------------------------------------------------------------------------------------------
  // setup:

  /** Sets the sample-rate on which the filter should operate. */
  void setSampleRate(double newSampleRate);

  /** Selects the mode for the filter. */
  void setMode(int newMode);

  /** Sets the characteristic frequency of the filter. For lowpass and highpass filters, this is 
  the cutoff frequency, for bandpass, bandreject and peak filters the center frequency. */
  void setFrequency(double newFrequency);

  /** Sets the bandwidth of the filter in octaves (for bandpass and bandreject filters). */
  void setBandwidth(double newBandwidth);

  /** Sets the length of the impulse-response. */
  void setImpulseResponseLength(int newLength);

  /** Sets the window that is to be used for the design. @see WindowDesigner::windowTypes */
  void setWindowType(int newWindow);

  /** Passes an impulse-response to be used by this filter. */
  //void setImpulseResponse(double *newImpulseResponse, int newLength);

  // \todo setMagnitudeResponse - create a linear-phase FIR filter with the desired magnitude 
  // response
  // \todo maybe introduce a flag to the setters to avoid instant update of the coeffs 
  //  -> facilitates avoidance multiple re-calculations when several parameters change at once
  //  ..or maybe work with a "dirty" flag

  //-----------------------------------------------------------------------------------------------
  // inquiry (todo: make them all const):

  /** Returns true if the currently selected mode supports a bandwidth parameter. */
  bool hasCurrentModeBandwidthParameter() { return designer.hasCurrentModeBandwidthParameter(); }

  /** Returns the value of the complex transfer-function H(z) at the given z. */
  Complex getTransferFunctionAt(Complex z);

  void getMagnitudeResponse(double* frequencies, double* magnitudes, int numBins, 
    bool inDecibels = false, bool accumulate = false);


  int getKernelLength() const { return kernelLength; }

  const double* getKernelPointer() const { return h; }


  //-----------------------------------------------------------------------------------------------
  // audio-processing:

  /** Computes one output-sample. */
  INLINE double getSample(double in)
  {
    return convolver.getSample(in);
  }

  //===============================================================================================

protected:

  /** Triggers a re-calculation of the filter coefficients (i.e. the impulse response). */
  void updateCoefficients();

  FiniteImpulseResponseDesigner designer;
  ConvolverPartitioned convolver;

  double* h;        // the impulse-response
  int kernelLength; // length of the filter kernel in samples (aka impulse response)

  // maybe maxLength, to be set in the constructor

};

}

#endif 
