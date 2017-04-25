#ifndef DirectFormIir_h
#define DirectFormIir_h

#include "AudioModule.h"
#include "IirDesigner.h"

/**

This class implements a direct form IIR filter of arbitrary order (up to
maxOrder as defined in IirDesigner.h). The coefficients of the filter are
calculated by the embedded IirDesigner object.

*/

class DirectFormIir : public AudioModule  
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
 // construction/destruction:

	         DirectFormIir();  ///< Constructor.
	virtual ~DirectFormIir();  ///< Destructor.

 //---------------------------------------------------------------------------
 //parameter settings:

 virtual void setSampleRate(double newSamplerate);   
 ///< Overrides the setSampleRate() method of the AudioModule base-class.
 
 virtual void setSlope(int newSlope);        
 /**< Determines the order of the filter - 1 results in 6 dB/oct slope, 
      2->12dB/oct 3->18, 4->24 and so on. The actual order of the filter
      depends on its characteristic: in LPF's and HPF's the order is
      identical to the value of "slope", in BPFs and BRFs the order is
      twice as high. */

 virtual void setMode(int newMode);
 /**< Sets the mode of the filter. 4 modes are available: Lowpass = 1, 
      Highpass = 2, Bandpass=3, Bandreject=4. These modes are also 
      enumerated, such that the enum-labels can be used instead of raw 
      numbers. */

 virtual void setApproximationMethod(int newMethod);       
 /**< Chooses the approximation method. At the moment only the Butterworth
      approximation is supported. */

 virtual void setFreq1(double newFreq1);        
 /**< Sets the cutoff frequency for lowpass- and highpass-filters, or the 
      lower cutoff frequency for bandpass- and bandreject-filters. */

 virtual void setFreq2(double newFreq2);        
 /**< Sets the upper cutoff frequency for bandpass- and bandreject-filters. */

 //---------------------------------------------------------------------------
 //audio processing:

 virtual INLINE double getSample(double in);     
 ///< Calculates one onutput sample at a time.

 //---------------------------------------------------------------------------
 //others:

 virtual void reset(); 
 ///< Resets the internal buffers to zeros.

 virtual void getMagnitudeResponse(float *frequencies, float *magnitudes, 
                                   int numBins, bool inDecibels);
 /** Calculates the magnitudes of the frequency-response at the frequencies
     given in the array "frequencies" (in Hz) and stores them in the 
     array "magnitudes". Both arrays are assumed to be "numBins" long. 
     "inDecibels" indicates, if the frequency response should be returned in
     decibels. */

protected:

 static const int maxOrder = 24; // 24 is the maximum number of poles - the
                                   // slope-parameter will be restricted to 
                                   // the range 1 - 12

 //embedded audio-modules:

 IirDesigner filterDesigner;
 /**< The embedded filter designer object, used to calculate the filter 
      coefficients. */

 //parameters:
 intA slope, order, mode, method;
 doubleA freq1, freq2;

 //the filter coefficient arrays:
 doubleA ffCoeffs[maxOrder+1];
 doubleA fbCoeffs[maxOrder+1];

 //buffer for the past input and output samples:
 doubleA xBuffer[maxOrder+1];
 doubleA yBuffer[maxOrder+1];

};

//-----------------------------------------------------------------------------
//from here: definitions of the functions to be inlined, i.e. all functions
//which are supposed to be called at audio-rate (they can't be put into
//the .cpp file):

//----------------------------------------------------------------------------
//audio processing:
INLINE double DirectFormIir::getSample(double in)
{
 static doubleA out;
 static intA i;

 //shift input samples in the input buffer one sample back:
 for(i=order; i>=1; i--)
  xBuffer[i] = xBuffer[i-1];

 //put the current sample into the first position (x[n]) of the buffer:
 xBuffer[0] = in;

 //apply the feedforward part of the filter:
 out = 0;
 for(i=0; i<=order; i++)
  out += ffCoeffs[i] * xBuffer[i];

 //apply the feedback part of the filter:
 for(i=1; i<=order; i++)
  out -= fbCoeffs[i] * yBuffer[i];

 //shift the samples in the y-buffer one sample back and store the
 //recently calculated sample in the y[n-1] position for the next call:
 for(i=order; i>=2; i--)
  yBuffer[i] = yBuffer[i-1];
 yBuffer[1] = out;

 return out;
}


#endif // DirectFormIir_h
