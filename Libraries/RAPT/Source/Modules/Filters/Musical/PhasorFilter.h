#ifndef RAPT_PHASORFILTER_H_INCLUDED
#define RAPT_PHASORFILTER_H_INCLUDED

/** This is a filter based on a spiraling phasor in the xy-plane... 

This class is just a stub - it's still under construction..
*/

template<class TSig, class TPar> // signal, parameter types
class PhasorFilter
{

public:

  /** \name Construction/Destruction */

  /** Constructor. */
  PhasorFilter();


  /** \name Setup */

  /** Sets the sample rate in Hz. */
  void setSampleRate(TPar newSampleRate);

  /** Sets the characteristic frequency in Hz. */
  void setFrequency(TPar newFrequency);

  /** Sets the decay time in seconds. You can also pass infinity in which case the filter will 
  produce a steady sinusoid as impulse response - so you can use it as oscillator as well. */
  void setDecayTime(TPar newDecay);


  /** \name Audio Processing */

  /** Processes a single sample frame... */
  inline void processFrame(TSig *x, TSig *y);

  /** Computes an output sample from a given input sample */
  inline TSig getSample(TSig in);


  /** \name Misc */

  /** Resets the internal states to zero. To be overriden in subclasses. */
  void reset();


protected:

  /** \name Internal Functions */
                        
  /** Updates the internal coefficients from the parameters. */
  void updateCoefficients(); 

  /** \name Data */

  // state variables:
  TSig xOld, yOld;           // previous values of x- and y-coordinate

  // internal coefficients:
  TPar Axx, Axy, Ayx, Ayy;   // spiraling matrix elements

  // user parameters:
  TPar sampleRate;           // samplerate in Hz
  TPar frequency;            // filter frequency in Hz
  TPar decay;                // decay time constant in seconds

};

#endif
