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

  /** Sets the nonlinear mapper object that is used internally to mess with the filter state after
  the linear part of the filter was applied. If you pass a nullptr, the filter will revert to the 
  linear case. */
  void setStateMapper(Mapper2D<TSig> *newMapper);


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

  // embedded objects:
  Mapper2D<TSig> *mapper;

};

//=================================================================================================

/**

\todo make different kinds of mappers, for example one that expresses the new phase explicitly as
function of the old phase and the magnitude.
*/

template<class T>
class PhasorStateMapper : public Mapper2D<T>
{

public:

  PhasorStateMapper();


  /** \name Setup */

  void setSameSquare(T newCoeff);

  void setOtherSquare(T newCoeff);

  void setCrossProduct(T newCoeff);

  void setAddedConstant(T newCoeff);

  void setPreNormalizeSaturation(T newCoeff);

  void setPostNormalizeSaturation(T newCoeff);


  /** \name Mapping */

  virtual void map(T *x, T *y) override;

protected:

  // add parameters
  T a, b, c, d, sat1, sat2;

};

#endif
