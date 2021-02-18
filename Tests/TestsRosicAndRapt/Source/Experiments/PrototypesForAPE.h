#ifndef PROTOTYPES_FOR_APE_H_INCLUDED
#define PROTOTYPES_FOR_APE_H_INCLUDED

/* This is the header for some prototype DSP classes that need to be accessible from the 
TestsRosicAndRapt project and also from APE scripts. This is necessarry because APE 
currently can't use rosic and therefore also not the rs_testing module due to lack of 
support for some standard C/C++ headers that are used in rosic (intrin.h, etc.). So, currently we 
are limited to use only RAPT classes in APE but not rosic and also not the Prototypes defined in
rs_testing.  */



//=================================================================================================

template<class T>
class rsShepardToneGenerator
{

public:

  //-----------------------------------------------------------------------------------------------
  // \name Setup

  /** Sets the sample rate at which this object runs. This determines the time increment for our
  phasor */
  inline void setSampleRate(T newRate) { dt = T(1) / newRate; }

  /** Sets the center of the bell-shaped spectral envelope as a midi pitch value. */
  inline void setCenterPitch(T newPitch);

  /** Sets the width of the bell-shaped spectral envelope in semitones. */
  inline void setPitchWidth(T newWidth);

  /** Sets the floor of the Gaussian bell curve, below which it is considered to be close enough to
  zero to switch the sine off. The actual bell is shifted and scaled, so as to actually hit zero at
  the boundaries and one in the middle, but this setting here changes the overall shape. Typical 
  values are 0.1...0.0001, maybe log-scaled - perhaps a user parameter should set this up in dB. 
  This here sets the raw value, however. */
  inline void setBellFloor(T newFloor);

  /** Sets the start pitch. This is the pitch of the center sinusoid that is heard immediately 
  after a reset. */
  inline void setStartPitch(T newPitch) { startPitch = newPitch; }

  /** Sets the ascension speed in semitones per second. Negative values will result in descending 
  pitch. */
  inline void setSpeed(T newSpeed)
  {
    pitchDelta = newSpeed;
    //pitchDelta = newSpeed * T(12); // is that correct?
  }



  //-----------------------------------------------------------------------------------------------
  // \name Inquiry

  /** Computes the desired gain factor for a given frequency in Hz. */
  inline T getGainForFrequency(T f) { return getGainForPitch(rsFreqToPitch(f)); }
  // avoid using that in realtime code - may be useful for plotting

  inline T getGainForPitch(T pitch);

  //-----------------------------------------------------------------------------------------------
  // \name Processing

  /** Produces one output sample at a time. */
  inline T getSample();

  /** Resets the internal state. */
  inline void reset() { t = T(0); pitch = startPitch; }

  /** Updates the coefficients for the spectral envelope computation. Called from getSample if 
  coeffsDirty is true. */
  inline void updateCoeffs();


protected:

  //inline void updateCutoffs()

  /** Conversion of physical frequency in Hz to radian frequency omega = 2*pi*f/fs. */
  inline T freqToOmega(T freq)  { return T(2.0*PI) * dt * freq; }
  // maybe use a const member "tau" for 2*pi - use 36 digits - that should be good enough for
  // quad precision - in sage: n(2*pi, digits=36), see
  // https://en.wikipedia.org/wiki/Quadruple-precision_floating-point_format

  //inline T omegaToFreq(T omega) { return omega / (T(2.0*PI) * dt); }

  T t     = T(0);  // Current time, normalized to 0..1. Our phasor for the reference sine.
  T dt    = T(0);  // Time increment per sample for t. Equal to 1/sampleRate
  T phase = T(0);  // Current phase of the refenrece sinusoid

  T argScale;      // Scales the argument of the Gaussian
  T resShift;      // Shifts the result of the Gaussian
  T resScale;      // Scales the shifted result of the Gaussian

  T bellFloor   = 0.01; // floor of the bell curve
  T centerPitch = 64;   // center of the bell
  T pitchWidth  = 72;   // width of the bell from left to right boundary
  T startPitch  = 64;   // Pitch of reference sine after reset
  T pitch       = 64;   // Current pitch of reference sine
  T pitchDelta  = 0;    // maybe this should be in semitones or octaves per second

  bool coeffsDirty = true;

  //T wRef;      // current reference omega - should later be subject to automatic raising/falling
  // change name

};

template<class T>
inline void rsShepardToneGenerator<T>::setCenterPitch(T newPitch) 
{ 
  if(newPitch != centerPitch) {
    centerPitch = newPitch;
    coeffsDirty = true; }
}

template<class T>
inline void rsShepardToneGenerator<T>::setPitchWidth(T newWidth) 
{ 
  if(newWidth != pitchWidth) {
    pitchWidth = newWidth;
    coeffsDirty = true; }
}

template<class T>
inline void rsShepardToneGenerator<T>::setBellFloor(T newFloor)
{
  if(newFloor != bellFloor) {
    bellFloor = newFloor;
    coeffsDirty = true; }
}

template<class T>
inline T rsShepardToneGenerator<T>::getGainForPitch(T pitch)
{
  T loCutoffPitch = centerPitch - T(0.5)*pitchWidth;
  T hiCutoffPitch = centerPitch + T(0.5)*pitchWidth;
  if(pitch < loCutoffPitch || pitch > hiCutoffPitch)
    return T(0);
  // This is relevant only for plotting. In getSample, we already ensure that pitch is within the
  // range at the call site -> factor out a function without the check and call that from the 
  // audio code "getGainForPitchNoCheck"  or something

  T x = pitch - centerPitch;
  return resScale * (exp(x*x * argScale) + resShift);

  // ToDo: maybe taper off the "feet" of the Gaussian by multiplying with a function that has zero
  // derivative - this helps to reduce clicks even more when sines are turned on and off

}

template<class T>
inline T rsShepardToneGenerator<T>::getSample()
{
  if(coeffsDirty)
    updateCoeffs();

  static const T tau = T(2.0*PI);   // the circle constant == circumfence/radius

  T halfWidth = T(0.5)*pitchWidth;
  T p = pitch - halfWidth;               // lowest pitch to generate
  T factor = T(1);

  T fLo = rsPitchToFreq(p);              // lowest frequency to generate
  T wLo = freqToOmega(fLo);              // lowest normalized radian frequency to generate

  /*
  if(p < centerPitch - halfWidth)
  {
    p += T(12);
    factor *= T(2);
  }
  */

  T y = T(0);                            // output accumulator
  while(p <= centerPitch + halfWidth)
  {
    y += getGainForPitch(p) * sin(factor*phase);
    factor *= T(2);
    p += T(12);
  }

  // Increment time, phase and pitch and return output:
  t += dt;
  if(t >= T(1)) 
    t -= T(1);
  phase += wLo; 
  if(phase >= tau)
    phase -= tau;
  pitch += dt*pitchDelta;
  if(pitch >= startPitch + T(12)) pitch -= T(12);  // for upward sweeps
  if(pitch <= startPitch - T(12)) pitch += T(12);  // for downward sweeps
  return y;
}
// todo: optimize using double- and half-angle formulas, maybe try to use a single loop (compute
// the lowest pitch that needs to be generated, start there and only go upward - only 
// double-angle formula needed
// -there are discontinuities after every second - i think we need separate phasors for the time
//  after which to wrap around the pitch ascension and the reference sinusoid

template<class T>
inline void rsShepardToneGenerator<T>::updateCoeffs()
{
  // To figure out the factor "a" inside the exponential, given a desired floor "b", we need to 
  // solve exp(a * x^2) = b, where x = pitchWidth/2. After finding the scaler for the argument of 
  // the exponential, we need to shift the bell curve down by b to actually hit zero at the 
  // boundaries. This in turn reduces the value at the center so we must scale the whole curve up 
  // to hit one at the center again. This leads to:
  T b = bellFloor;
  T x = T(0.5) * pitchWidth; 
  argScale = log(b) / (x*x);
  resShift = -b;
  resScale = T(1) / (T(1) - b);
  coeffsDirty = false;
}
// maybe use a normalized range, say -1...+1 and integrate into rsBellFunctions - make a 
// function: tweakedGaussian(T x, T xFloor) and 
// tweakedGaussianCoeffs(T xFloor, T* argScale, T* resScale, T* resShift)
// but using it here would increase the amount of computation to do in getSample, so maybe not



#endif