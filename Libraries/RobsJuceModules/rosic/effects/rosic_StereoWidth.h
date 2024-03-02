#ifndef rosic_StereoWidth_h
#define rosic_StereoWidth_h

namespace rosic
{

/** This is a simple Mid/Side mixer which can be used to widen or narrowing the stereo-image. L 
and R are converted to M and S channels, individual gain is applied to those M and S signals and 
then they are re-converted to L and R. A global gain is applied afterwards. So it has 3 Parameters:
midGain, sideGain and globalGain.

\todo: extend this class with a ChannelRemixer, parameters:
gain, pan, mid/side ratio, inverter switches for left, right, mid side  */

class StereoWidth
{

public:

  //-----------------------------------------------------------------------------------------------
  // \name Lifetime

  /** Constructor. */
  StereoWidth();

  /** Destructor. */
  ~StereoWidth();


  //-----------------------------------------------------------------------------------------------
  // \name Setup

  /** Adjusts the ratio between the mid and side signal where 0.0 means mid only and 1.0 means
  side only. */
  void setMidSideRatio(double newRatio)
  {
    RAPT::rsEqualPowerGainFactors(newRatio, &midGain, &sideGain, 0.0, 1.0);
  }


  /** Adjusts the gain for the mid-signal - value is expected in dB. */
  //void setMidGain(double newMidGain) { midGain = RAPT::rsDbToAmp(newMidGain); }

  /** Adjusts the gain for the side-signal - value is expected in dB. */
  //void setSideGain(double newSideGain) { sideGain = RAPT::rsDbToAmp(newSideGain); }

  // The setMidGain/setSideGain functions are deprecated and are replaced by setMidSideRatio. It 
  // doesn't make sense to keep them because using them together with setMidSideRatio will 
  // interfere. Either the user should use setMidGain/setSideGain *or* setMidSideRatio - but not a
  // mix of them. So we enforce this here to reduce possible misuse. There doesn't seem to be any
  // code that uses setMidGain/setSideGain anyway.


  /** Adjusts the global gain for the signal - value is expected in dB. */
  void setGlobalGain(double newGlobalGain) { globalGain = SQRT2_INV*RAPT::rsDbToAmp(newGlobalGain); }

  /** Selects, that the final output should be mixed to mono. Useful to quickly check mono
  compatibility. */
  void setMixToMono(bool shouldMixToMono) { mixToMono = shouldMixToMono; }

  /** Selects whether the output should be polarity inverted. This can be useful in mixing to 
  quickly check, which polarity sounds better. */
  void setInvertPolarity(bool shouldInvert) { invertPolarity = shouldInvert; }


  //-----------------------------------------------------------------------------------------------
  // \name Processing

  /** Calculates a stereo-ouput frame. */
  INLINE void getSampleFrameStereo(double* inOutL, double* inOutR);


  //-----------------------------------------------------------------------------------------------
  // \name Internals

protected:

  doubleA midGain        = 1.0;     // The user sets up the gain values in dB
  doubleA sideGain       = 1.0;     // ...but we keep them here as raw factors
  doubleA globalGain     = 1.0;
  bool    mixToMono      = false;
  bool    invertPolarity = false;

};


INLINE void StereoWidth::getSampleFrameStereo(double* inOutL, double* inOutR)
{
  double mid  = midGain    * (*inOutL + *inOutR);
  double side = sideGain   * (*inOutL - *inOutR);
  *inOutL     = globalGain * (mid + side);
  *inOutR     = globalGain * (mid - side);

  if(mixToMono == true)
    * inOutL = *inOutR = SQRT2_INV * (*inOutL + *inOutR);

  if(invertPolarity == true)
  {
    *inOutL = -(*inOutL);
    *inOutR = -(*inOutR);
  }

  // ToDo:
  // -Optimize - everything can be collapsed into matrix-vector multiplication:
  //    outL = gain_LL * inL  +  gain_LR * inR;
  //    outR = gain_RL * inL  +  gain_RR * inR;
}

} // end namespace rosic

#endif // rosic_StereoWidth_h
