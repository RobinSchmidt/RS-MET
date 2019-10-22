#ifndef rosic_StereoWidth_h
#define rosic_StereoWidth_h

//// rosic-indcludes:
//#include "../math/rosic_ElementaryFunctionsReal.h"

namespace rosic
{

  /**

  This is a simple Mid/Side mixer which can be used to widen or narrowing the stereo-image. L and R 
  are converted to M and S channels, individual gain is applied to those M and S signals and then 
  they are re-converted to L and R. A global gain is applied afterwards. So it has 3 Parameters: 
  midGain, sideGain and globalGain.

  \todo: extend this class with a ChannelRemixer, parameters: 
  gain, pan, mid/side ratio, inverter switches for left, right, mid side

  */

  class StereoWidth
  {

  public:

    //---------------------------------------------------------------------------------------------
    // construction/destruction:

    /** Constructor. */
    StereoWidth(); 

    /** Destructor. */
    ~StereoWidth();  

    //---------------------------------------------------------------------------------------------
    // parameter settings:

    /** Adjusts the ratio between the mid and side signal where 0.0 means mid only and 1.0 means 
    side only. */
    void setMidSideRatio(double newRatio) 
    { RAPT::rsEqualPowerGainFactors(newRatio, &midGain, &sideGain, 0.0, 1.0); }

    /** Adjusts the gain for the mid-signal - value is expected in dB. */
    void setMidGain(double newMidGain) { midGain = RAPT::rsDbToAmp(newMidGain); }

    /** Adjusts the gain for the side-signal - value is expected in dB. */
    void setSideGain(double newSideGain) { sideGain = RAPT::rsDbToAmp(newSideGain); }

    /** Adjusts the global gain for the signal - value is expected in dB. */
    void setGlobalGain(double newGlobalGain) { globalGain = SQRT2_INV*RAPT::rsDbToAmp(newGlobalGain); }

    /** Selects, that the final output should be mixed to mono to quickly check mono 
    compatibility. */
    void setMixToMono(bool shouldMixToMono) { mixToMono = shouldMixToMono; }

    //---------------------------------------------------------------------------------------------
    // audio processing:

    /** Calculates a stereo-ouput frame. */
    INLINE void getSampleFrameStereo(double* inOutL, double* inOutR);

    //=============================================================================================

  protected:

    // gain values (as raw factors):
    doubleA midGain;
    doubleA sideGain;
    doubleA globalGain;
    bool    mixToMono;

  };

  //-----------------------------------------------------------------------------------------------
  // from here: definitions of the functions to be inlined, i.e. all functions which are supposed 
  // to be called at audio-rate (they can't be put into the .cpp file):

  INLINE void StereoWidth::getSampleFrameStereo(double* inOutL,  double* inOutR)
  {
    double mid  = midGain    * ( *inOutL + *inOutR );
    double side = sideGain   * ( *inOutL - *inOutR ); 
    *inOutL     = globalGain * (  mid + side );
    *inOutR     = globalGain * (  mid - side );

    if( mixToMono == true )
      *inOutL = *inOutR = SQRT2_INV * (*inOutL + *inOutR);

    // \todo: optimize (maybe)
  }

} // end namespace rosic

#endif // rosic_StereoWidth_h
