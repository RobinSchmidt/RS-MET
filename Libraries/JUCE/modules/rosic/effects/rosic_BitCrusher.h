#ifndef rosic_BitCrusher_h
#define rosic_BitCrusher_h

//// rosic-indcludes:
//#include "../math/rosic_ElementaryFunctionsReal.h"

namespace rosic
{

  /**

  This is a digital BitCrusher/distortion unit with sample rate decimation and re-quantization 
  (the latter is better known as bitcrushing).

  \todo: introduce companding transfer curve

  */

  class BitCrusher
  {

  public:

    //---------------------------------------------------------------------------------------------
    // construction/destruction:

    /** Constructor. */
    BitCrusher();   

    /** Destructor. */
    ~BitCrusher();  

    //---------------------------------------------------------------------------------------------
    // parameter settings:

    /** Sets the decimation factor for the sample rate. */
    void setDecimationFactor(int newDecimationFactor)
    { 
      if( newDecimationFactor >= 1 ) 
        decimationFactor = newDecimationFactor; 
    }

    /** Sets the quatization interval for the signal value. */
    void setQuantizationInterval(double newInterval)
    {
      if( newInterval >= 0.0 )
        quantizationInterval = newInterval;
    }

    /** Sets the amount of the effect in percent (scales the difference between the original and 
    crushed signal by a factor). */
    void setAmount(double newAmount) { amount = newAmount; }

    //---------------------------------------------------------------------------------------------
    // inquiry:

    /** Returns the decimation factor for the sample rate. */
    int getDecimationFactor() const { return decimationFactor; }

    /** Returns the quatization interval for the signal value. */
    double getQuantizationInterval() const { return quantizationInterval; }

    /** Returns the amount of the effect in percent. */
    double getAmount() const { return amount; }

    //---------------------------------------------------------------------------------------------
    // audio processing:

    /** Processes a block of samples at a time. */
    inline void processBlock(double *inOutL, double *inOutR, double numSampleFrames);

    /** Calculates a stereo-ouput frame at a time. */
    INLINE void getSampleFrameStereo(double* inOutL, double* inOutR);

    //---------------------------------------------------------------------------------------------
    // others:

    /** Resets the internal buffers to zero. */
    void reset();

    //=============================================================================================

  protected:

    double outSampleL, outSampleR;
    double quantizationInterval, amount;
    int    decimationFactor;
    int    sampleCounter;

  };

  //-----------------------------------------------------------------------------------------------
  // from here: definitions of the functions to be inlined, i.e. all functions which are supposed 
  // to be called at audio-rate (they can't be put into the .cpp file):

  inline void BitCrusher::processBlock(double *inOutL, double *inOutR, double numSampleFrames)
  {
    for(int n=0; n<numSampleFrames; n++)
      getSampleFrameStereo(&inOutL[n], &inOutR[n]);
  }

  INLINE void BitCrusher::getSampleFrameStereo(double* inOutL,  double* inOutR)
  {
    sampleCounter++;
    if( sampleCounter >= decimationFactor )
    {
      // todo: apply compander here...

      outSampleL = RAPT::rsQuant(*inOutL, quantizationInterval);
      outSampleR = RAPT::rsQuant(*inOutR, quantizationInterval);
      sampleCounter = 0;
    }

    double diffL = outSampleL - *inOutL;
    double diffR = outSampleR - *inOutR;

    *inOutL += (0.01*amount)*diffL;
    *inOutR += (0.01*amount)*diffR;

  }

} // end namespace rosic

#endif // rosic_BitCrusher_h
