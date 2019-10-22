#ifndef rosic_EqualizerStereo_h
#define rosic_EqualizerStereo_h

//// rosic-indcludes:
//#include "rosic_Equalizer.h"

namespace rosic
{

  /**

  This class puts together two rosic::Equalizer objects in order to implement different equalization settings for two channels. It provides
  the follwing channel modes: stereo linked, left/right, mid/side, mono where the stereo linked corresponds to the operation of the 
  embedded Equalizer class.
  
  */

  class EqualizerStereo
  {

  public:

    enum stereoModes
    {
      STEREO_LINKED = 0,
      LEFT_RIGHT,
      MID_SIDE,
      MONO
    };

    //-------------------------------------------------------------------------------------------------------------------------------------
    // construction/destruction:

    /** Constructor. */
    EqualizerStereo();

    /** Destructor. */
    ~EqualizerStereo();

    //-------------------------------------------------------------------------------------------------------------------------------------
    // setup:

    /** Switches the whole equalizer into (or out of) bypass mode. */
    void setBypass(bool shouldBeBypassed) 
    { 
      bypass = shouldBeBypassed; 
    }

    /** Selects one of the available modes for stereo processing. @see stereoModes */
    void setStereoMode(int newStereoMode) 
    { 
      stereoMode = newStereoMode; 
      if( stereoMode == MONO )
        equalizers[0].setMono(true);
      else
        equalizers[0].setMono(false);
    }

    void setSampleRate(double newSampleRate)
    { 
      equalizers[0].setSampleRate(newSampleRate); 
      equalizers[1].setSampleRate(newSampleRate); 
    }

    int addBand(int channel, int newMode, double newFrequency, double newGain, double newBandwidth = 2.0*asinh(1.0/sqrt(2.0))/log(2.0))
    { 
      return equalizers[channel].addBand(newMode, newFrequency, newGain, newBandwidth); 
    }

    bool removeBand(int channel, unsigned int index)
    { 
      return equalizers[channel].removeBand(index);  
    }

    bool removeAllBands()
    {
      bool result  = equalizers[0].removeAllBands();
      result      &= equalizers[1].removeAllBands();
      return result;
    }

    bool setBandMode(int channel, int index, int newMode)   
    { 
      return equalizers[channel].setBandMode(index, newMode); 
    }

    bool setBandFrequency(int channel, int index, double newFrequency)  
    { 
      return equalizers[channel].setBandFrequency(index, newFrequency); 
    }

    bool setBandGain(int channel, int index, double newGain)
    { 
      return equalizers[channel].setBandGain(index, newGain); 
    }

    bool setBandBandwidth(int channel, int index, double newBandwidth)
    { 
      return equalizers[channel].setBandBandwidth(index, newBandwidth); 
    }

    bool setLowerBandedgeFrequency(int channel, unsigned int index, double newLowerBandedgeFrequency)
    { 
      return equalizers[channel].setLowerBandedgeFrequency(index, newLowerBandedgeFrequency); 
    }

    bool setUpperBandedgeFrequency(int channel, unsigned int index, double newUpperBandedgeFrequency)
    { 
      return equalizers[channel].setUpperBandedgeFrequency(index, newUpperBandedgeFrequency); 
    }

    void setGlobalGain(double newGainInDB) 
    { 
      equalizers[0].setGlobalGain(newGainInDB);
      equalizers[1].setGlobalGain(newGainInDB);
    }

    //-------------------------------------------------------------------------------------------------------------------------------------
    // inquiry:

    /** Returns the selected mode for stereo processing. @see stereoModes */
    int getStereoMode()
    {
      return stereoMode;
    }

    int getNumBands(int channel)
    {
      return equalizers[channel].getNumBands();
    }

    int getBandMode(int channel, unsigned int index)
    {
      return equalizers[channel].getBandMode(index);
    }

    double getBandFrequency(int channel, unsigned int index)
    {
      return equalizers[channel].getBandFrequency(index);
    }

    double getBandGain(int channel, unsigned int index)
    {
      return equalizers[channel].getBandGain(index);
    }

    double getBandBandwidth(int channel, unsigned int index)
    {
      return equalizers[channel].getBandBandwidth(index);
    }

    double getLowerBandedgeFrequency(int channel, unsigned int index)
    {
      return equalizers[channel].getLowerBandedgeFrequency(index);
    }

    double getUpperBandedgeFrequency(int channel, unsigned int index)
    {
      return equalizers[channel].getUpperBandedgeFrequency(index);
    }

    bool doesModeSupportBandwidth(int channel, unsigned int index)
    {
      return equalizers[channel].doesModeSupportBandwidth(index);
    }

    bool doesModeSupportGain(int channel, unsigned int index)
    {
      return equalizers[channel].doesModeSupportGain(index);
    }

    double getGlobalGain() const 
    {
      return equalizers[0].getGlobalGain();
    }

    void getMagnitudeResponse(int channel, double* frequencies, double* magnitudes, int numBins);

    //-------------------------------------------------------------------------------------------------------------------------------------
    // audio processing:

    /** Calculates a mono output sample - we need this for EchoLab. */
    INLINE double getSample(double in);

    /** Calculates a stereo output sampleframe. */
    INLINE void getSampleFrameStereo(double *inOutL, double *inOutR);

    /** Calculates an entire block of samples at once. */
    INLINE void processBlock(float *inOutL, float *inOutR, int numSamples);

    //-------------------------------------------------------------------------------------------------------------------------------------
    // others:

    /** Resets the internal states of the filters. */
    void reset();

    //=====================================================================================================================================

    Equalizer equalizers[2]; // temporarily moved to public during development - need access in EchoLabStateToXml etc.

  protected:

    int  stereoMode;
    bool bypass;

  };

  //---------------------------------------------------------------------------------------------------------------------------------------
  // inlined functions:

  INLINE double EqualizerStereo::getSample(double in)
  {
    if( bypass == true )
      return in;

    return equalizers[0].getSample(in);
  }

  INLINE void EqualizerStereo::getSampleFrameStereo(double *inOutL, double *inOutR)
  {
    if( bypass == true )
      return;

    if( stereoMode == LEFT_RIGHT )
    {
      *inOutL = equalizers[0].getSample(*inOutL);
      *inOutR = equalizers[1].getSample(*inOutR);
    }
    else  // stereoMode == MID_SIDE
    {
      double mid  = SQRT2_INV * (*inOutL + *inOutR);
      double side = SQRT2_INV * (*inOutL - *inOutR);
      mid         = equalizers[0].getSample(mid);
      side        = equalizers[1].getSample(side);
      *inOutL     = SQRT2_INV * (mid + side);
      *inOutR     = SQRT2_INV * (mid - side);
    }
  }

  INLINE void EqualizerStereo::processBlock(float *inOutL, float *inOutR, int numSamples)
  {
    if( bypass == true )
      return;

    if( stereoMode == STEREO_LINKED || stereoMode == MONO )
    {
      equalizers[0].processBlock(inOutL, inOutR, numSamples);
      return;
    }

    double tmpL, tmpR;
    int n;

    for(n=0; n<numSamples; n++)
    {
      tmpL = (double) inOutL[n];
      tmpR = (double) inOutR[n];
      getSampleFrameStereo(&tmpL, &tmpR);
      inOutL[n] = (float) tmpL;
      inOutR[n] = (float) tmpR;
    }
  }

} 

#endif
