#ifndef rosic_RingModulator_h
#define rosic_RingModulator_h

//// rosic-indcludes:
//#include "../filters/rosic_EllipticSubBandFilterDirectForm.h"
//#include "../filters/rosic_NyquistBlocker.h"
//#include "../generators/rosic_SineOscillator.h"
//#include "../infrastructure/rosic_MutexLock.h"

namespace rosic
{

  //===============================================================================================
  // class RingModulator:

  /**

  This is a ringmodulator effect that ringmodulates the incoming signal with a sinusoidal 
  modulation signal.

  */

  class RingModulator
  {

  public:

    //---------------------------------------------------------------------------------------------
    // construction/destruction:
 
    /** Constructor. */
    RingModulator(); 

    /** Destructor. */
    ~RingModulator(); 

    //---------------------------------------------------------------------------------------------
    // parameter settings:

    /** Sets the sample-rate. */
    void setSampleRate(double newSampleRate);

    /** Sets the modulation frequency in Hz. */
    void setModulatorFrequency(double newFrequency);

    /** Sets the factor for feeding back the output to the input. */
    void setFeedbackFactor(double newFactor) { feedbackFactor = newFactor; }

    /** Selects whether or not an oversampling of factor 2 should be used to avoid aliasing. */
    void setAntiAliasing(bool shouldAntiAlias);

    //---------------------------------------------------------------------------------------------
    // audio processing:

    /** Calculates output sample at a time. */
    INLINE double getSample(double in);

    //---------------------------------------------------------------------------------------------
    // others:

    /** Resets the internal state. */
    void reset();

    //=============================================================================================

  protected:

    /** Sets up the oscillators according to the desired frequency shift. */
    //void setupOscillators();

    double yOld;
    double modulatorFrequency, sampleRate, feedbackFactor;
    bool   antiAlias;

    rsSubBandFilterMono upsamplingFilter, downsamplingFilter;
    SineOscillator                  sineOscillator;
    //NyquistBlocker   nyquistBlocker;
    MutexLock        mutex;

  };

  //-----------------------------------------------------------------------------------------------
  // definitions of inlined functions:

  INLINE double RingModulator::getSample(double in)
  { 
    mutex.lock();

    in += feedbackFactor * yOld;

    if( antiAlias == true )
    {
      yOld  = upsamplingFilter.getSample(2*in);     // yOld is used for temporary results as well
      yOld *= sineOscillator.getSample();
      yOld  = downsamplingFilter.getSample(yOld);
      yOld  = upsamplingFilter.getSample(0.0);
      yOld *= sineOscillator.getSample();
      yOld  = downsamplingFilter.getSample(yOld);
    }
    else
      yOld = in * sineOscillator.getSample();

    mutex.unlock();

    return yOld;
  }


  //===============================================================================================
  // class RingModulatorStereo:

  /**

  This is a stereo-version of the RingModulator class.

  */

  class RingModulatorStereo
  {

  public:

    //-----------------------------------------------------------------------------------------------
    // construction/destruction:
 
    /** Constructor. */
    RingModulatorStereo(); 

    /** Destructor. */
    ~RingModulatorStereo(); 

    //---------------------------------------------------------------------------------------------
    // parameter settings:

    /** Sets the sample-rate. */
    void setSampleRate(double newSampleRate) 
    { modulatorL.setSampleRate(newSampleRate); modulatorR.setSampleRate(newSampleRate); }

    /** Sets the frequency shift in Hz. */
    void setModulatorFrequency(double newFrequency) 
    { 
      modulatorFrequency = newFrequency; 
      modulatorL.setModulatorFrequency(modulatorFrequency+0.5f*stereoOffset); 
      modulatorR.setModulatorFrequency(modulatorFrequency-0.5f*stereoOffset); 
    }

    /** Sets a feedback factor by which the frequency shifted output sample from the previous call
    will be fed back to the input. */
    void setFeedbackFactor(double newFactor) 
    { modulatorL.setFeedbackFactor(newFactor); modulatorR.setFeedbackFactor(newFactor); }

    /** Sets the feedback in percent. */
    void setFeedbackInPercent(double newFeedback) { setFeedbackFactor(0.01*newFeedback); }

    /** Selects whether or not an oversampling of factor 2 should be used to avoid aliasing. */
    void setAntiAliasing(bool shouldAntiAlias)
    { modulatorL.setAntiAliasing(shouldAntiAlias); modulatorR.setAntiAliasing(shouldAntiAlias); }

    /** Sets the offset of the frequency shift between left and right channel in Hz. */
    void setStereoOffset(double newOffset) 
    { 
      stereoOffset = newOffset; 
      modulatorL.setModulatorFrequency(modulatorFrequency+0.5f*stereoOffset); 
      modulatorR.setModulatorFrequency(modulatorFrequency-0.5f*stereoOffset); 
    }

    /** Sets the ratio between dry and wet between 0...1. */
    void setDryWetRatio(double newDryWet) 
    { RAPT::rsEqualPowerGainFactors(newDryWet, &dry, &wet, 0.0, 1.0); }

    //---------------------------------------------------------------------------------------------
    // audio processing:

    /** Calculates one stereo output frame at a time. */
    INLINE void getSampleFrameStereo(double* inOutL, double* inOutR);

    //---------------------------------------------------------------------------------------------
    // others:

    /** Resets the internal states. */
    void reset() { modulatorL.reset(); modulatorR.reset(); }

    //=============================================================================================

  protected:

    double dry, wet;
    double modulatorFrequency, stereoOffset;

    RingModulator modulatorL, modulatorR;

  };

  //-----------------------------------------------------------------------------------------------
  // definitions of inlined functions:

  INLINE void RingModulatorStereo::getSampleFrameStereo(double* inOutL, double* inOutR)
  { 
    double tmpL = modulatorL.getSample(*inOutL);
    double tmpR = modulatorR.getSample(*inOutR);

    *inOutL     = dry*(*inOutL) + wet*tmpL;
    *inOutR     = dry*(*inOutR) + wet*tmpR;
  }

} // end namespace rosic

#endif // rosic_RingModulator_h