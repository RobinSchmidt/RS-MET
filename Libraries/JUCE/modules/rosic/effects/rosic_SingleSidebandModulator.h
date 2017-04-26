#ifndef rosic_SingleSidebandModulator_h
#define rosic_SingleSidebandModulator_h

// rosic-indcludes:
#include "../filters/rosic_EllipticSubBandFilterDirectForm.h"
#include "../filters/rosic_NyquistBlocker.h"
#include "../filters/rosic_QuadratureNetwork.h"
#include "../generators/rosic_SineOscillatorStereo.h"
#include "../infrastructure/rosic_MutexLock.h"

namespace rosic
{

  //===============================================================================================
  // class SingleSidebandModulator:

  /**

  This is a single-sideband (SSB) modulator effect that SSB-modulates the incoming signal with a 
  sinusoidal modulation signal.

  */

  class SingleSidebandModulator
  {

  public:

    //---------------------------------------------------------------------------------------------
    // construction/destruction:
 
    /** Constructor. */
    SingleSidebandModulator(); 

    /** Destructor. */
    ~SingleSidebandModulator(); 

    //---------------------------------------------------------------------------------------------
    // parameter settings:

    /** Sets the sample-rate. */
    void setSampleRate(double newSampleRate);

    /** Sets the modulation frequency in Hz. */
    void setModulatorFrequency(double newFrequency);

    /** Sets the volume level (in dB) of the upper sideband. */
    void setUpperSidebandLevel(double newLevel) { usbFactor = dB2amp(newLevel); }

    /** Sets the volume level (in dB) of the lower sideband. */
    void setLowerSidebandLevel(double newLevel) { lsbFactor = dB2amp(newLevel); }

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

    //** Acquires the mutex-lock for the sine oscillator - should be done whenever
    //void acquireLock() { mutex.lock(); }

    /** Resets the internal state. */
    void reset();

    //=============================================================================================

  protected:

    double yOld;
    double modulatorFrequency, sampleRate, feedbackFactor, usbFactor, lsbFactor;
    bool   antiAlias;

    EllipticSubBandFilterDirectForm upsamplingFilter, downsamplingFilter;
    QuadratureNetwork               quadratureNetwork;
    SineOscillatorStereo            sineOscillator;
    MutexLock                       mutex;

  };

  //-----------------------------------------------------------------------------------------------
  // definitions of inlined functions:

  INLINE double SingleSidebandModulator::getSample(double in)
  { 
    mutex.lock();

    in += feedbackFactor * yOld;

    double x, xHilbert;  // input and its Hilbert-transform
    double m, mHilbert;  // modulator and its Hilbert-transform
    double tmp1, tmp2;
    double usb, lsb;     // upper and lower sideband

    if( antiAlias == true )
    {
      /*
      yOld  = upsamplingFilter.getSample(2*in);     // yOld is used for temporary results as well
      yOld *= sineOscillator.getSample();
      yOld  = downsamplingFilter.getSample(yOld);
      yOld  = upsamplingFilter.getSample(0.0);
      yOld *= sineOscillator.getSample();
      yOld  = downsamplingFilter.getSample(yOld);
      */
      yOld = 0.0;
    }
    else
    {
      quadratureNetwork.getOutputSamplePair(in, &x, &xHilbert);
      sineOscillator.getSampleFrameStereo(&m, &mHilbert); 
      tmp1 = x*m;
      tmp2 = xHilbert*mHilbert;
      usb  = tmp1-tmp2;
      lsb  = tmp1+tmp2;
      yOld = lsbFactor*usb + usbFactor*lsb;
        // it actually should be the other way around...maybe the quadrature network shifts x by 90 
        // degrees with respect to to xHilbert instead of vice versa? or the osc is wrongly phased?
    }

    mutex.unlock();

    return yOld;
  }


  //===============================================================================================
  // class SingleSidebandModulatorStereo:

  /**

  This is a stereo-version of the SingleSidebandModulator class.

  */

  class SingleSidebandModulatorStereo
  {

  public:

    //-----------------------------------------------------------------------------------------------
    // construction/destruction:
 
    /** Constructor. */
    SingleSidebandModulatorStereo(); 

    /** Destructor. */
    ~SingleSidebandModulatorStereo(); 

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

    /** Sets the volume level (in dB) of the upper sideband. */
    void setUpperSidebandLevel(double newLevel) 
    { modulatorL.setUpperSidebandLevel(newLevel); modulatorR.setUpperSidebandLevel(newLevel); }

    /** Sets the volume level (in dB) of the lower sideband. */
    void setLowerSidebandLevel(double newLevel) 
    { modulatorL.setLowerSidebandLevel(newLevel); modulatorR.setLowerSidebandLevel(newLevel); }

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
    { equalPowerGainFactors(newDryWet, &dry, &wet, 0.0, 1.0); }

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

    SingleSidebandModulator modulatorL, modulatorR;

  };

  //-----------------------------------------------------------------------------------------------
  // definitions of inlined functions:

  INLINE void SingleSidebandModulatorStereo::getSampleFrameStereo(double* inOutL, double* inOutR)
  { 
    double tmpL = modulatorL.getSample(*inOutL);
    double tmpR = modulatorR.getSample(*inOutR);

    *inOutL     = dry*(*inOutL) + wet*tmpL;
    *inOutR     = dry*(*inOutR) + wet*tmpR;
  }

} // end namespace rosic

#endif // rosic_SingleSidebandModulator_h