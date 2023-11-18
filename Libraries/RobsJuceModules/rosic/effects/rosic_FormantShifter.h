#ifndef rosic_FormantShifter_h
#define rosic_FormantShifter_h

//// rosic-indcludes:
//#include "../others/rosic_SpectralEnvelopeProcessor.h"

namespace rosic
{

  //===============================================================================================
  // class FormantShifter:

  /**

  This class implements a spectrum based formant shifter. It can scale the spectral envelope by
  a factor and shift the spectral envelope by some constant offset in Hz.

  ToDo:
  -Make this class realtime ready by not allocating in processSpectrum.


  */

  class FormantShifter : public SpectralEnvelopeProcessor
  {

  public:

    //---------------------------------------------------------------------------------------------
    // construction/destruction:

    /** Constructor. You must pass the maximum blocksize here and may optionally pass a maximum
    overlap- and zero-padding factor. For spectral processing with a cosine^2 window on the input,
    an overlap and zero-padding of 2 is usually a good choice. */
    FormantShifter(int maxBlockSize, int maxOverlapFactor = 2,
      int maxPaddingFactor = 2);

    /** Destructor */
    virtual ~FormantShifter() {}

    //---------------------------------------------------------------------------------------------
    // parameter settings:

    /** Sets a scale factor for the formant frequencies. */
    void setFormantScale(double newScale) { scale = newScale; }

    /** Sets a frequency offset (in Hz) for the formant frequencies. */
    void setFormantOffset(double newOffset) { offset = newOffset; }

    //=============================================================================================

  protected:

    /** This function does the actual filtering of the spectrum. */
    virtual void processSpectrum(Complex *spectrum, int spectrumSize);

    double scale, offset;
    bool   energyCompensation;

  };

  //===============================================================================================
  // class FormantShifterStereo:

  /**

  This is a stereo-version of the FormantShifter class.

  */

  class FormantShifterStereo
  {

  public:

    //-----------------------------------------------------------------------------------------------
    // construction/destruction:

    /** Constructor. */
    FormantShifterStereo(int maxBlockSize, int maxOverlapFactor = 2, int maxPaddingFactor = 2);

    /** Destructor. */
    ~FormantShifterStereo();

    //---------------------------------------------------------------------------------------------
    // parameter settings:

    /** Sets the sample-rate. */
    void setSampleRate(double newSampleRate)
    { shifterL.setSampleRate(newSampleRate); shifterR.setSampleRate(newSampleRate); }

    /** Sets a scale factor for the formant frequencies. */
    void setFormantScale(double newScale)
    { shifterL.setFormantScale(newScale); shifterR.setFormantScale(newScale); }

    /** Sets a frequency offset (in Hz) for the formant frequencies. */
    void setFormantOffset(double newOffset)
    { shifterL.setFormantOffset(newOffset);  shifterR.setFormantOffset(newOffset); }

    /** Sets the ratio between dry and wet between 0...1. */
    void setDryWetRatio(double newDryWet)
    { RAPT::rsEqualPowerGainFactors(newDryWet, &dry, &wet, 0.0, 1.0); }

    // setBlockSize, setMono

    //---------------------------------------------------------------------------------------------
    // audio processing:

    /** Calculates one stereo output frame at a time. */
    INLINE void getSampleFrameStereo(double* inOutL, double* inOutR);

    //---------------------------------------------------------------------------------------------
    // others:

    /** Resets the internal states. */
    void reset() { shifterL.clearBuffers(); shifterR.clearBuffers(); }

    //=============================================================================================

  protected:

    double dry, wet;
    FormantShifter shifterL, shifterR;

  };

  //-----------------------------------------------------------------------------------------------
  // definitions of inlined functions:

  INLINE void FormantShifterStereo::getSampleFrameStereo(double* inOutL, double* inOutR)
  {
    double tmpL = shifterL.getSample(*inOutL);
    double tmpR = shifterR.getSample(*inOutR);
    *inOutL     = dry*(*inOutL) + wet*tmpL;
    *inOutR     = dry*(*inOutR) + wet*tmpR;
  }

} // end namespace rosic

#endif // rosic_FormantShifter_h
