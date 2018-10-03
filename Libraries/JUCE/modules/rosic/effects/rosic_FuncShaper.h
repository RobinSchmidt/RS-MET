#ifndef rosic_FuncShaper_h
#define rosic_FuncShaper_h

//// rosic-indcludes:
//#include "../basics/rosic_TabulatedFunction.h"
//#include "../filters/rosic_LowpassHighpass.h"
//#include "../filters/rosic_EllipticSubBandFilter.h"
//#include "../filters/rosic_OnePoleFilterStereo.h"

namespace rosic
{

  /**

  This is a waveshaping  distortion unit which lets the user specify the mapping between the input 
  and output as a mathematical formula.

  */

  class FuncShaper ///: public PresetRememberer
  {

  public:

    //-----------------------------------------------------------------------------------------------
    // construction/destruction:
 
    /** Constructor. */
    FuncShaper(); 

    /** Destructor. */
    ~FuncShaper(); 

    //-----------------------------------------------------------------------------------------------
    // parameter settings:

    /** Sets the sample-rate for this module */
    void setSampleRate(double newSampleRate);

    /** Sets a new function string in the embedded TabulatedFunction object and returns true, if 
    the string is valid. Optionally, it promptly re-calculates the table but sometimes this is 
    undesirable because the not all variables are properly set up yet. If false is passed as second 
    argument, you can avoid prompt re-calculation - in this case you should call calculateTable 
    yorself manually at some later stage. */
    bool setFunctionString(const char* newFunctionString, bool reCalculateTable);

    /** Assigns the a-constant in the expression for the distortion curve to a numeric value */
    void setA(double newA, bool reCalculateTable);

    /** Sets the minimum value for the a-parameter (this is only for the GUI). */
    //void setMinA(double newMinA);

    /** Sets the maximum value for the a-parameter (this is only for the GUI). */
    //void setMaxA(double newMaxA);

    /** Assigns the b-constant in the expression for the distortion curve to a numeric value */
    void setB(double newB, bool reCalculateTable);

    /** Sets the minimum value for the b-parameter (this is only for the GUI). */
    //void setMinB(double newMinB);

    /** Sets the maximum value for the b-parameter (this is only for the GUI). */
    //void setMaxB(double newMaxB);

    /** Assigns the c-constant in the expression for the distortion curve to a numeric value */
    void setC(double newC, bool reCalculateTable);

    /** Sets the minimum value for the c-parameter (this is only for the GUI). */
    //void setMinC(double newMinC);

    /** Sets the maximum value for the c-parameter (this is only for the GUI). */
    //void setMaxC(double newMaxC);

    /** Assigns the d-constant in the expression for the distortion curve to a numeric value */
    void setD(double newD, bool reCalculateTable);

    /** Sets the minimum value for the d-parameter (this is only for the GUI). */
    //void setMinD(double newMinD);

    /** Sets the maximum value for the d-parameter (this is only for the GUI). */
    //void setMaxD(double newMaxD);

    /** Switches the input filter on or off. */
    void useInputFilter(bool shouldBeUsed);

    /** Sets the cutoff frequency of the input lowpass filter (expected in Hz). */
    void setInLowpassCutoff(double newInLowpassCutoff);

    /** Sets the cutoff frequency of the input highpass filter (expected in Hz). */
    void setInHighpassCutoff(double newInHighpassCutoff);

    /** Sets the drive (gain) for the input (expected in dB). */
    void setDrive(double newDrive);

    /** Sets the DC offset to be added to the signal (expected as raw value). */
    void setDcOffset(double newDcOffset);

    /** Sets the oversampling factor for the internal signal processing - a factor of one (or 
    lower) will result in no oversampling. */
    void setOversampling(int newOversamplingFactor);

    /** Switches the output filter on or off. */
    void useOutputFilter(bool shouldBeUsed);

    /** Sets the cutoff frequency of the output lowpass filter (expected in Hz). */
    void setOutLowpassCutoff(double newOutLowpassCutoff);

    /** Sets the cutoff frequency of the output highpass filter (expected in Hz). */
    void setOutHighpassCutoff(double newOutHighpassCutoff);

    /** Sets the global output volume (expected in dB). */
    void setOutVol(double newOutVol);

    /** Sets the ratio between the original and the effect signal (expected in % wet). */
    void setDryWet(double newDryWet);

    //---------------------------------------------------------------------------------------------
    // inquiry:

    /** Returns the current function-string as c-string. */
    const char* getFunctionString() const { return distortionCurve.getExpressionString(); }

    /** Returns the value of the a-parameter. */
    //double getA() const { return a; }

    /** Returns the minimum value of the a-parameter. */
    //double getMinA() const { return minA; }

    /** Returns the maximum value of the a-parameter. */
    //double getMaxA() const { return maxA; }

    /** Returns the value of the b-parameter. */
    //double getB() const { return b; }

    /** Returns the minimum value of the b-parameter. */
    //double getMinB() const { return minB; }

    /** Returns the maximum value of the b-parameter. */
    //double getMaxB() const { return maxB; }

    /** Returns the value of the c-parameter. */
    //double getC() const { return c; }

    /** Returns the minimum value of the c-parameter. */
    //double getMinC() const { return minC; }

    /** Returns the maximum value of the c-parameter. */
    //double getMaxC() const { return maxC; }

    /** Returns the value of the d-parameter. */
    //double getD() const { return d; }

    /** Returns the minimum value of the d-parameter. */
    //double getMinD() const { return minD; }

    /** Returns the maximum value of the d-parameter. */
    //double getMaxD() const { return maxD; }

    /** Returns true, when the input filter is active, false otherwise. */
    //bool isInputFilterUsed() const { return inFilterActive; }

    /** Returns the cutoff frequency of the input lowpass filter. */
    //double getInLowpassCutoff() const { return inputFilterL.getLowpassCutoff(); }

    /** Returns the cutoff frequency of the input highpass filter. */
    //double getInHighpassCutoff() const { return inputFilterL.getHighpassCutoff(); }

    /** Returns the value of the drive parameter in dB. */
    //double getDrive() const { return drive; }

    /** Returns the value of the dc-offset parameter. */
    //double getDcOffset() const { return dcOffset; }

    /** Returns the oversampling factor. */
    //int getOversampling() const { return oversampling; }

    /** Returns true, when the output filter is active, false otherwise. */
    //bool isOutputFilterUsed() const { return outFilterActive; }

    /** Returns the cutoff frequency of the output lowpass filter. */
    //double getOutLowpassCutoff() const { return outputFilterL.getLowpassCutoff(); }

    /** Returns the cutoff frequency of the output highpass filter. */
    //double getOutHighpassCutoff() const { return outputFilterL.getHighpassCutoff(); }

    /** Returns the output volume of the wet signal in dB. */
    //double getOutVol() const { return outVol; }

    /** Returns the value of the dry/wet parameter in % wet. */
    //double getDryWet() const { return dryWet; }

    //---------------------------------------------------------------------------------------------
    // audio processing:

    /** Calculates one output sample at a time. */
    INLINE double getSample(double in);

    /** Calculates one streo output frame at a time. */
    INLINE void getSampleFrameStereo(double* inL, double* inR, double* outL, double* outR);

    //---------------------------------------------------------------------------------------------
    // others:

    /** Triggers (re-) calculation of the function table. */
    void calculateTable();

    //---------------------------------------------------------------------------------------------
    // embedded objects (in the order of the signal flow):

    LowpassHighpass           inputFilterL, inputFilterR;
    TabulatedFunction         distortionCurve;
    LowpassHighpass           outputFilterL, outputFilterR;
    OnePoleFilterStereo       deClickingFilter;
    rsSubBandFilterMonoBQ upsamplerL, upsamplerR;
    rsSubBandFilterMonoBQ antiAliasFilterL, antiAliasFilterR;

    //=============================================================================================

  protected:

    // parameter variables:
    doubleA a, b, c, d;
    doubleA driveFactor, dcOffset;
    doubleA outVolFactor;
    doubleA dryVol, wetVol;

    doubleA sampleRate;

    doubleA minA, maxA, minB, maxB, minC, maxC, minD, maxD;

    doubleA drive, outVol, dryWet;

    bool inFilterActive, outFilterActive;
    int  oversampling;

    // for smoother preset switching (reduce clicks for presets with different DC):
    int numFadeSamples;
    int fadeCountDown;
  };

  //-----------------------------------------------------------------------------------------------
  // from here: definitions of the functions to be inlined, i.e. all functions which are supposed 
  // to be called at audio-rate (they can't be put into the .cpp file):

  INLINE double FuncShaper::getSample(double in)
  {
    return distortionCurve.getValueLinear(in);
  }

  INLINE void FuncShaper::getSampleFrameStereo(double* inL,  double* inR, 
    double* outL, double* outR)
  {
    static doubleA tmpL, tmpR, oldOutL, oldOutR;
    static intA i;

    // apply the drive, input filter and dc-offset:
    if( inFilterActive )
    {
      tmpL = driveFactor * inputFilterL.getSample(*inL) + dcOffset;
      tmpR = driveFactor * inputFilterR.getSample(*inR) + dcOffset;
    }
    else
    {
      tmpL = driveFactor * (*inL) + dcOffset;
      tmpR = driveFactor * (*inR) + dcOffset;
    }

    distortionCurve.acquireLock();
    if( oversampling > 1 )
    {
      // calculate n frames of the oversampled, distorted, and anti-alias filtered signal
      // (we do oversampling by a factor of n):
      tmpL = upsamplerL.getSampleDirect1(oversampling*tmpL);
      tmpR = upsamplerR.getSampleDirect1(oversampling*tmpR);

      tmpL = distortionCurve.getValueLinearNoLock(tmpL); 
      tmpR = distortionCurve.getValueLinearNoLock(tmpR);

      tmpL = antiAliasFilterL.getSampleDirect1(tmpL);
      tmpR = antiAliasFilterR.getSampleDirect1(tmpR);

      for(i=1; i<oversampling; i++)
      {
        tmpL = upsamplerL.getSampleDirect1(0.0);
        tmpR = upsamplerR.getSampleDirect1(0.0);

        tmpL = distortionCurve.getValueLinearNoLock(tmpL);
        tmpR = distortionCurve.getValueLinearNoLock(tmpR);

        tmpL = antiAliasFilterL.getSampleDirect1(tmpL);
        tmpR = antiAliasFilterR.getSampleDirect1(tmpR);
      }
    }
    else
    {
      tmpL = distortionCurve.getValueLinearNoLock(tmpL);
      tmpR = distortionCurve.getValueLinearNoLock(tmpR);
    }
    distortionCurve.releaseLock();

    // apply the output filter:
    if( outFilterActive )
    {
      tmpL = outputFilterL.getSample(tmpL);
      tmpR = outputFilterR.getSample(tmpR);
    }

    // mix dry and wet signal:
    *outL = outVolFactor*wetVol*tmpL + dryVol*(*inL);
    *outR = outVolFactor*wetVol*tmpR + dryVol*(*inR);

    // apply a kind of crossfade/smoothing if the preset was switched to avoid clicks:
    if( fadeCountDown > 0 )
    {
      // warm up the deClickingFilter:
      if( fadeCountDown == numFadeSamples )
        deClickingFilter.setState(0.0, 0.0, oldOutL, oldOutR);

      deClickingFilter.getSampleFrameStereo(outL, outR, outL, outR);
      fadeCountDown--;
    }

    oldOutL = *outL;
    oldOutR = *outR;
  }

} // end namespace rosic

#endif // rosic_FuncShaper_h