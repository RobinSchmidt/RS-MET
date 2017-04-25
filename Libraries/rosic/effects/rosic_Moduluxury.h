#ifndef rosic_Moduluxury_h
#define rosic_Moduluxury_h

// rosic-indcludes:
//#include "../filters/rosic_TwoStageBiquad.h"
#include "../filters/rosic_MultiModeFilter.h"
#include "../modulators/rosic_Modulator.h"
#include "../delaylines/rosic_FractionalDelayLineStereo.h"

namespace rosic
{

  /**

  This is a effect that can apply tremolo, vibrato, wah-wah or pan-modulation to
  an incoming signal. Alternatively it can output its modulating signal for 
  control-purposes.

  // todo: don't derive from Modulator but embedd two Modulators for left and right channel

  */

  class Moduluxury : public Modulator
  {

  public:

    enum modulationTargets
    {
      NO_TARGET = 0,
      AMPLITUDE,
      PANORAMA,
      FILTER_FREQUENCY,
      DELAY,
      //PITCH,
      MIDI_CONTROLLER
    };

    enum panRules
    {
      LINEAR_PAN = 1,
      CONSTANT_POWER_PAN
    };

    enum delayModulationModes
    {
      ADDITIVE = 1,
      MULTIPLICATIVE
    };

    //---------------------------------------------------------------------------------------------
    // construction/destruction:

    Moduluxury();   
    ///< Constructor.

    ~Moduluxury();  
    ///< Destructor.

    //---------------------------------------------------------------------------------------------
    // parameter settings (set-functions):

    void setSampleRate(double newSampleRate);
    ///< Sets up the sample-rate.

    void setModulationTarget(int newModulationTarget);
    /**< Sets the target-parameter which will be modulated. 
         @see modulationTargets. */

    void setPanRule(int newPanRule);
    /**< Sets the panning rule for panorama modulation. 
         @see panRules. */

    void setFilterFrequency(double newFilterFrequency);
    ///< Sets the nominal (not yet modulated) cutoff-/center-frequency for the filter.

    void setDelayModulationMode(int newDelayModulationMode);
    /**< Selects the mode for delay modulation. 
         @see delayModulationModes */

    void setDelay(double newDelay);
    ///< Sets the nominal (not yet modulated) delay for the delayline (in seconds).

    void setMidiControllerNumber(int newControllerNumber);
    /**< Sets the controller number which should be modulated. */

    //void setScaleFactor(double newScaleFactor);
    /**< Sets a scale factor for the modulating signal. */

    //void setOffset(double newOffset);
    /**< Sets an offset for the modulating signal. */

    //---------------------------------------------------------------------------------------------
    // inquiry (get-, is-, etc. functions):

    int getModulationTarget();
    /**< Returns the currently chosen the target-parameter which will be modulated. 
         @see modulationTargets. */

    int getPanRule();
    /**< Returns the panning rule for panorama modulation. 
         @see panRules. */

    double getFilterFrequency();
    ///< Returns the nominal (not yet modulated) cutoff-/center-frequency for the filter.

    int getDelayModulationMode();
    /**< Returns the mode for delay modulation. 
         @see delayModulationModes */

    double getDelay();
    ///< Returns the nominal (not yet modulated) delay for the delayline (in seconds).

    int getMidiControllerNumber();
    /**< Returns the controller number which is modulated. */

    //double getScaleFactor();
    /**< Returns the scale factor for the modulating signal. */

    //double getOffset();
    /**< Returns the offset for the modulating signal. */

    //---------------------------------------------------------------------------------------------
    // audio processing:

    INLINE void getSampleFrameStereo(double* inL,  double* inR, 
      double* outL, double* outR, int* midiControllerOut);
    /**< Calculates a stereo-ouput frame. */

    //---------------------------------------------------------------------------------------------
    // embedded public modules:
    FractionalDelayLineStereo delayLine;
    //TwoStageBiquad filter;
    MultiModeFilter filter;

    //Modulator modulationGenerator;

    //=============================================================================================

  protected:

    int    modulationTarget;
    int    panRule;
    double filterFrequency;
    int    delayModulationMode;
    double delay;
    int    midiControlNumber;
    //double scaleFactor;
    //double offset;

  };

  //-----------------------------------------------------------------------------------------------
  // from here: definitions of the functions to be inlined, i.e. all functions
  // which are supposed to be called at audio-rate (they can't be put into
  // the .cpp file):

  INLINE void Moduluxury::getSampleFrameStereo(double* inL,  double* inR, 
    double* outL, double* outR, int* midiControllerOut)
  {
    //*modValue = BreakpointModulator::getSample();

    // get the raw modulator output:
    double mod;
    if( modulationSource == BREAKPOINT_MODULATOR )
      mod = breakpointModulator.getSample();
    else
      mod = sampleModulator.getSample();

    // apply scaling and offset:
    mod = scaleFactor * mod + offset;

    // apply the SlewRateLimiter:
    mod        = upDownSlewRateLimiter.getSample(mod);

    // calculate the instantaneous depth:
    double instDepth = 0.0;
    if( fadeInIsOn && currentKey != -1 ) 
      instDepth = fadeInOutSlewRateLimiter.getSample(depth);  // fade in the depth:
    else if( fadeOutIsOn && currentKey == -1 )
      instDepth = fadeInOutSlewRateLimiter.getSample(0.0);    // fade out the depth:
    else 
      instDepth = depth;
    


    // initialize the midi control change output to a 'neutral' value:
    *midiControllerOut = 64;



    switch( modulationTarget )
    {
    case AMPLITUDE:
      {
        //*outL = (*inL) * (*modValue);
        //*outR = (*inR) * (*modValue);
        *outL = instDepth * (*inL) * mod + (1.0-instDepth) * (*inL);
        *outR = instDepth * (*inR) * mod + (1.0-instDepth) * (*inR);
      }
      break;

    case PANORAMA:
      {
        // map the range -1...+1 to 0...1:
        //*modValue = 0.5 * (*modValue + 1.0);
        mod = 0.5 * (instDepth * mod + 1.0);

        switch( panRule )
        {
        case LINEAR_PAN:
          {
            *outL = (*inL) * mod;
            *outR = (*inR) * (1.0-mod);
          }
          break;
        case CONSTANT_POWER_PAN:
          {
            double s, c;
            sinCosApprox(0.5 * PI * mod, &s, &c);
            *outL = (*inL) * s;
            *outR = (*inR) * c;
          }
          break;
        default: // pan-rule setting is messed up - output silence
          {
            *outL = 0.0;
            *outR = 0.0; 
          }
        }  // end of  switch( panRule )
      }
      break;

    case FILTER_FREQUENCY:
      {
        // calculate the instantaneous filter frequency
        double f = instDepth * filterFrequency * mod + (1.0-instDepth) * filterFrequency; 
        //double f = filterFrequency * (*modValue); // instantaneous filter frequency
        // setup and apply the filter:
        filter.setFrequencyInstantaneous(f, true);
        //filter.updateFilterCoefficients(); // is done by passing 'true' to the former call
        filter.getSampleFrameStereo(inL, inR, outL, outR);
      }
      break;

    case DELAY:
      {
        if( delayModulationMode == ADDITIVE )
          //delayLine.setDelayTime( delay + *modValue * delay );
          delayLine.setDelayTime( delay + mod * delay * instDepth );
        else
          //delayLine.setDelayTime( delay * (*modValue) );
          delayLine.setDelayTime( instDepth * delay * mod + (1.0-instDepth) * delay );

        *outL = *inL;
        *outR = *inR;
        delayLine.getSampleFrameStereo(outL, outR);
      }
      break;

    case MIDI_CONTROLLER:
      {
        //*modValue = *modValue * scaleFactor + offset;
        //*modValue = round(*modValue);
        mod = mod * scaleFactor * instDepth + offset;
        *midiControllerOut = roundToInt(127.0*mod);
        if( *midiControllerOut < 0 )
          *midiControllerOut = 0;
        else if( *midiControllerOut > 127 )
          *midiControllerOut = 127;
      }
      break;

    } // end of switch( modulationTarget )
  }

} // end namespace rosic

#endif // rosic_Moduluxury_h
