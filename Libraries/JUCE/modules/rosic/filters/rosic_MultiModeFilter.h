#ifndef rosic_MultiModeFilter_h
#define rosic_MultiModeFilter_h

//// rosic-indcludes:
//#include "rosic_LadderFilter.h"
//#include "rosic_FourPoleFilter.h"

namespace rosic
{

  /**

  This class defines the user parameters for the MultiModeFilter class.

  */

  class MultiModeFilterParameters
  {

  public:

    enum MultiModeFilterClasses
    {
      NO_FILTER = 0,
      LADDER_FILTER,
      TWO_STAGE_BIQUAD
      // ALLPASS_CHAIN
      // COMB
      // HIGH_ORDER
    };

    enum multiModeFilterModes
    {
      BYPASS = 0,
      MOOGISH_LOWPASS,
      LOWPASS_6,
      LOWPASS_RBJ,
      HIGHPASS_6,
      HIGHPASS_RBJ,
      BANDPASS_RBJ,
      BANDREJECT_RBJ,
      PEAK_OR_DIP_RBJ,
      LOW_SHELV_1ST,
      LOW_SHELV_RBJ,
      HIGH_SHELV_1ST,
      HIGH_SHELV_RBJ,
      ALLPASS_1ST,
      ALLPASS_RBJ,

      MORPH_LP_BP_HP,
      MORPH_LP_PK_HP,

      NUM_FILTER_MODES
    };

    double freqNominal;       // characterisic frequency without modulations applied
    double freqByKey;         // dependency of the frequency on the key
    double freqByVel;         // dependency of the frequency on the velocity
    double resonance;         // resonance parameter in percent
    int    mode;              // mode of the filter
    int    order;             // order / number of 1st order stages
    int    filterClass;       // filter class to be used (MoogyFiler, FourPoleFilter, etc.)

    MultiModeFilterParameters()
    {
      freqNominal       = 1000.0;
      freqByKey         = 0.0;
      freqByVel         = 0.0;
      resonance         = 0.0;
      mode              = BYPASS;
      order             = 4;
      filterClass       = NO_FILTER;
      //mode              = MORPHABLE1;
    }

  };

  /**

  This class is a Moog style filter similar to the MoogFilter class. This class here, however
  mnakes several tweaks to make it more CPU-friendly and customize the characteristic sound.

  */

  class MultiModeFilter //: public PresetRememberer
  {

  public:

    //---------------------------------------------------------------------------------------------
    // construction/destruction:

    /** Constructor. */
    MultiModeFilter();

    /** Destructor. */
    ~MultiModeFilter();

    //---------------------------------------------------------------------------------------------
    // parameter settings:

    /** Sets the sample-rate for this filter. */
    void setSampleRate(double newSampleRate);

    /** Sets the mode of the filter. @see: MultiModeFilterModes. */
    void setMode(int newMode);

    /** Sets the nominal characteristic frequency for this filter. */
    void setFrequencyNominal(double newFrequency);

    /** Sets the keytracking for the filter, that is: the dependency of the characteristic
    frequency on the note key. */
    void setFrequencyByKey(double newFrequencyByKey);

    /** Sets the velocity tracking for the filter, that is: the dependency of the characteristic
    frequency on the note velocity. */
    void setFrequencyByVel(double newFrequencyByVel);

    /** Sets the current key. */
    void setKey(double newKey);

    /** Sets the current key and velocity. */
    void setKeyAndVel(double newKey, double newVel);

    /** Switches the filter into glide mode in which the filter's cutoff frequency will not
    immediately jump to new value when the note changes (and key-tracking is non-zero), but slide
    to the new target cutoff frequency. */
    //void setGlideMode(bool shouldGlide);

    /** Sets up the time it takes to glide to a new cutoff frequency after a new note was received
    (only relevant when glide is active and frequency-keytracking is nonzero). */
    //void setGlideTime(double newGlideTime);

    /** Sets the amount of resonance in percent where 100% leads to self oscillation. */
    void setResonance(double newResonance, bool updateCoefficients = true);

    /** Sets the Q-factor (applies only to the RBJ filter modes) */
    void setQ(double newQ, bool updateCoefficients = true);

    /** Sets the gain (in dB) for shelving and peaking modes. */
    void setGain(double newGain);

    /** Sets the drive for the nonlinear filters. */
    void setDrive(double newDrive);

    /** Sets the the order of the filter for the Moogish filter. */
    void setOrder(int newOrder);

    /** Decides, whether only the first stage or both should be used for the FourPoleFilter. */
    void useTwoStages(bool shouldUseTwoStages);

    /** Sets the frequency of the allpass filter for the Moog filter. */
    void setAllpassFreq(double newAllpassFreq);

    /** Sets the gain compensation for certain types of filters (in percent). */
    void setMakeUp(double newMakeUp);

    /** Sets the cutoff frequency for this filter - the actual coefficient calculation may be
    supressed by passing 'false' as second parameter, in this case, it should be triggered manually
    later by calling calculateCoefficients. */
    INLINE void setFrequencyInstantaneous(double newFrequency, bool updateCoefficients = true);

    /** Sets the morph parameter for the morphable modes. */
    INLINE void setMorph(double newMorph, bool updateCoefficients = true);

    /** Sets the morph-transition parameter for the morphable modes. */
    INLINE void setTransition(double newTransition, bool updateCoefficients = true);

    /** Causes the filter to re-calculate the coeffiecients. */
    INLINE void calculateCoefficients();

    //---------------------------------------------------------------------------------------------
    // inquiry:

    /** Returns the mode of the filter. @see: MultiModeFilterModes. */
    int getMode();

    /** Returns the nominal (unmodulated) characteristic frequency of this filter. */
    double getFrequencyNominal();

    /** Returns the keytracking for the filter, that is: the dependency of the characteristic
    frequency on the note key. */
    double getFrequencyByKey();

    /** Sets the velocity tracking for the filter, that is: the dependency of the characteristic
    frequency on the note velocity. */
    double getFrequencyByVel();

    /** Returns the characteristic frequency of this filter scaled by key and velocity.
    DEPRECATE AND REPLACE with getFreqWithKeyVelAndGlide */
    double getFrequencyWithKeyAndVel();

    /** Informs, whether or not this filter will glide it's cutoff frequency when a new note comes
    in and keytracking is nonzero. */
    //bool isInGlideMode();

    /** Sets up the time it takes to glide to a new cutoff frequency after a new note was received
    (only relevant when glide is active and frequency-keytracking is nonzero). */
    //double getGlideTime();

    /** Returns the resonance parameter of this filter. */
    double getResonance();

    /** Returns the Q-factor of this filter. */
    double getQ();

    /** Returns the gain (in dB) for shelving and peaking modes. */
    double getGain();

    /** Returns the drive for the nonlinear filters. */
    double getDrive();

    /** Returns the order of the filter for the Moogish mode. */
    int getOrder();

    /** Informs, whether only the first stage or both are used of the FourPoleFilter. */
    bool usesTwoStages();

    /** Returns the frequency of the allpass filter for the Moog filter. */
    double getAllpassFreq();

    /** Returns the make-up parameter of the Moog filter. */
    double getMakeUp();

    /** Returns the morph parameter for the morphable modes. */
    double getMorph();

    /** Informs whether or not the currently selected filter mode supports the Q-parameter. */
    bool currentModeSupportsQ();

    /** Informs whether or not the currently selected filter mode supports the gain-parameter. */
    bool currentModeSupportsGain();

    /** Informs whether or not the currently selected filter mode supports the two-stages
    parameter. */
    bool currentModeSupportsTwoStages();

    /** Returns the value of the filters transfer function (including resonance) at the complex
    value z (the filter is viewed as a linear filter for that matter, which it is not). */
    Complex getTransferFunctionAt(Complex z);

    /** Returns the value of the magnitude response of the filter at the given frequency. */
    double getMagnitudeAt(double frequency);

    /** Calculates the magnitudes of the frequency-response at the frequencies given in the array
    "frequencies" (in Hz) and stores them in the array "magnitudes". Both arrays are assumed to be
    "numBins" long. "inDecibels" indicates, if the frequency response should be returned in
    decibels. If "accumulate" is true, the magnitude response of this biquad-cascade will be
    multiplied with (or added to, when "inDecibels" is true) to the magnitudes which are already
    there in the "magnitudes"-array. This is useful for calculating the magnitude response
    of several biquad-cascades in series. */
    void getMagnitudeResponse(double *frequencies, double *magnitudes, int numBins,
      bool inDecibels = false, bool accumulate = false);

    //---------------------------------------------------------------------------------------------
    // audio processing:

    /** Calculates one output sample at a time. */
    INLINE double getSample(double in);

    /** Calculates one stereo output sample frame at a time. */
    INLINE void getSampleFrameStereo(double *inL, double *inR, double *outL, double *outR);

    /** Calculates one stereo output sample frame at a time. */
    INLINE void getSampleFrameStereo(double *inOutL, double *inOutR);

    //---------------------------------------------------------------------------------------------
    // matser/slave configuration:

    /** Adds a  slave instance to the vector of slaves - this will also cause the 'isMaster' flag
    of the slave to be set to false, redirect the slaves parameters-pointer to the one of this
    instance and delete the old (now unused) parameters pointer of the slave. */
    void addSlave(MultiModeFilter* newSlave);

    //---------------------------------------------------------------------------------------------
    // others:

    void   resetBuffers();

    //---------------------------------------------------------------------------------------------
    // public data members:

    LadderFilter ladderFilter;
    // A filter inspired by the famous Moog lowpass.

    FourPoleFilter twoStageBiquad;
    // A filter based on two stages of biquad filters.

    //=============================================================================================

  protected:

    /** Updates the characteristic frequency scaled by key and velocity taking into account the
    nominal frequency, the current key and velocity and the keyTrack and velocityTrack
    parameters.  DEPRECATE.... */
    void updateFreqWithKeyAndVel();


    MultiModeFilterParameters* parameters;
    // A pointer to the parameters which are potentially shared by among instances.

    std::vector<MultiModeFilter*> slaves;
    // A vector of pointers to other instances of this class which shall be kept in sync to this
    // instance with regard to their parameters.

    bool isMaster;
    // A flag which indicates whether or not this instance is a master which controls other
    // instances of this class - this will also determine whether or not this objects will delete
    // the pointer to the parameter set on destruction. By default, instances will be constructed
    // as master, later they can be re-configured as slaves by adding them as slaves via addSlave
    // to another instance. */

    double freqWithKeyAndVel;
    // Characteristic frequency with key and velocity tracking applied.

    double freqInstantaneous;
    // Characterisic frequency including all modulations.

    double currentKey, currentVel;
    // Current key and velocity.

  };

  //-----------------------------------------------------------------------------------------------
  // from here: definitions of the functions to be inlined, i.e. all functions
  // which are supposed to be called at audio-rate (they can't be put into
  // the .cpp file):

  INLINE void MultiModeFilter::setFrequencyInstantaneous(double newFrequency,
    bool updateCoefficients)
  {
    ladderFilter.setCutoff(newFrequency, updateCoefficients);
    twoStageBiquad.setFrequency(newFrequency);
    if( updateCoefficients == true )
      twoStageBiquad.updateFilterCoefficients();

    // \todo: optimize this - currently it seems to calcalute the coeffs for both filters even
    // though it should b done fo only one of them....nah not true....but bad design somehow
  }

  INLINE void MultiModeFilter::setMorph(double newMorph, bool updateCoefficients)
  {
    twoStageBiquad.setMorph(newMorph);
    if( updateCoefficients == true )
      twoStageBiquad.updateFilterCoefficients();
  }

  INLINE void MultiModeFilter::setTransition(double newTransition, bool updateCoefficients)
  {
    twoStageBiquad.setTransition(newTransition);
    if( updateCoefficients == true )
      twoStageBiquad.updateFilterCoefficients();
  }

  INLINE void MultiModeFilter::calculateCoefficients()
  {
    switch( parameters->filterClass )
    {
    case MultiModeFilterParameters::LADDER_FILTER:
      ladderFilter.calculateCoefficients();
      break;

    default:
      twoStageBiquad.updateFilterCoefficients();
      break;
    }
  }

  INLINE double MultiModeFilter::getSample(double in)
  {
    switch( parameters->filterClass )
    {
    case MultiModeFilterParameters::LADDER_FILTER:    return( ladderFilter.getSample(in) );
    case MultiModeFilterParameters::TWO_STAGE_BIQUAD: return( twoStageBiquad.getSample(in) );
    default:               return in;
    }
  }

  INLINE void MultiModeFilter::getSampleFrameStereo(double *inOutL, double *inOutR)
  {
    double tmpL = *inOutL;
    double tmpR = *inOutR;
    getSampleFrameStereo(&tmpL, &tmpR, inOutL, inOutR);
  }

  INLINE void MultiModeFilter::getSampleFrameStereo(double *inL, double *inR, double *outL, double *outR)
  {
    switch( parameters->filterClass )
    {
    case MultiModeFilterParameters::NO_FILTER:
      {
        *outL = *inL;
        *outR = *inR;
      }
      break;

    case MultiModeFilterParameters::LADDER_FILTER:
      {
        *outL = *inL;
        *outR = *inR;      
        ladderFilter.getSampleFrameStereo(outL, outR);
      }
      break;

    case MultiModeFilterParameters::TWO_STAGE_BIQUAD:
      {      
        *outL = *inL;
        *outR = *inR;  
        twoStageBiquad.getSampleFrameStereo(outL, outR);
      }
      break;

      //.....

    default: // bypass
      {
        *outL = *inL;
        *outR = *inR;
      }
      break;
    }
  }

}

#endif // rosic_MultiModeFilter_h
