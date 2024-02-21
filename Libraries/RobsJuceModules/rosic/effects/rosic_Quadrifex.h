#ifndef rosic_Quadrifex_h
#define rosic_Quadrifex_h

//// rosic-indcludes:
//#include "../infrastructure/rosic_EffectModules.h"
//#include "../others/rosic_RoutingMatrix.h"

namespace rosic
{

  /**

  This class implements a multieffect with 4 assignable effects slots and various routings between 
  these 4 slots and global feedback.

  */

  class Quadrifex
  {

  public:

    /** An enumeration of the different effect algorithms. */
    enum effectAlgorithms
    {
      MUTE = 0,
      BYPASS,

      BIT_CRUSHER,
      CHORUS,
      COMB_BANK, // still to do....
      COMB_RESONATOR,
      COMB_STEREOIZER,
      COMPRESSOR,
      //COMP_SHAPER,
      DUAL_TWO_POLE_FILTER,
      EQUALIZER,
      EXPANDER,
      FLANGER,
      FORMANT_SHIFTER,
      FOUR_POLE_FILTER,
      FREQUENCY_SHIFTER,
      HARMONICS,
      LADDER_FILTER,
      LIMITER,
      MODULATED_ALLPASS,
      NOISE_GATE,
      NOISIFIER,
      //ONE_POLE_FILTER,
      PHASER,
      PHASE_STEREOIZER,
      PINGPONG_ECHO,
      PITCH_SHIFTER,
      REVERB,
      RINGMODULATOR,
      SIMPLE_DELAY,
      SINE_OSCILLATOR,
      SSB_MODULATOR,
      SLEWRATE_LIMITER,
      SLOPE_FILTER,
      STEREO_PAN,
      STEREO_WIDTH,
      TREMOLO,
      TWO_POLE_FILTER,
      VIBRATO,
      WAH_WAH,
      WAVESHAPER,

      NUM_EFFECT_ALGORITHMS
    };
    // ToDo: use an enum class

    enum slotRoutings
    {
      R_BYPASS = 0,
      R_1TO2TO3TO4,
      R_1TO2TO3_PLUS4,
      R_1TO2_PLUS_3TO4,
      R_1PLUS2PLUS3PLUS4,
      R_1PLUS2PLUS3_TO_4,
      R_1_TO_2PLUS3_TO_4,
      R_1PLUS2_TO_3TO4,
      R_1TO2_TO_3PLUS4,
      R_1PLUS2_TO_3PLUS4,
      R_1_TO_2PLUS3PLUS4,
      MATRIX,

      NUM_SLOT_ROUTINGS
    };
    // ToDo: use an enum class

    //---------------------------------------------------------------------------------------------
    // construction/destruction:

    /** Constructor. */
    Quadrifex();

    /** Destructor. */
    ~Quadrifex();

    //---------------------------------------------------------------------------------------------
    // setup:

    /** Sets up the sample-rate. */
    void setSampleRate(double newSampleRate);

    /** Sets up the tempo in beats per minute (relevant for tempo sync). */
    void setTempoInBPM(double newTempoInBPM);

    /** Chooses one of the routings between the slots @see: slotRoutings. */
    void setSlotRouting(int newRoutingIndex);

    /** Sets up the given slot to the given algorithm @see effectAlgorithms. */
    void setEffectAlgorithm(int slotIndex, int newAlgorithmIndex);

    /** Sets up the dry/wet ratio (between 0...1). */
    void setDryWet(double newDryWet);

    /** Sets a level (in dB) for the wet signal. */
    void setWetLevel(double newLevel);

    //---------------------------------------------------------------------------------------------
    // inquiry:

    /** Returns the selected routing index for the routing between the slots. @see slotRoutings. */
    int getSlotRouting() const { return slotRouting; }

    /** Returns the index of the effect algorithm that is used in the effect-slot. 
    @see: effectAlgorithms. */
    int getEffectAlgorithmIndex(int slotIndex) const;

    /** Retruns the dry/wet ratio (between 0...1). */
    double getDryWet() const { return dryWet; }

    /** Retruns the level (in dB) for the wet signal. */
    double getWetLevel() const { return wetLevel; }

    //---------------------------------------------------------------------------------------------
    // others:

    /** Returns a pointer to the EffectModule that sits in the given slotIndex. 
    WARNING: this is not thread-safe - it will cause access violations to use this pointer when 
    some other thread changes the effectModule via setEffectAlgorithm */
    Module* getEffectModule(int slotIndex) { return effectModules[slotIndex]; }

    /** Acquires the mutex lock for accessing the pointer-array to the EffectModules. */
    void acquireLock() { mutex.lock(); }

    /** Releases the mutex lock for accessing the pointer-array to the EffectModules. */
    void releaseLock() { mutex.unlock(); }

    /** Resets the internal states of all the embedded effects. */
    void reset();

    /** Triggers LFOs etc. in all effects where applicable. */
    void trigger();

    //---------------------------------------------------------------------------------------------
    // audio processing:

    /** Calculates a stereo output sampleframe. */
    //void getSampleFrameStereo(double *inOutL, double *inOutR);

    /** Calculates an entire block of samples at once. */
    void processBlock(double* inOutL, double* inOutR, int numFrames);

    //---------------------------------------------------------------------------------------------
    // the embedded modules:

    static const int numEffectSlots = 4;
    Module           *effectModules[numEffectSlots];

    RoutingMatrix mixMatrix;

    //=============================================================================================

  protected:

    double inputsL[5],  inputsR[5];
    double outputsL[5], outputsR[5];

    double sampleRate, bpm;                
    double dryWet, wetLevel, dryFactor, wetFactor;

    int slotRouting;
    int effectAlgorithmIndices[numEffectSlots];
      // these integers hold the index of the effect-algorithm  (see enum effectAlgorithms)
      // for each of the four effect slots

    MutexLock mutex; // for accessing the pointer-array 'effectModules' 
                     // ...obsolete..now handled in jura?

  };

} // end namespace rosic

#endif // rosic_Quadrifex_h
