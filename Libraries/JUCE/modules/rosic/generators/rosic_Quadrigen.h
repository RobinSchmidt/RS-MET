#ifndef rosic_Quadrigen_h
#define rosic_Quadrigen_h

// rosic-indcludes:
#include "../infrastructure/rosic_GeneratorModules.h"
#include "../infrastructure/rosic_EffectModules.h"
#include "../others/rosic_RoutingMatrix.h" // maybe move to basics
#include "../modulators/rosic_BreakpointModulator.h"

namespace rosic
{

  /**

  This class implements a soundgenerator with 4 assignable slots for actual sound-generation
  algorithms.

  */

  class Quadrigen
  {

  public:

    /** An enumeration of the different soundgenerator algorithms. */
    enum generatorAlgorithms
    {
      MUTE = 0,
      OSCILLATOR_STEREO,

      NUM_GENERATOR_ALGORITHMS
    };

    /** An enumeration of the different modulator algorithms. */
    enum modulatorAlgorithms
    {
      ALWAYS_ZERO = 0,
      BREAKPOINT_MODULATOR,

      NUM_MODULATOR_ALGORITHMS
    };

    //---------------------------------------------------------------------------------------------
    // construction/destruction:

    /** Constructor. */
    Quadrigen();

    /** Destructor. */
    ~Quadrigen();

    //---------------------------------------------------------------------------------------------
    // setup:

    /** Sets up the sample-rate. */
    void setSampleRate(double newSampleRate);

    /** Sets up the tempo in beats per minute (relevant for tempo sync). */
    void setTempoInBPM(double newTempoInBPM);

    /** Sets up the given slot to the given algorithm @see generatorAlgorithms. */
    void setGeneratorAlgorithm(int slotIndex, int newAlgorithmIndex);

    /** Sets up the given slot to the given algorithm @see modulatorAlgorithms. */
    void setModulatorAlgorithm(int slotIndex, int newAlgorithmIndex);

    //---------------------------------------------------------------------------------------------
    // inquiry:

    /** Returns the index of the generator algorithm that is used in the generator-slot. 
    @see: generatorAlgorithms. */
    int getGeneratorAlgorithmIndex(int slotIndex) const;

    /** Returns the index of the modulator algorithm that is used in the modulator-slot. 
    @see: modulatorAlgorithms. */
    int getModulatorAlgorithmIndex(int slotIndex) const;

    //---------------------------------------------------------------------------------------------
    // others:

    /** Returns a pointer to the GeneratorModule that sits in the given slotIndex. 
    WARNING: this is not thread-safe - it will cause access violations to use this pointer when 
    some other thread changes the generatorModule via setGeneratorAlgorithm */
    Module* getGeneratorModule(int slotIndex) { return generatorModules[slotIndex]; }

    /** Returns a pointer to the ModulationSource that sits in the given slotIndex. 
    WARNING: this is not thread-safe - it will cause access violations to use this pointer when 
    some other thread changes the modulatorModule via setModulatorAlgorithm */
    ModulationSource* getModulatorModule(int slotIndex) { return modulatorModules[slotIndex]; }

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
    INLINE void getSampleFrameStereo(double *inOutL, double *inOutR);

    /** Calculates an entire block of samples at once. */
    void processBlock(float *inOutL, float *inOutR, int numFrames);

    //---------------------------------------------------------------------------------------------
    // the embedded modules:

    static const int maxNumVoices      = 32;
    static const int numGeneratorSlots = 4;
    static const int numModulatorSlots = 8;
    static const int numEffectSlots    = 2;    // for the filters

    Module           *generatorModules[numGeneratorSlots];    // ...times maxNumVoices
    ModulationSource *modulatorModules[numModulatorSlots];    // ...times maxNumVoices
    //Module           *effectModules[numEffectSlots];       // ...times maxNumVoices

    // ...4x4 RM-Matrix

    RoutingMatrix mixMatrix;

    //=============================================================================================

  protected:

    double inputsL[5],  inputsR[5];
    double outputsL[5], outputsR[5];

    double sampleRate, bpm;                
    double dryWet, wetLevel, dryFactor, wetFactor;

    int slotRouting;
    int generatorAlgorithmIndices[numGeneratorSlots];
    int modulatorAlgorithmIndices[numModulatorSlots];

    MutexLock mutex; 
      // for accessing the pointer-array 'generatorModules' and/or 'modulatorModules'

  };

  //-----------------------------------------------------------------------------------------------
  // from here: definitions of the functions to be inlined, i.e. all functions which are supposed 
  // to be called at audio-rate (they can't be put into the .cpp file):

  INLINE void Quadrigen::getSampleFrameStereo(double *inOutL, double *inOutR)
  {
    *inOutL = *inOutR = 0.0;
  }

} // end namespace rosic

#endif // rosic_Quadrigen_h
