#ifndef rosic_ModulationEffect_h
#define rosic_ModulationEffect_h

// rosic-indcludes:
#include "../modulators/rosic_LowFrequencyOscillator.h"

namespace rosic
{

  /**

  This class serves a subclass for simple modulation effects like tremolo, vibrato, etc. It 
  provides two sinusoidal LFOs - one for the left and one for the right channel.

  todo: make a bigger version of this class with waveform loading (left and right channels 
  separately), AR-slew-rate limiters on the LFOs,
  let an AR-envelope follower control the depth - vibrato only for sustain, etc...

  */

  class ModulationEffect 
  {

  public:

    //---------------------------------------------------------------------------------------------
    // construction/destruction:

    /** Constructor. */
    ModulationEffect();

    /** Destructor */
    ~ModulationEffect();

    //---------------------------------------------------------------------------------------------
    // parameter settings (set-functions):

    /** Sets the sample-rate. */
    void setSampleRate(double newSampleRate);

    /** Selects one of the standard waveforms. @see WaveTable waveforms. */
    //void setWaveform(int newWaveform) { lfo.setWaveform(newWaveform); }

    /** Passes a custom waveform to be used for the LFO. */
    //void setWaveform(double *newWaveformBuffer, int newLength, char* name)
    //{ lfo.setWaveform(newWaveformBuffer, newLength, name); }

    /** Passes a custom waveform to be used for the LFO. */
    //void setWaveform(float *newWaveformBuffer, int newLength, char* name)
    //{ lfo.setWaveform(newWaveformBuffer, newLength, name); }

    /** Sets the cycle-length in seconds or beats (depending on whether sync is active). */
    void setCycleLength(double newCycleLength);

    /** Sets the depth of the effect bewteen 0...1, where 1 is full depth. */
    void setDepth(double newDepth) { depth = newDepth; }

    /** Sets the depth of the effect in percent. */
    void setDepthInPercent(double newDepthInPercent) { setDepth(newDepthInPercent); }

    /** Sets the start phase of the LFOs. That the phase value to which the LFOs will be reset 
    on calls to triggerLFOs. */
    void setStartPhaseInDegrees(double newPhase) 
    { lfo.setStartPhase(newPhase); }

    /** Sets the attack-time (in milliseconds) for the slewrate limiter - this is time which it 
    takes to rise 63% for upward jumps in the waveform. */
    void setUpwardSlewRate(double newSlewRate) 
    { lfo.setUpwardSlewRate(newSlewRate); }

    /** Sets the release-time (in milliseconds) for the slewrate limiter - this is time which it 
    takes to fall to 37% for downward jumps in the waveform. */
    void setDownwardSlewRate(double newSlewRate) 
    { lfo.setDownwardSlewRate(newSlewRate); }

    /** Sets the phase offset between left and right channel LFO. Left channels phase offset will 
    be shifted downwards by half this value, right channels phase offset will be shifted upwards
    by half this value. */
    void setStereoPhaseOffsetInDegrees(double newOffset) { stereoPhaseOffset = newOffset; }

    /** Switches the tempo-sync on or off. */
    void setTempoSync(bool shouldTempoSync);

    /** Sets up the tempo in  beats per minute. */
    void setTempoInBPM(double newTempoInBPM);

    //---------------------------------------------------------------------------------------------
    // inquiry:

    /** Returns the cycle-length in seconds or beats (depending on whether sync is active). */
    double getCycleLength() const { return lfo.getCycleLength(); }

    /** Returns true when tempo-sync is active, false otherwise. */
    int isInSyncMode() const { return lfo.isInSyncMode(); }

    /** Returns the name of the waveform as zero-terminated c-string - if an audiofile is used, 
    this 'name' will represent the relative path of the file. */
    //char* getWaveformName() const { return lfo.getWaveformName(); }

    /** Fills the targetBuffer with values suitable for displaying the current waveform in a
    display. The targetBuffer is assumed to be of size [numSamplesToShow] where the index
    represents the sample-number on the display (which is the x-coordinate in pixels). */
    //void getWaveformForDisplay(double* targetBuffer, int numSamplesToShow) 
    //{ lfo.getWaveformForDisplay(targetBuffer, numSamplesToShow); }

    //---------------------------------------------------------------------------------------------
    // others:

    /** Resets the phases of the oscillators to their start phases with an optional additional 
    phase-shift (given in radiant). */
    void resetOscillatorPhases(double phaseShift = 0.0);

    //---------------------------------------------------------------------------------------------
    // embedded modules:

    LowFrequencyOscillator lfo;

    //=============================================================================================

  protected:

    double sampleRate;
    double bpm;
    double depth;               // modulation depth
    double stereoPhaseOffset;   // phase-offset between left and right channel LFO in degrees

  };

} // end namespace rosic

#endif // #ifndef rosic_ModulationEffect_h
