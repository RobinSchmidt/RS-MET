#ifndef rosic_EnvelopeGenerator_h
#define rosic_EnvelopeGenerator_h

// rosic-indcludes:
#include "rosic_BreakpointModulator.h"

namespace rosic
{

  /**

  This class specializes the very general BreakpointModulator into a kind of ADSR envelope.

  */

  class EnvelopeGenerator : public BreakpointModulator
  {

  public:

    //---------------------------------------------------------------------------------------------
    // construction/destruction:
    
    /** Constructor. */
    EnvelopeGenerator();   

    //---------------------------------------------------------------------------------------------
    // parameter settings:

    /** Sets the point where the envelope starts (as raw value). */
    void setStartLevel(double newStart) { setBreakpointLevel(0, newStart); }

    /** Sets the point where the envelope starts (in dB). */
    void setStartInDecibels(double newStart) { setStartLevel(dB2amp(newStart)); }

    /** Sets the point where the envelope starts (in semitones). */
    void setStartInSemitones(double newStart) { setStartLevel(pitchOffsetToFreqFactor(newStart)); }

    /** Sets the highest point of the envelope (as raw value). */
    void setPeakLevel(double newPeak) { setBreakpointLevel(1, newPeak); }

    /** Sets the highest point of the envelope (in dB). */
    void setPeakInDecibels(double newPeak) { setPeakLevel(dB2amp(newPeak)); }

    /** Sets the highest point of the envelope (in semitones). */
    void setPeakInSemitones(double newPeak) { setPeakLevel(pitchOffsetToFreqFactor(newPeak)); }

    /** Sets the sustain level (as raw value). */
    void setSustainLevel(double newSustain) 
    { 
      setBreakpointLevel(2, newSustain); 
      setBreakpointLevel(3, newSustain); 
    }

    /** Sets the sustain level (in dB). */
    void setSustainInDecibels(double newSustain) { setSustainLevel(dB2amp(newSustain)); }

    /** Sets the sustain level (in semitones). */
    void setSustainInSemitones(double newSustain) 
    { setSustainLevel(pitchOffsetToFreqFactor(newSustain)); }

    /** Sets the end point of the envelope (as raw value). */
    void setEndLevel(double newEnd) { setBreakpointLevel(4, newEnd); }

    /** Sets the end point of the envelope (in dB). */
    void setEndInDecibels(double newEnd) { setEndLevel(dB2amp(newEnd)); }

    /** Sets the end point of the envelope (in semitones). */
    void setEndInSemitones(double newEnd) { setEndLevel(pitchOffsetToFreqFactor(newEnd)); }

    /** Sets the length of attack phase (in seconds). */
    void setAttack(double newAttackTime);    

    /** Sets the hold time (in seconds). */
    void setHold(double newHoldTime);      

    /** Sets the length of decay phase (in seconds). */
    void setDecay(double newDecayTime);     
 
    /** Sets the length of release phase (in seconds). */
    void setRelease(double newReleaseTime);   

  };

} // end namespace rosic

#endif // rosic_EnvelopeGenerator_h
