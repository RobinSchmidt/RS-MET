#ifndef rosic_Constants_h
#define rosic_Constants_h

namespace rosic
{

  // This file defines various constants and enumerations that are used throughout the rosic library
 
  /** The classical standard waveforms that are found in synthesizers. */
  enum standardWaveforms
  {
    SINE = 0,
    SAW,
    SQUARE,
    TRIANGLE,

    NUM_STANDARD_WAVEFORMS
  };
  // make enum class, rename to rsStandardWaveforms ...or just rsWaveForms

  // \todo: move mathematical constants here too (PI, SQRT2, etc.)


} // end namespace rosic

#endif 