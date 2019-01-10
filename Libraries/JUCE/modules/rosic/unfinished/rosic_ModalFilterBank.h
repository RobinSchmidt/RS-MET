#pragma once
namespace rosic
{

/** Implements a whole bank of filters of type rsModalFilterFloatSSE2. */

class rsModalBankFloatSSE2
{

public:

  /** \name Setup */

  //void setModeParameters(int numModes,
  //  double* omegas,   double* amplitudes, double* phases, 
  //  double* attacks1, double* attacks2,   double* attackBlends,
  //  double* decays1,  double* decays2,    double* decayBlends);

  /** Sets the parameters for the filter with given index. */
  void setModalFilterParameters(int index,
    double omega, double amplitude, double phase, double attack, double decay,
    double deltaOmega = 0, double phaseDelta = 0, double blend = 0.5,
    double attackScale = 1.0, double decayScale = 1.0)
  {
    modeFilters[index].setParameters(omega, amplitude, phase, attack, decay,
      deltaOmega, phaseDelta, blend, attackScale, decayScale);
  }

  /*
  void setNumModes(int newNumModes) 
  { 
    RAPT::rsAssert(newNumModes <= maxNumModes);
    numModes = newNumModes;  
  }
  */

  /** Index of the lowest mode to be produced (range: 0..maxNumModes-1). */
  void setLowModeIndex(int newIndex)
  {
    RAPT::rsAssert(newIndex >= 0 && newIndex < maxNumModes);
    lowModeIndex = newIndex;
  }


  /** Index of the highest mode to be produced (range: 0..maxNumModes-1). */
  void setHighModeIndex(int newIndex)
  {
    RAPT::rsAssert(newIndex >= 0 && newIndex < maxNumModes);
    highModeIndex = newIndex;
  }


  /** \name Inquiry */

  //int getTotalLength(double threshold = 0.000001); // cut off at -120 dB
  // see getRingingTimeEstimate in RAPT::rsBiquadCascade for hwo to compute this - but it may be
  // better to compute it in rsModalSynth where we have the decay-times available - it would be
  // wasteful to re-compute the decay time from the filter-coeffs



  /** \name Processing */

  inline rsFloat32x4 getSample(rsFloat32x4 in)
  {
    rsFloat32x4 y(0);
    for(int i = lowModeIndex; i <= highModeIndex; i++)
      y += modeFilters[i].getSample(in);
    return y;
  }

  inline float getSample(float in) 
  { 
    rsFloat32x4 y = getSample(rsFloat32x4(in));
    return y.getSum();
  }

  // maybe the getSample functions should be made virtual in a baseclass rsModalBank, such 
  // that we can use a pointer to a baseclass in rsModalSynth and instantiate either a SSE2 or
  // AVX2 subclass based on the available CPU. the function call overhead will be dwarfed by
  // the loop over the modes anyway

  void reset()
  {
    for(int i = 0; i < maxNumModes; i++)
      modeFilters[i].reset();
    // sampleCounter = 0;
  }
  // maybe allow a partial reset (scale all state variables by a given number between 0 and 1)

  // maybe have a function resetActiveModes that only runs the loop from lowModeIndex to
  // highModeIndex - might be more efficient - but actually it's not something to be called at 
  // sample-rate anyway, so who cares?


  static const int maxNumModes = 1024; // maybe factor out into a baseclass

protected:

  rsModalFilterFloatSSE2 modeFilters[maxNumModes];
  int lowModeIndex  = 0;
  int highModeIndex = maxNumModes-1;

  //int numModes;
  // int sampleCounter = 0;
  // maybe need a system to switch off modes that have gone silent - maybe sort the modes by their
  // decayTimes (descending) and at each sample (or block) decrease the upper loop limit for the
  // loop over the modes in getSample. but we would also have to re-activate them when we get a new
  // nonzero input/excitation - the whole thing would make sense only for impulse excitation - for
  // continuous excitation (like noise), all modes must remain active all the time
};




/*
class rsModalBank
{

protected:

  static const int maxNumVoices = 16;
  static const int maxNumModes  = 1024;

  rsModalFilterFloatSSE2 modeFilters[maxNumVoices][maxNumModes];

  struct ModeParameters
  {
    float frequency, amplitude, phase, attack1, attack2, attackBlend, decay1, decay2, decayBlend;
  };
  ModeParameters modeParams[maxNumModes];


  double sampleRate;
  int numModes;


  //std::vector<rsModalFilterFloatSSE2> modeFilters;

};
*/

}