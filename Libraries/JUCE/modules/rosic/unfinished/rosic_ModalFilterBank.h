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

  void setNumModes(int newNumModes) 
  { 
    RAPT::rsAssert(newNumModes <= maxNumModes);
    numModes = newNumModes;  
  }


  /** \name Processing */

  inline rsFloat32x4 getSample(rsFloat32x4 in)
  {
    rsFloat32x4 y(0);
    for(int i = 0; i < numModes; i++)
      y += modeFilters[i].getSample(in);
  }

  inline float getSample(float in) 
  { 
    rsFloat32x4 y = getSample(rsFloat32x4(in));
    return y.getSum();
  }

  void reset()
  {
    for(int i = 0; i < numModes; i++)
      modeFilters[i].reset();
    // sampleCounter = 0;
  }
  // maybe allow a partial reset (scale all state variables by a given number between 0 and 1)



  static const int maxNumModes = 1024; // maybe factor out into a baseclass

protected:


  rsModalFilterFloatSSE2 modeFilters[maxNumModes];
  int numModes;
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