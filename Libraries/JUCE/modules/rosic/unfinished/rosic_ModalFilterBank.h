#pragma once
namespace rosic
{

/** Implements a whole bank of filters of type rsModalFilterFloatSSE2. */

class rsModalBankFloatSSE2
{

public:

  void setModeParameters(int numModes,
    double* omegas,   double* amplitudes, double* phases, 
    double* attacks1, double* attacks2,   double* attackBlends,
    double* decays1,  double* decays2,    double* decayBlends);


  inline float getSample(float in) 
  { 
    rsFloat32x4 x(in);
    rsFloat32x4 y(0);
    for(int i = 0; i < numModes; i++)
      y += modeFilters[i].getSample(x);
    return y.getSum();
  }

  void reset()
  {
    for(int i = 0; i < numModes; i++)
      modeFilters[i].reset();
    // sampleCounter = 0;
  }



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