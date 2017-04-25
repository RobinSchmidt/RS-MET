#ifndef rosic_Tremolo_h
#define rosic_Tremolo_h

// rosic-indcludes:
#include "rosic_ModulationEffect.h"

namespace rosic
{

  /**

  This class implements a tremolo effect with sinusoidal modulators for left and right channel
  separately such that it can also generate an auto-pan effect.

  */

  class Tremolo : public ModulationEffect
  {

  public:

    //---------------------------------------------------------------------------------------------
    // audio processing:

    INLINE void getSampleFrameStereo(double *inOutL, double *inOutR);

    //=============================================================================================

  };

  //-----------------------------------------------------------------------------------------------
  // from here: definitions of the functions to be inlined, i.e. all functions which are supposed
  // to be called at audio-rate (they can't be put into the .cpp file):

  INLINE void Tremolo::getSampleFrameStereo(double *inOutL, double *inOutR)
  {
    double left, right;
    lfo.getSampleFrameStereo(&left, &right);
    *inOutL *= (1.0 + depth * left);
    *inOutR *= (1.0 + depth * right);
  }

} // end namespace rosic

#endif // #ifndef rosic_Tremolo_h
