#ifndef rosic_MusiciansFilter_h // maybe rename to ResoWave - also the files
#define rosic_MusiciansFilter_h

namespace rosic
{


class rsResoWave // mayb rename to rsResoWaveFilter
{

public:

  //-----------------------------------------------------------------------------------------------
  /** \name Lifetime */

  rsResoWave();

  //-----------------------------------------------------------------------------------------------
  /** \name Setup */

  void setCutoff(float newCutoff);
  // maybe use setOmega instead?

  void setResonance(float newResonance);

  // ...etc...




  //-----------------------------------------------------------------------------------------------
  /** \name Processing */


protected:

  RAPT::rsLadderFilter<rsFloat32x4, rsFloat32x4> ladder;
  // ToDo: 
  // -use rsSimdVector<float, 4> instead of rsFloat32x4 when the class is ready for production
  // -factor out a class rsLadderFilterCore that is set of in terms of omega - get rid of the 
  //  sample-rate
  // -the 4 elements of the simd vector are used for:: 
  //  0: resonant left, 1: resonant right 2: non-resonant left, 3: non-resonant right
};





} // rosic



#endif