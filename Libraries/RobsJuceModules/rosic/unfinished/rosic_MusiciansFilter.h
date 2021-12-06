#ifndef rosic_MusiciansFilter_h // maybe rename to ResoWave - also the files
#define rosic_MusiciansFilter_h

namespace rosic
{


/** A filter based on 2 instances of RAPT::rsLadderFilter (concptually - in reality there's just one 
instance using simd for creating the 2 signal paths), one with resonance and one without. A pure 
resonance signal is obtained by subtracting the nonresonant output from the resonant. The 
purpose of this splitting of the output into the 2 components "filtered input" and "resonance" is 
to make the resonance signal availabe for post-processing in order to achieve filters with 
different character. The post processing itself is delegated to subclasses or driver classes that
embed an object of class rsResoSplitFilter. ...tbc... */

class rsResoSplitFilter
{

public:

  //-----------------------------------------------------------------------------------------------
  /** \name Setup */

  void setup(float newOmega, float newResonance);


  //-----------------------------------------------------------------------------------------------
  /** \name Processing */

  void processFrame(float inL, float inR, 
    float* fltOutL, float* fltOutR, float* resOutL, float* resOutR);



protected:

  RAPT::rsLadderFilter<rsFloat32x4, rsFloat32x4> ladder;
  // ToDo: 
  // -use rsSimdVector<float, 4> instead of rsFloat32x4 when the class is ready for production
  // -factor out a class rsLadderFilterCore that is set of in terms of omega - get rid of the 
  //  sample-rate
  // -the 4 elements of the simd vector are used for:: 
  //  0: resonant left, 1: resonant right 2: non-resonant left, 3: non-resonant right
  // -hmm...actually, we could perhaps optimize the memory usage by using indeed 2 instances but
  //  with smaller datatypes and a special nonresonant implementation..we'll see

};

//=================================================================================================

class rsResoWaveFilter // mayb rename to rsResoWaveFilter
{

public:

  //-----------------------------------------------------------------------------------------------
  /** \name Lifetime */

  rsResoWaveFilter();

  //-----------------------------------------------------------------------------------------------
  /** \name Setup */

  void setCutoff(float newCutoff);
  // maybe use setOmega instead?

  void setResonance(float newResonance);
  // maybe set up the resonance in terms of a decay time in ms

  // ...etc...




  //-----------------------------------------------------------------------------------------------
  /** \name Processing */


protected:

  rsResoSplitFilter resoSplitter;

};





} // rosic



#endif