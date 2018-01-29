#ifndef rosic_NewSynth_h
#define rosic_NewSynth_h

namespace rosic
{


class rsNewSynth
{

public:




  INLINE rsFloat64x2 getSample(rsFloat64x2 x)
  {
    return filter.getSample(source.getSample(x));
  }

  INLINE void getSampleFrameStereo(double* left, double* right)
  {
    rsFloat64x2 tmp(*left, *right);
    *left  = tmp[0];
    *right = tmp[1];
  }

protected:

  rsQuadSource source;
  rsDualFilter filter;
  //rsPolyModulator polyMod;

};


}

#endif