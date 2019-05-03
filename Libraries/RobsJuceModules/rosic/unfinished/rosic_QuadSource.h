#ifndef rosic_QuadSource_h
#define rosic_QuadSource_h

namespace rosic
{

class rsQuadSource
{

public:


  INLINE rsFloat64x2 getSample(const rsFloat64x2& x)
  {
    return a1*s1->getSample(x) + a2*s2->getSample(x) + a3*s3->getSample(x) + a4*s4->getSample(x);
  }

protected:

  // the sources:
  rsPolyModule *s1 = nullptr, *s2 = nullptr, *s3 = nullptr, *s4 = nullptr;

  double a1, a2, a3, a4; // amplitudes of the 4 sources (maybe factor out into VectorMixer)

};


}

#endif