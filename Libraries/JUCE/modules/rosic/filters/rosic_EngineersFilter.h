#ifndef rosic_EngineersFilter_h
#define rosic_EngineersFilter_h

// todo: maybe put all RAPT instantiations into a common file pair (h/cpp) which should be included
// as one of the first things in rosic.h/cpp

namespace rosic
{


/** Specialization of RAPT::rsEngineersFilter for double-precision mono signals. */

class rsEngineersFilterMono : public RAPT::rsEngineersFilter<double, double>
{

public:

  INLINE double getSample(double in) { return getSampleDirect2(in); }
  // maybe implement in RAPT::rsEngineersFilter and don't subclass here but use a typedef'd
  // explicit instantiation

};


/** Specialization of RAPT::rsEngineersFilter for double-precision stereo signals. */

class rsEngineersFilterStereo : public RAPT::rsEngineersFilter<rsFloat64x2, double>
{

public:

  INLINE void getSampleFrameStereo(double* left, double* right)
  {
    //rsFloat64x2 tmp = getSample(rsFloat64x2(*left, *right));
    rsFloat64x2 tmp = getSampleDirect2(rsFloat64x2(*left, *right));
    *left  = tmp[0];
    *right = tmp[1];
  }

};

// rename to rsEngineersFilterD2D, or maybe rename the other to rsEngineersFilterMono

}

#endif