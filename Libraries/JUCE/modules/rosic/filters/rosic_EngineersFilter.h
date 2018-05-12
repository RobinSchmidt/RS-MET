#ifndef rosic_EngineersFilter_h
#define rosic_EngineersFilter_h

namespace rosic
{

//=================================================================================================

/** 

There seems to be a bug when using elliptic and chebychev-2 bandpass/bandreject/peak filters - 
any kind of filter with stopband zeros and order-doubling - seems like the zeros don't get mapped
properly or something. It works correctly in the old rosic::EngineersFilterOld
*/

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

}

#endif