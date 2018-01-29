#ifndef rosic_DualFilter_h
#define rosic_DualFilter_h

namespace rosic
{

class rsDualFilter
{

public:



  INLINE rsFloat64x2 getSample(rsFloat64x2 x)
  {
    if(delayConfig == DELAY_BOTH)
    {
      rsFloat64x2 y1Old = y1;
      y1 = filter1->getSample(a*x + d*y2);
      y2 = filter2->getSample(b*x + c*y1Old);
    }
    else if(delayConfig == DELAY_SECOND)
    {
      y1 = filter1->getSample(a*x + d*y2);
      y2 = filter2->getSample(b*x + c*y1);
    }
    else  // delayConfig == DELAY_FIRST
    {
      y2 = filter2->getSample(b*x + c*y1);
      y1 = filter1->getSample(a*x + d*y2);
    }
    return e*y1 + f*y2;
  }

  INLINE void getSampleFrameStereo(double* left, double* right)
  {
    rsFloat64x2 tmp(*left, *right);
    *left  = tmp[0];
    *right = tmp[1];
  }

protected:

  rsPolyModule *filter1, *filter2;
  rsFloat64x2 y1, y2;
  double a, b, c, d, e, f;

};


}

#endif