#ifndef rosic_DualFilter_h
#define rosic_DualFilter_h

namespace rosic
{

class rsDualFilterPoly : public rsPolyModule
{

public:

  /** Determines whether the output of the 1st filter or the output of the 2nd filter (or both) 
  has a unit delay in the feedback path. */
  enum delayConfigs
  {
    DELAY_FIRST = 0,
    DELAY_SECOND,
    DELAY_BOTH
  };


  //virtual rsFloat64x2 getSample(const rsFloat64x2& in, int voice)


  rsFloat64x2 getSample(const rsFloat64x2& x, int v) override
  {
    if(delayConfig == DELAY_BOTH)
    {
      rsFloat64x2 y1Old = y1;
      y1 = filter1->getSample(a*x + d*y2,    v);
      y2 = filter2->getSample(b*x + c*y1Old, v);
    }
    else if(delayConfig == DELAY_SECOND)
    {
      y1 = filter1->getSample(a*x + d*y2, v);
      y2 = filter2->getSample(b*x + c*y1, v);
    }
    else  // delayConfig == DELAY_FIRST
    {
      y2 = filter2->getSample(b*x + c*y1, v);
      y1 = filter1->getSample(a*x + d*y2, v);
    }
    return e*y1 + f*y2;
  }
  // needs test

  /*
  INLINE void getSampleFrameStereo(double* left, double* right)
  {
    rsFloat64x2 tmp(*left, *right);
    *left  = tmp[0];
    *right = tmp[1];
  }
  */

protected:

  /** Given x and y (in the range 0..1), this function computes our a,..,f coefficients (gains for
  the inputs, feedback and outputs. ...or maybe the range should be -1..+1 */
  void computeCoeffs(double x, double y);

  int delayConfig = DELAY_BOTH;
  rsPolyModule *filter1 = nullptr, *filter2 = nullptr;

  // factor out - these need to be polyphonic, too
  rsFloat64x2 y1, y2;
  double a, b, c, d, e, f;

};


}

#endif