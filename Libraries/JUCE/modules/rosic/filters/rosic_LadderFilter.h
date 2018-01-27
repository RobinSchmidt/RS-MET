#ifndef rosic_LadderFilter_h
#define rosic_LadderFilter_h

namespace rosic
{

/*
class rsLadderFilterMono : public RAPT::rsLadderFilter<double, double>
{

};

class rsLadderFilterStereo : public RAPT::rsLadderFilter<rsFloat64x2, double>
{

};
*/

class rsLadderFilter : public RAPT::rsLadderFilter<rsFloat64x2, rsFloat64x2>
{

public:

  // we need to replicate the interface of RAPT::rsLadderFilter but with double for the parameter
  // types and convertd to rsFloat64x2 as necessarry, also take into account stereo spread for the
  // cutoff frequency

  void setStereoSpread(double newSpread) { stereoSpread = newSpread; }
    // more to do

protected:

  double stereoSpread = 0.0;

};

}

#endif