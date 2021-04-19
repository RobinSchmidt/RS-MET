#ifndef rosic_LadderFilter_h
#define rosic_LadderFilter_h

namespace rosic
{

/*
class rsLadderFilterMono : public RAPT::rsLadderFilter<double, double>
{

};

*/

//=================================================================================================

/** A ladder filter based on instantiating the template ladder filter with the "rsFloat64x2" SSE2 
wrapper for the signal and "double" for the parameters (cutoff, etc). This gives two parallel
ladder filters with the same settings for the price of one. */

class rsLadderFilterStereo : public RAPT::rsLadderFilter<rsFloat64x2, double>
{

public:

  INLINE void getSampleFrameStereo(double* left, double* right)
  {
    rsFloat64x2 tmp = getSample(rsFloat64x2(*left, *right));
    *left  = tmp[0];
    *right = tmp[1];
  }

};

//=================================================================================================

/** A ladder filter based on instantiating the template ladder filter with the rsFloat64x2 SSE2 
wrapper for both, the signal and the parameters (cutoff, etc.). This gives two parallel ladder
filters with potentially different settings for the price of one, as far as the audio processing 
functions (getSample, etc.) are concerned. Setting paramaters (setCutoff, etc.) will still cost as 
much as if there would be two separate filters because the coefficient compuations can't be 
vectorized (yet?) and require to fall back to scalar code. */

/*
class rsLadderFilter : public RAPT::rsLadderFilter<rsFloat64x2, rsFloat64x2> // maybe rename
{

public:

  // we need to replicate the interface of RAPT::rsLadderFilter but with double for the parameter
  // types and convertd to rsFloat64x2 as necessarry, also take into account stereo spread for the
  // cutoff frequency
  // maybe do not inherit publically but protected

  void setStereoSpread(double newSpread) { stereoSpread = newSpread; }
    // more to do

protected:

  double stereoSpread = 0.0;

};
*/





}

#endif