namespace rosic {
namespace Sampler {

//=================================================================================================

void AmplifierCore::setup(float volume, float pan, float width, float pos)
{
  using Mat = RAPT::rsMatrix2x2<float>;

  // Construct gain/scale matrix:
  float s = RAPT::rsDbToAmp(volume);
  Mat S(s, 0.f, 0.f, s);

  // Construct pan matrix:
  float p = (0.005f*pan) + 0.5f;  // -100..+100 -> 0..1
  Mat P(2*(1-p), 0, 0, 2*p);    // verify! maybe use panRule

  // Construct width (i.e. mid/side mixing) matrix:
  float w = 0.01f * width;      // -100..+100 -> -1.0..+1.0
  float m = 2.f - fabs(w);      // gain for mid signal (side gain is w itself)
  m = sqrt(m);                  // test..seems ok
  Mat C(1, 1, 1, -1);           // conversion R/L to M/S (or back), except for factor 1/sqrt(2)
  Mat A(0.5f*m, 0, 0, 0.5f*w);  // gains for mid/side, compensation for missing 1/sqrt(2)
  Mat W = C*A*C;                // convert L/R to M/S, apply M/S gains, convert M/S to L/R
  // ...hmmm...i think we should perhpas use a formula based on sin/cos for the M/S gains
  // -> check what we do in other M/S processors
  // rosic::StereoWidth uses 
  // RAPT::rsEqualPowerGainFactors(newRatio, &midGain, &sideGain, 0.0, 1.0);
  // but this uses a parametrization via mid/side ratio, not via width in %

  // This is the current mid- and side gain as function of w (excluding the 0.5 factor):
  //   https://www.desmos.com/calculator/iycuh1lvmg   
  //   y1 = sqrt(2-|x|) and y1 = x
  // and it's not very nicely behaved. Maybe try this instead:
  //   https://www.desmos.com/calculator/mzu474syud  
  //   y1 = s*cos(a), y2 = s*sin(a) where a = pi*x/4, s = sqrt(2)
  // it's much smoother and will repeat periodically when values go beyond the range. Maybe wrap 
  // into an rsStereoWidthToMidSideGain function. Maybe provide both variants. The sin/cos based
  // formula has the power-presverving property, iff L and R are uncorrelated. The squares of both
  // gains sum to unity, see:
  //   https://www.desmos.com/calculator/go5awdnxo7
  // could we achieve polarity inversion by going beyond +-2? At +-4, we would hear only the 
  // inverted mid-signal on both channels. Can we perhaps think about it in analogy to SVD where
  // we rotate -> scale along x,y -> rotate again. Could we perhaps view the twe pannings together
  // with the uniform scaling as a decomposition of the scaling matrix in SVD? I think, the W
  // matrix is a pure rotation (up to the constant factor of sqrt(2)). So geometrically, we do:
  // scale uniformly -> pan (non-uniform scale) -> rotate -> pan
  // where the "pan" is a combined stretch-x-squeeze-y (may let's call it "thinstretching") 
  // operation or vice versa. Maybe that should be power-preserving, too - maybe a sin/cos based 
  // rule should be used for that, too. So, leaving out the scaling for a moment, our matrix is
  // Q*W*P where Q,P are "thinstretchers" and W is a rotation (all power preserving). Maybe we can
  // equate that to the k*U*S*V of a singular value composition (where k is an additional scalar
  // scaling factor that accounts for the fact that S is in general not energy-preserving), we have
  // found a new and potentially interesting matrix decomposition?

   
  // Construct position matrix:
  float q = (0.005f*pos) + 0.5f;  // -100..+100 -> 0..1
  Mat Q(2*(1-q), 0, 0, 2*q);      // verify! maybe use panLaw

  // Combine all 4 matrices and extract coeffs:
  Mat M = Q*W*P*S;              // application order is right-to-left (scale->pan->width->pos)
  gLL = M.a; gLR = M.b; gRL = M.c; gRR = M.d;


  int dummy = 0;

  // ToDo: 
  // -Verify, if the formulas are implemented correctly, i.e. correctly resemble the behavior
  //  of reference sfz implementations. They represent my first guess, i.e. the way I would 
  //  probably do it.
  //  -Check the behavior for mono and stereo samples. Maybe we need an isStereo flag to switch 
  //   between different formulas?
  // -Optimize the calculations: obtain expressions for gLL etc. and get rid of using rsMatrix2x2,
  //  but maybe keep the matrix-based implementation somewhere as prototype for documentation 
  //  purposes
  // -Maybe add a panLaw parameter that influences the computation of the P and/or the Q matrix.
  //  ...or maybe have two seperate parameters - one for each of the matrices.
  // -Maybe provide a different parametrization where the volume is expressed as linear gain, 
  //  so we may also model polarity inversions. Or maybe realize this here via an additional 
  //  boolean flag "invertPolarity"

  // See also:
  // https://en.wikipedia.org/wiki/Matrix_decoder
}

//=================================================================================================

void FilterCore::setupCutRes(FilterCore::Type type, float w, float resoGainDb)
{
  if(type == Type::Unknown || w == 0.f)
    type = Type::Bypass;
  // When the cutoff is not defined, it defaults to zero. In this case, we switch into bypass
  // mode. We also default to bypass if mode is not set...hmm...maybe we should default to
  // lpf_12 in this case? ...but only if cutoff is nonzero...we'll see...

  using namespace RAPT;

  this->type = type;
  using FO = rsOnePoleFilter<float, float>;
  using BQ = rsBiquadDesigner;  // maybe it should have a template parameter?

  static const float s = float(1/(2*PI));
  // Preliminary to cater for the API of rsBiquadDesigner - ToDo: change API (maybe write a new 
  // class fo that and deprecate the old)

  // Compute the desired filter quality factor Q from the resonance gain in dB:
  float A = rsDbToAmp(resoGainDb);  // Raw resonance amplitude
  float Q = A;                      // This is correct for 2nd order bandpass filters...
  if(type == Type::BQ_Lowpass || type == Type::BQ_Highpass) // ..low- and highpass filters need..
    Q = rsBandwidthConverter::lowpassResoGainToQ(A);        // ..a more complicated formula

  FilterImpl& i = impl;  // as abbreviation
  switch(type)
  {
  case Type::Bypass: FO::coeffsBypass(&i.fo.b0, &i.fo.b1, &i.fo.a1); return;
    // maybe use this as default branch

  case Type::FO_Lowpass:  FO::coeffsLowpassIIT( w, &i.fo.b0, &i.fo.b1, &i.fo.a1); return;
  case Type::FO_Highpass: FO::coeffsHighpassMZT(w, &i.fo.b0, &i.fo.b1, &i.fo.a1); return;

    // This API sucks - fix it! The functions should take w for the frequency (not freq in Hz and 
    // sample-rate), their names should be much shorter (e.g. coeffsLowpassRBJ or just lowpassRBJ),
    // output params should come last and be passed as pointers.
  case Type::BQ_Lowpass:  BQ::calculateCookbookLowpassCoeffs(
    i.bqd.b0, i.bqd.b1, i.bqd.b2, i.bqd.a1, i.bqd.a2, 1.f, s*w, Q); return;
  case Type::BQ_Highpass: BQ::calculateCookbookHighpassCoeffs(
    i.bqd.b0, i.bqd.b1, i.bqd.b2, i.bqd.a1, i.bqd.a2, 1.f, s*w, Q); return;
  case Type::BQ_Bandpass_Skirt: BQ::calculateCookbookBandpassConstSkirtCoeffsViaQ(
    i.bqd.b0, i.bqd.b1, i.bqd.b2, i.bqd.a1, i.bqd.a2, 1.f, s*w, Q); return;
  case Type::BQ_Bandstop: BQ::calculateCookbookBandrejectCoeffsViaQ(
    i.bqd.b0, i.bqd.b1, i.bqd.b2, i.bqd.a1, i.bqd.a2, 1.f, s*w, Q); return;

  }
  RAPT::rsError("Unknown filter type in rsSamplerFilter::setupCutRes");

  // ToDo:
  // -Make a consistent choice for all of RAPT whether recursion coeffs of filters should have a 
  //  minus sign or not and update all code accordingly. Be careful - this change ripples through 
  //  all products - maybe introduce new names for the functions and deprecate the old ones instead
  //  of just changing their code. Make benchmarks what is faster, ask at KVR what others do.
  // -Optimize: Compute resonance related stuff only when applicable. ...but maybe we should have 
  //  a separate class for first order filters anyway to save memory
  // -Figure out, if sfz+ and other implementations also use the exact formula for Q for lowpass 
  //  and highpass filters. It's quite expensive to compute and the simple identity function that 
  //  works  perfectly for bandpass gives a reasonable approximation for low- and highpass, too, 
  //  especially as Q gets larger. At low Q, it would give a little bit of extra resonance.
  // -Figure out, if the Q = A formula should also be used for bandreject filters. At the moment,
  //  we just do it, but i'm not sure, if that's the right thing to do. It seems plausible, though.
}

void FilterCore::setupGainFreqBw(Type type, float gainDb, float w, float bw)
{
  using namespace RAPT;
  rsAssert(type == Type::BQ_Bell);
  rsAssert(bw > 0.f, "Bandwidth must be positive in setupGainFreqBw");
  rsAssert( w > 0.f, "Omega must be positive in setupGainFreqBw");
  // ToDo: allow w=0 later

  this->type = type;
  static const float s = float(1/(2*PI));
  float rawGain = rsDbToAmp(gainDb);
  FilterImpl& i = impl;  // as abbreviation
  rsBiquadDesigner::calculatePrescribedNyquistGainEqCoeffs(
    i.bqd.b0, i.bqd.b1, i.bqd.b2, i.bqd.a1, i.bqd.a2, 1.f, s*w, bw, rawGain, 1.f);


  /*
  float k = powf(2.f, 0.5f*bw);
  float Q = k / (k*k - 1.f);       // Q = 2^(bo/2) / (2^bo - 1)
  // ToDo: verify formula - if correct, move to RAPT::rsBandwidthConverter and call it like:
  //float Q = rsBandwidthConverter::octavesToQ(bw);


  using BQ = rsBiquadDesigner;

  rsBiquadDesigner::calculateCookbookPeakFilterCoeffsViaQ(
  i.bqd.b0, i.bqd.b1, i.bqd.b2, i.bqd.a1, i.bqd.a2, 1.f, s*w, Q, rawGain);
  */

  //rsError("Unknown filter type in rsSamplerFilter::setupGainFreqBw");
}

/*
void FilterCore::initCoeffs()
{

}

void FilterCore::updateCoeffs()
{

}
*/

void FilterCore::processFrame(float* L, float* R)
{
  TSig io(*L, *R);
  FilterImpl& i = impl;
  switch(type)
  {
  case Type::Bypass: break;

    // 1st order filters:
  case Type::FO_Lowpass:        io = i.fo.getSample(io); break;
  case Type::FO_Highpass:       io = i.fo.getSample(io); break;

    // Biquads:
  case Type::BQ_Lowpass:        io = i.bqd.getSample(io); break;
  case Type::BQ_Highpass:       io = i.bqd.getSample(io); break;
  case Type::BQ_Bandpass_Skirt: io = i.bqd.getSample(io); break;
  case Type::BQ_Bandstop:       io = i.bqd.getSample(io); break;
  case Type::BQ_Bell:           io = i.bqd.getSample(io); break;

  };
  *L = io.x; // Preliminary - as long as we are abusing rsVector2D for the signal
  *R = io.y;

  // ...later, we want to use a simd type and retrieve the elements like so:
  //L = io[0]; R = io[1];

  // ToDo:
  // -Organize the Type enum in such a way that we can retrieve the filter topology from it via 
  //  bitmasking such that we do not need a branch for every type but only one for every topology.
  //  This will reduce the boilerplate a lot.
}

void FilterCore::resetState()
{
  FilterImpl& i = impl;
  switch(type)
  {
  case Type::Bypass: return;

  case Type::FO_Lowpass:  i.fo.resetState();  return;
  case Type::FO_Highpass: i.fo.resetState();  return;

  case Type::BQ_Lowpass:        i.bqd.resetState(); return;
  case Type::BQ_Highpass:       i.bqd.resetState(); return;
  case Type::BQ_Bandpass_Skirt: i.bqd.resetState(); return;
  case Type::BQ_Bandstop:       i.bqd.resetState(); return;
  case Type::BQ_Bell:           i.bqd.resetState(); return;

  }
  RAPT::rsError("Unknown filter type in rsSamplerFilter::resetState");
}




}
}