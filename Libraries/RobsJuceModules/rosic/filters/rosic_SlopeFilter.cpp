//-------------------------------------------------------------------------------------------------
// construction/destruction:

SlopeFilter::SlopeFilter()
{
  sampleRate = 44100.0;
  slope      = RAPT::rsAmpToDb(sqrt(0.5));  // pinking by default: slope = -3.01 dB/oct
  reset();
  updateCoefficients();
}

SlopeFilter::~SlopeFilter()
{

}

//-------------------------------------------------------------------------------------------------
// setup:

void SlopeFilter::setSampleRate(double newSampleRate)
{
  sampleRate = newSampleRate;
  updateCoefficients();
}

void SlopeFilter::setSlope(double newSlope)
{
  slope = newSlope;
  updateCoefficients();
}

//-------------------------------------------------------------------------------------------------
// inquiry:

double SlopeFilter::getMagnitudeAt(double frequency)
{
  double w, g;
  w  = 2*PI*frequency/sampleRate;
  g  = RAPT::biquadMagnitudeAt(s1b0, s1b1, s1b2, -s1a1, -s1a2, w);
  g *= RAPT::biquadMagnitudeAt(s2b0, s2b1, s2b2, -s2a1, -s2a2, w);
  return g;
}
// todo: make sign-convention for a-coeffs consistent throughout library

//-------------------------------------------------------------------------------------------------
// others:

void SlopeFilter::reset()
{
  s1x1L = s1x2L = s1y1L = s1y2L = 0.0;  
  s2x1L = s2x2L = s2y1L = s2y2L = 0.0;
  s1x1R = s1x2R = s1y1R = s1y2R = 0.0;
  s2x1R = s2x2R = s2y1R = s2y2R = 0.0;
}

//-------------------------------------------------------------------------------------------------
// internal functions:

void SlopeFilter::updateCoefficients()
{
  static const double freq1 =  250.0;
  static const double freq2 = 4000.0;
  static const double q     =    0.40406101782088438; // bandwidth is 3 octaves

  double gainFactor = RAPT::rsDbToAmp(-2.0*slope);

  BiquadDesigner::calculateCookbookLowShelvCoeffs(s1b0, s1b1, s1b2, s1a1, s1a2, 1.0/sampleRate, freq1, q, gainFactor);
  BiquadDesigner::calculateCookbookLowShelvCoeffs(s2b0, s2b1, s2b2, s2a1, s2a2, 1.0/sampleRate, freq2, q, gainFactor);

  double normalizer = 1.0 / (gainFactor*gainFactor);
  s1b0 *= normalizer;
  s1b1 *= normalizer;
  s1b2 *= normalizer;

  // What is the rationale for the formula for the normalizer? We may need a variant that 
  // normalizes to unit gain at DC (for downward slopes) and another one that normalizes at Nyquist
  // (for upward slopes). I think, it currently normalizes the gain to unity at the center/pivot 
  // frequency which is 1000 = sqrt(250 * 4000). Verify and document this!
}
