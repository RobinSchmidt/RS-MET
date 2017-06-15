//#include "rosic_TwoPoleFilter.h"
//using namespace rosic;

//-----------------------------------------------------------------------------------------------------------------------------------------
// construction/destruction:

TwoPoleFilter::TwoPoleFilter()
{
  sampleRate = 44100.0;
  mode       = PEAK;
  frequency  = 1000.0;
  gain       = 0.0;
  bandwidth  = 1.0;
  radius     = 0.9;
  updateCoeffs();
}

TwoPoleFilter::~TwoPoleFilter()
{

}

//-----------------------------------------------------------------------------------------------------------------------------------------
// setup:

void TwoPoleFilter::setSampleRate(double newSampleRate)
{
  if( newSampleRate <= 0.0 )
  {
    DEBUG_BREAK;
    return;
  }
  sampleRate = newSampleRate;
  updateCoeffs();
}

void TwoPoleFilter::setMode(int newMode)
{
  if( newMode >= BYPASS && newMode < NUM_FILTER_MODES )
    mode = newMode;
  updateCoeffs();
}

void TwoPoleFilter::setFrequency(double newFrequency)
{
  frequency = clip(newFrequency, 2.0, 20000.0);
  updateCoeffs();
}

void TwoPoleFilter::setBandwidth(double newBandwidth)
{
  bandwidth = clip(newBandwidth, 0.25, 6.0);
  updateCoeffs();
}

void TwoPoleFilter::setRadius(double newRadius)
{
  radius = newRadius;
  updateCoeffs();
}

void TwoPoleFilter::setGain(double newGain)
{
  gain = newGain;
  updateCoeffs();
}

void TwoPoleFilter::setParameters(int newMode, double newFrequency, double newGain,
                                  double newBandwidth, bool /*updateCoefficients*/)
{
  if( newMode >= BYPASS && newMode < NUM_FILTER_MODES )
    mode = newMode;
  frequency = clip(newFrequency, 2.0, 20000.0);
  gain = newGain;
  bandwidth = clip(newBandwidth, 0.25, 6.0);
  updateCoeffs();
}

void TwoPoleFilter::setLowerBandedgeFrequency(double newLowerBandedgeFrequency)
{
  //setBandwidth(2.0*log2(frequency/newLowerBandedgeFrequency));
  setBandwidth(lowerBandedgeFrequencyToBandwdith(newLowerBandedgeFrequency, frequency));
}

void TwoPoleFilter::setUpperBandedgeFrequency(double newUpperBandedgeFrequency)
{
  //setBandwidth(2.0*log2(newUpperBandedgeFrequency/frequency));
  setBandwidth(upperBandedgeFrequencyToBandwdith(newUpperBandedgeFrequency, frequency));
}

//-----------------------------------------------------------------------------------------------------------------------------------------
// inquiry:

bool TwoPoleFilter::doesModeSupportBandwidth() const
{
  if(  mode == PEAK || mode == BANDREJECT || mode == LOW_SHELF || mode == HIGH_SHELF
    || mode == BANDPASS )
    return true;
  else
    return false;
}

bool TwoPoleFilter::doesModeSupportGain() const
{
  if(  mode == PEAK || mode == LOW_SHELF || mode == HIGH_SHELF
    || mode == LOWPASS12  || mode == HIGHPASS12)
    return true;
  else
    return false;
}

bool TwoPoleFilter::doesModeSupportRadius() const
{
  if( mode == REAL_POLE || mode == REAL_ZERO || mode == POLE_PAIR || mode == ZERO_PAIR )
    return true;
  else
    return false;
}

double TwoPoleFilter::lowerBandedgeFrequencyToBandwdith(double lowerBandEdgefrequency, double centerFrequency)
{
  return 2.0*log2(centerFrequency/lowerBandEdgefrequency);
}

double TwoPoleFilter::upperBandedgeFrequencyToBandwdith(double upperBandEdgefrequency, double centerFrequency)
{
  return  2.0*log2(upperBandEdgefrequency/centerFrequency);
}

//-----------------------------------------------------------------------------------------------------------------------------------------
// others:

void TwoPoleFilter::updateCoeffs()
{
  //double w0 = 2.0*PI*frequency/sampleRate;
  //double q2  = 1.0 / (2.0*sinh( 0.5*log(2.0) * bandwidth * w0/sin(w0) ) );
  //q = sqrt(0.5);

  double q2  = 1.0 / (2.0*sinh( 0.5*log(2.0) * bandwidth ));  // why q2 and not just q?

  switch(mode)
  {
  case BYPASS:
    BiquadDesigner::makeBypassBiquad(b0, b1, b2, a1, a2);
    break;
  case PEAK:
    BiquadDesigner::calculatePrescribedNyquistGainEqCoeffs(b0, b1, b2, a1, a2, 1.0/sampleRate, frequency, bandwidth, dB2amp(gain), 1.0);
    break;
  case LOW_SHELF:
    BiquadDesigner::calculateCookbookLowShelvCoeffs(b0, b1, b2, a1, a2, 1.0/sampleRate, frequency, q2, dB2amp(0.5*gain));
    // factor  0.5 because we define the characteristic frequency as the frequency where the dB-gain is one half of the dB-gain at DC

    //BiquadDesigner::calculateFirstOrderLowShelvCoeffsPrescribedNyQuist(b0, b1, b2, a1, a2, sampleRate, frequency, dB2amp(gain));

    break;
  case HIGH_SHELF:
    BiquadDesigner::calculateCookbookHighShelvCoeffs(b0, b1, b2, a1, a2, 1.0/sampleRate, frequency, q2, dB2amp(0.5*gain));
    //BiquadDesigner::calculateFirstOrderHighShelvCoeffsPrescribedNyQuist(b0, b1, b2, a1, a2, sampleRate, frequency, dB2amp(gain));
    break;
  case LOWPASS6:
    {
      //BiquadDesigner::calculateFirstOrderLowpassCoeffsBilinear(b0, b1, b2, a1, a2, 1.0/sampleRate, frequency);
      BiquadDesigner::calculateFirstOrderLowpassCoeffsPrescribedNyquist(b0, b1, b2, a1, a2, sampleRate, frequency);
      //BiquadDesigner::calculateFirstOrderLowpassCoeffs(b0, b1, b2, a1, a2, 1.0/sampleRate, frequency);
    }
    break;
  case LOWPASS12:
    BiquadDesigner::calculateCookbookLowpassCoeffs(b0, b1, b2, a1, a2, 1.0/sampleRate, frequency, dB2amp(gain));
    break;
  case HIGHPASS6:
    //BiquadDesigner::calculateFirstOrderHighpassCoeffsBilinear(b0, b1, b2, a1, a2, 1.0/sampleRate, frequency);
    BiquadDesigner::calculateFirstOrderHighpassCoeffsPrescribedNyquist(b0, b1, b2, a1, a2, sampleRate, frequency);
    break;
  case HIGHPASS12:
    BiquadDesigner::calculateCookbookHighpassCoeffs(b0, b1, b2, a1, a2, 1.0/sampleRate, frequency, dB2amp(gain));
    break;
  case BANDREJECT:
    BiquadDesigner::calculateCookbookBandrejectCoeffsViaBandwidth(b0, b1, b2, a1, a2, 1.0/sampleRate, frequency, bandwidth);
    break;
  case BANDPASS:
    BiquadDesigner::calculateCookbookBandpassConstSkirtCoeffsViaBandwidth(b0, b1, b2, a1, a2, 1.0/sampleRate, frequency, bandwidth);
    break;
      /*
  case LOWHIGH_PASS_CHAIN:
    {

    }
    break;
    */
  case REAL_POLE:
    BiquadDesigner::calculateOnePoleCoeffs(b0, b1, b2, a1, a2, radius);
    break;
  case REAL_ZERO:
    BiquadDesigner::calculateOneZeroCoeffs(b0, b1, b2, a1, a2, radius);
    break;
  case POLE_PAIR:
    {
      double tmpRadius = radius;
      if( fabs(radius) > 0.999999 )
        tmpRadius = sign(radius) * 0.999999;
      BiquadDesigner::calculatePolePairCoeffs(b0, b1, b2, a1, a2, tmpRadius, 2*PI*frequency/sampleRate);
    }
    break;
  case ZERO_PAIR:
    BiquadDesigner::calculateZeroPairCoeffs(b0, b1, b2, a1, a2, radius, 2*PI*frequency/sampleRate);
    break;
  default:
    BiquadDesigner::calculatePrescribedNyquistGainEqCoeffs(b0, b1, b2, a1, a2, 1.0/sampleRate, frequency, bandwidth, dB2amp(gain), 1.0);
  }
}



