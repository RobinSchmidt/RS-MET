#include "rosic_CombResonator.h"
using namespace rosic;

//-------------------------------------------------------------------------------------------------
// construction/destruction:

CombResonator::CombResonator(int bufferLengthToAllocate) : CombFilter(bufferLengthToAllocate)
{
  decayTime      = 3.0;   
  highDecayScale = 0.25; 
  lowDecayScale  = 2.0;
  highCrossOver  = 4000.0;
  lowCrossOver   = 250;
  oddOnlyMode    = false;

  setupDelayInSamples();
  //setupFilters();
}

CombResonator::~CombResonator()
{

}

//-------------------------------------------------------------------------------------------------
// parameter settings:

void CombResonator::setSampleRate(double newSampleRate)
{
  CombFilter::setSampleRate(newSampleRate);
  feedbackFilter.setSampleRate(newSampleRate);
  correctionFilter.setSampleRate(newSampleRate);
  setupDelayInSamples();
}

void CombResonator::setFrequency(double newFrequency)
{
  CombFilter::setFrequency(newFrequency);
  setupDelayInSamples();
}

void CombResonator::setDecayTime(double newDecayTime)
{
  if( newDecayTime >= 0.0 )
  {
    decayTime = newDecayTime;
    setupFilters();
  }
}

void CombResonator::setHighDecayScale(double newScale)
{
  if( newScale >= 0.0 )
  {
    highDecayScale = newScale;
    setupFilters();
  }
}

void CombResonator::setLowDecayScale(double newScale)
{
  if( newScale >= 0.0 )
  {
    lowDecayScale = newScale;
    setupFilters();
  }
}

void CombResonator::setHighCrossoverFreq(double newFreq)
{
  if( newFreq >= 0.0 )
  {
    highCrossOver= newFreq;
    setupFilters();
  }
}

void CombResonator::setLowCrossoverFreq(double newFreq)
{
  if( newFreq >= 0.0 )
  {
    lowCrossOver= newFreq;
    setupFilters();
  }
}

void CombResonator::setOddOnlyMode(bool shouldCreateOnlyOddHarmonics)
{
  oddOnlyMode = shouldCreateOnlyOddHarmonics;
  setupDelayInSamples();
  //setupFilters();
}

//-------------------------------------------------------------------------------------------------
// others:

void CombResonator::clearBuffer()
{
  CombFilter::clearBuffer();
  feedbackFilter.reset();
}

void CombResonator::setupDelayInSamples()
{
  double delayInSeconds;
  if( !oddOnlyMode )
    delayInSeconds = 1.0 / frequency;
  else
    delayInSeconds = 0.5 / frequency;

  delayInSamples = sampleRate*delayInSeconds;
  delayInSamples = clip(delayInSamples, (double)(interpolatorMargin-1), 
    (double) (length-1-interpolatorMargin));

  // calculate the integer and fractional parts of the delay:
  double tmp   = floor(delayInSamples);
  int    dInt  = (int) tmp;
  double dFrac = delayInSamples - tmp;
  frac         = 1.0 - dFrac; // because we look backwards

  // adjust tapOut-pointer:
  tapOut = tapIn - dInt - 1;
  if( frac >= 1.0 )
  {
    frac    = 0.0;
    tapOut += 1;
  }
  tapOut = wrapAround(tapOut, length);


  setupFilters();
}

void CombResonator::setupFilters()
{
  // todo: move these computations into class DampingFilter

  double T_l, T_lb, T_m, T_hb, T_h;    // desired decay-times at 5 frequencies
  double g_l_abs, g_m_abs, g_h_abs;    // absolute gain-factors for the 3 bands
  double g_lb_abs, g_hb_abs;           // absolute gains at which the bandwidths are measured
                                       // (using the terminology from the Orfanidis-papaer)
  double g_l_rel, g_h_rel;             // gains for low and high band relative to the mid-band
  double g_lb_rel, g_hb_rel;           // relative gains at the bandwidth measurement frequencies

  T_l  = lowDecayScale*decayTime;      // desired decay-time at low frequencies
  T_m  = decayTime;                    // desired decay-time at mid frequencies
  T_h  = highDecayScale*decayTime;     // desired decay-time at high frequencies
  T_lb = sqrt(T_l*T_m);                // desired decay-time at the low crossover frequency
  T_hb = sqrt(T_h*T_m);                // desired decay-time at the high crossover frequency

  // calculate the required absolute gain-factors at 5 frequencies:
  g_l_abs  = pow(10.0, -3.0*delayInSamples/(T_l* sampleRate) );
  g_lb_abs = pow(10.0, -3.0*delayInSamples/(T_lb*sampleRate) );
  g_m_abs  = pow(10.0, -3.0*delayInSamples/(T_m* sampleRate) );
  g_hb_abs = pow(10.0, -3.0*delayInSamples/(T_hb*sampleRate) );
  g_h_abs  = pow(10.0, -3.0*delayInSamples/(T_h* sampleRate) );

  // the desired absolute gains will be approximately realized by a low- and high-shelving filter 
  // and a global gain, where the global gain is equal to the desired mid-frequency gain and the 
  // shelves are set up to account for the relative gain-deviation:
  g_l_rel  = g_l_abs  / g_m_abs;
  g_lb_rel = g_lb_abs / g_m_abs;
  g_h_rel  = g_h_abs  / g_m_abs;
  g_hb_rel = g_hb_abs / g_m_abs;

  // set up the damping-filter for delayline d:
  if( oddOnlyMode == true )
    feedbackFilter.setGlobalGainFactor(-g_m_abs);
  else
    feedbackFilter.setGlobalGainFactor(g_m_abs);

  feedbackFilter.setLowCrossoverFreq(lowCrossOver);
  feedbackFilter.setLowCrossoverGainFactor(g_lb_rel);
  feedbackFilter.setLowGainFactor(g_l_rel);
  feedbackFilter.setHighCrossoverFreq(highCrossOver);
  feedbackFilter.setHighCrossoverGainFactor(g_hb_rel);
  feedbackFilter.setHighGainFactor(g_h_rel);

  // set up the correction filters which decouple the overall frequency response from the frequency 
  // dependent decay times:
  correctionFilter.setGlobalGainFactor(1.0/sqrt(T_m));
  correctionFilter.setLowGainFactor(1.0/sqrt(T_l/T_m));
  correctionFilter.setLowCrossoverGainFactor(1.0/sqrt(T_lb/T_m));
  correctionFilter.setHighGainFactor(1.0/sqrt(T_h/T_m));
  correctionFilter.setHighCrossoverGainFactor(1.0/sqrt(T_hb/T_m));
  correctionFilter.setLowCrossoverFreq(lowCrossOver);
  correctionFilter.setHighCrossoverFreq(highCrossOver);
}