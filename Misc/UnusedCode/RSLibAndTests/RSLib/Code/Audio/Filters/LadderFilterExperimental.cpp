using namespace RSLib;


rsLadderFilter::rsLadderFilter()
{
  sampleRate = 44100.0;
  cutoff     = 1000.0;
  resonance  = 0.0;
  setMode(LP_24);
  calcCoeffs();
  reset();
}

void rsLadderFilter::setCutoff(double newCutoff)
{
  cutoff = newCutoff;
  calcCoeffs();
}

void rsLadderFilter::setSampleRate(double newSampleRate)
{
  sampleRate = newSampleRate;
  calcCoeffs();
}

void rsLadderFilter::setResonance(double newResonance)
{
  resonance = newResonance;
  calcCoeffs();
}

void rsLadderFilter::setMode(int newMode)
{
  if( newMode >= 0 && newMode < NUM_MODES )
  {
    mode = newMode;
    switch(mode)
    {
    case FLAT:     c[0] = 1; c[1] =  0; c[2] =  0; c[3] =  0; c[4] =  0;  break;
    case LP_6:     c[0] = 0; c[1] =  1; c[2] =  0; c[3] =  0; c[4] =  0;  break;
    case LP_12:    c[0] = 0; c[1] =  0; c[2] =  1; c[3] =  0; c[4] =  0;  break;
    case LP_18:    c[0] = 0; c[1] =  0; c[2] =  0; c[3] =  1; c[4] =  0;  break;
    case LP_24:    c[0] = 0; c[1] =  0; c[2] =  0; c[3] =  0; c[4] =  1;  break;
    case HP_6:     c[0] = 1; c[1] = -1; c[2] =  0; c[3] =  0; c[4] =  0;  break;
    case HP_12:    c[0] = 1; c[1] = -2; c[2] =  1; c[3] =  0; c[4] =  0;  break;
    case HP_18:    c[0] = 1; c[1] = -3; c[2] =  3; c[3] = -1; c[4] =  0;  break;
    case HP_24:    c[0] = 1; c[1] = -4; c[2] =  6; c[3] = -4; c[4] =  1;  break;  
    case BP_6_6:   c[0] = 0; c[1] =  1; c[2] = -1; c[3] =  0; c[4] =  0;  break;
    case BP_6_12:  c[0] = 0; c[1] =  0; c[2] =  1; c[3] = -1; c[4] =  0;  break;
    case BP_6_18:  c[0] = 0; c[1] =  0; c[2] =  0; c[3] =  1; c[4] = -1;  break;
    case BP_12_6:  c[0] = 0; c[1] =  1; c[2] = -2; c[3] =  1; c[4] =  0;  break;
    case BP_12_12: c[0] = 0; c[1] =  0; c[2] =  1; c[3] = -2; c[4] =  1;  break;
    case BP_18_6:  c[0] = 0; c[1] =  1; c[2] = -3; c[3] =  3; c[4] = -1;  break;  
    default:       c[0] = 1; c[1] =  0; c[2] =  0; c[3] =  0; c[4] =  0;  // flat
    }
  }
}

void rsLadderFilter::getState(double *state)
{
  rsCopyBuffer(y, state, 5);
}

void rsLadderFilter::reset()
{
  rsFillWithZeros(y, 5);
}

double rsLadderFilter::computeCompensationGain(double a, double b, double k)
{
  double b4 = b*b*b*b;     // b^4
  return ((((a+4) * a+6) * a+4) * a + k*b4 + 1) / b4;
}

double rsLadderFilter::computeFeedbackFactor(double fb, double cosWc, double a, double b)
{
  double g2 = b*b / (1 + a*a + 2*a*cosWc);
  return fb / (g2*g2);
}

double rsLadderFilter::resonanceDecayToFeedbackGain(double decay, double cutoff)
{
  if(decay > 0.0)
    return exp(-1/(decay*cutoff));
  else
    return 0.0;
  // The time tr for a sinusoid at the cutoff frequency fc to complete a single roundtrip around 
  // the filter loop is given by tr = 1/fc. The amplitude of the sinusoid as function of t is given
  // by a(t) = r^(t/tr) = r^(t*fc). The decaytime is defined as the time instant where a(t) = 1/e, 
  // so we need to solve 1/e = r^(t*fc) leading to r = (1/e)^(1/(t*fc)).
  // can be expressed as exp(-1/(decay*cutoff)) - avoid expensive pow
}

void rsLadderFilter::computeCoeffs(double wc, double fb, double *a, double *b, double *k, 
  double *g)
{
  computeCoeffs(wc, fb, a, b, k);
  *g = computeCompensationGain(*a, *b, *k);
}

void rsLadderFilter::computeCoeffs(double wc, double fb, double *a, double *b, double *k)
{
  double s, c, t;                          // sin(wc), cos(wc), tan((wc-PI)/4)
  rsSinCos(wc, &s, &c);             
  t  = tan(0.25*(wc-PI));
  *a = t / (s-c*t);
  *b = 1.0 + *a;
  *k  = computeFeedbackFactor(fb, c, *a, *b);
}

void rsLadderFilter::calcCoeffs()
{
  double wc = 2 * PI * cutoff / sampleRate;
  computeCoeffs(wc, resonance, &a, &b, &k, &g);
}

//-------------------------------------------------------------------------------------------------
// This is the zero delay feedback implementation of the ladder filter. Let z1,z2,.. denote the old 
// values of y1,y2..., then the equations are:
// y0 =    x - k*y4
// y1 = b*y0 - a*z1 =          b*(x-k*y4) - a*z1
// y2 = b*y1 - a*z2 =       b*(b*(x-k*y4) - a*z1) - a*z2
// y3 = b*y2 - a*z3 =    b*(b*(b*(x-k*y4) - a*z1) - a*z2) - a*z3
// y4 = b*y3 - a*z4 = b*(b*(b*(b*(x-k*y4) - a*z1) - a*z2) - a*z3) - a*z4
//
// solve(y4 = b*(b*(b*(b*(x-k*y4) - a*z1) - a*z2) - a*z3) - a*z4, y4)
// gives:
// y4 = (-a b^3 z1-a b^2 z2-a b z3-a z4+b^4 x) / (b^4 k+1)
// we can write this as:
// y4 = p0*x + p1*z1 + p2*z2 + p3*z3 + p4*z4
// where:
// q = 1/(b^4*k + 1), p0 = q*b^4, p1 = -q*a*b^3, p2 = -q*a*b^2, p3 = -q*a*b, p4 = -q*a
//
// These p-coefficients are used for predicting the value of the 4th filterstage. Compared to the
// UDF (unit-delay-feedback) case, the computation of the a,b coefficients is a bit simpler because
// we dont have the additional delay into account (which simplifies the formula), but we have to
// additionaly compute the prediction coefficients and the per-sample calculations are more 
// expensive due to the evaluation of the prediction equation.

rsLadderFilterZDF::rsLadderFilterZDF()
{
  calcCoeffs();
}

void rsLadderFilterZDF::calcCoeffs()
{
  double wc, s, c, b4, q;

  wc = 2 * PI * cutoff / sampleRate;         
  rsSinCos(wc, &s, &c);               // s = sin(wc), c = cos(wc)
  a = -1 / (s+c);
  b =  1 + a;
  k = computeFeedbackFactor(resonance, c, a, b);
  g = computeCompensationGain(a, b, k);

  // prediction coefficients:
  b4   =  b*b*b*b;         // b^4
  q    =  1 / (b4*k + 1);  // scaler
  p[4] = -q*a;
  p[3] =  b*p[4];
  p[2] =  b*p[3];
  p[1] =  b*p[2];
  p[0] =  q*b4;

  // for optimization, we may copy the code from computeFeedbackFactor, computeCompensationGain and
  // then avoid some redundant computations (b^4, for example) ... but i'm not sure if it's worth 
  // it.
}

//-------------------------------------------------------------------------------------------------

rsLadderFilterFeedbackSaturated::rsLadderFilterFeedbackSaturated()
{ 
  drive   = 1.0;
  loLimit = -16777216;
  hiLimit = +16777216;
   // +-2^24 - large values - i think, it's good when they have an exact binary representation but 
   // that may not be too critical

  sigmoid.setValueAt1(1.0);  // hardclipping mode by default
  mode    = 0;
  computeScaleAndShift();
}

void rsLadderFilterFeedbackSaturated::setLowerFeedbackLimit(double newLimit)
{
  loLimit = newLimit;
  computeScaleAndShift();
}

void rsLadderFilterFeedbackSaturated::setUpperFeedbackLimit(double newLimit)
{
  hiLimit = newLimit;
  computeScaleAndShift();
}

void rsLadderFilterFeedbackSaturated::setFeedbackDrive(double newDrive)
{
  drive = newDrive;
}

void rsLadderFilterFeedbackSaturated::setSaturationGainAt1(double newGain)
{
  sigmoid.setValueAt1(newGain);
}

void rsLadderFilterFeedbackSaturated::setSaturationMode(int newMode)
{
  mode = newMode;
}

void rsLadderFilterFeedbackSaturated::computeScaleAndShift()
{
  rsRangeConversionCoefficients(loLimit, hiLimit,    -1.0,    +1.0, &scaleX, &shiftX);
  rsRangeConversionCoefficients(-1.0,       +1.0, loLimit, hiLimit, &scaleY, &shiftY);
}

double rsLadderFilterFeedbackSaturated::getCompensationGain()
{
  //double wc = 2*PI*cutoff/sampleRate;
  //double b0_1, a1_1, k_1;
  //computeCoeffs(wc, fb, &a1_1, &b0_1, &k_1, &g_1);

  double gain;
  gain = computeCompensationGain(a, b, k);  // maybe drive*k?
  //computeCompensationGain(a1, b0, drive*k, &gain);
  return gain;

  // this is under construction
  // maybe we should assume the following:
  // fb: the net feedback loop gain coming from the decay setting - a value between 0...1
  // drive: the additional multiplier coming from the driven saturation function
  // if fb*drive <= 1: compute DC gain using regular formula (using a k computed from
  //                   fb*drive (instead of just fb))
  // else: compute DC-gain assuming the net-feedback loop gain to be unity - using a k computed
  //       using a normalized feedback of 1)

}

//-------------------------------------------------------------------------------------------------

rsLadderResoShaped::rsLadderResoShaped()
{
  sampleRate  = 44100;
  inGain      = 1.0;
  leak        = 0.0;
  cutoff      = 1000;
  decay       = 0;
  decayByFreq = 0;
  scaledDecay = 0;
  attack      = 0;
  resGain     = 1.0;
  phase       = 0.0;

  allpass.setMode(rsOnePoleFilter::ALLPASS);
}

void rsLadderResoShaped::setSampleRate(double newSampleRate)
{
  sampleRate = newSampleRate;
  resonant.setSampleRate(sampleRate);
  nonResonant.setSampleRate(sampleRate);
  allpass.setSampleRate(sampleRate);
  setupAttackSmoother();
}

void rsLadderResoShaped::setInputGain(double newGain)
{
  inGain = newGain;
}

void rsLadderResoShaped::setInputLeak(double newLeak)
{
  leak = newLeak;
}

void rsLadderResoShaped::setCutoff(double newCutoff)
{
  cutoff = newCutoff;
  resonant.setCutoff(cutoff);
  nonResonant.setCutoff(cutoff);
  allpass.setCutoff(cutoff);
  setupScaledDecay();
}

void rsLadderResoShaped::setResonanceDecay(double newDecay)
{
  decay = newDecay;
  setupScaledDecay();
}

void rsLadderResoShaped::setDecayByFrequency(double newDecayByFreq)
{
  decayByFreq = newDecayByFreq;
  setupScaledDecay();
}

void rsLadderResoShaped::setResonanceAttack(double newAttack)
{
  attack = newAttack;
  setupAttackSmoother();
}

void rsLadderResoShaped::setResonanceGain(double newGain)
{
  resGain = newGain;
}

void rsLadderResoShaped::setResonancePhase(double newPhase)
{
  phase = newPhase;
}

void rsLadderResoShaped::setFeedbackDrive(double newDrive)
{
  resonant.setFeedbackDrive(newDrive);
    // todo: we may have to recalculate the DC gain compensation factor here
}

void rsLadderResoShaped::setFeedbackLowerLimit(double newLimit)
{
  resonant.setLowerFeedbackLimit(newLimit);
}

void rsLadderResoShaped::setFeedbackUpperLimit(double newLimit)
{
  resonant.setUpperFeedbackLimit(newLimit);
}

void rsLadderResoShaped::setFeedbackSaturationGainAt1(double newGain)
{
  resonant.setSaturationGainAt1(newGain);
}

void rsLadderResoShaped::setFeedbackSaturationPlace(int newPlace)
{
  resonant.setSaturationMode(newPlace);
}

void rsLadderResoShaped::getSignalParts(double in, double *yf, double *yr)
{
  double r; // temporary variable for resonance signal - maybe, get rid of it

  *yf = nonResonant.getSample(in);
  r   = resonant.getSample(in) - *yf;

  // todo: use getSampleNoGain function here, apply gain here...

  // apply attack smoothing filter:
  r = attackSmoother.getSample(r);

  // apply phase shift:
  double yq  = allpass.getSample(r);  // quadrature signal
  double gd, gq;                      // direct gain and quadrature gain
  rsSinCos(phase, &gq, &gd);          // make gd, gq members - move calculation to setResonancePhase
  r = gd*r + gq*yq;
  // apply phase-modulation by input, maybe the modulator can be pre-shaped as well...

  *yr = r;
}

double rsLadderResoShaped::getSample(double in)
{
  in *= inGain;
  double yf, yr;
  getSignalParts(in, &yf, &yr);
  return leak*in + yf + resGain*yr;
}

void rsLadderResoShaped::reset()
{
  resonant.reset();
  nonResonant.reset();
  attackSmoother.reset();
}

void rsLadderResoShaped::setupScaledDecay()
{
  scaledDecay = decay * pow(2.0, -decayByFreq*rsLog2(0.001*cutoff)); // neutral at cutoff = 1kHz
    // i think, this may be optimized so as to use only exp/log instead of pow/log2 (the latter
    // pair of which is more expensive to compute) by using some rules for exponentials and 
    // logarithms

  double r = rsLadderFilter::resonanceDecayToFeedbackGain(scaledDecay, cutoff);
  resonant.setResonance(r);
  
  // old:
  //resonant.setResonanceDecay(scaledDecay);

  setupAttackSmoother();
}

void rsLadderResoShaped::setupAttackSmoother()
{
  double w = 2*PI*cutoff/sampleRate;
  double d = attack*scaledDecay*sampleRate;
  attackSmoother.setFrequencyAndDecay(w, d);

  // enforce a unit-gain at the resonance frequency:
  attackSmoother.setOutputGain(1.0/attackSmoother.getMagnitudeAt(w)); 
}

//-------------------------------------------------------------------------------------------------

rsLadderResoShaped2::rsLadderResoShaped2()
{
  drive     = 1.0;
  addConst  = 0.0;
  addIn     = 0.0;
  addFlt    = 0.0; 
  driveComp = 0.0;
  unDrive   = 1.0;
  gate      = 1.0;
  gateMix   = 0.0;

  //saturator.setValueAt1(1.0); // hardclipper by default
  saturator.setMode(rsSidechainSaturator::BYPASS);  // no waveshaping by default
}

void rsLadderResoShaped2::setSaturationMode(int newMode)
{
  saturator.setMode(newMode);
}

void rsLadderResoShaped2::setSaturationGainAt1(double newGain)
{
  saturator.setValueAt1(newGain);
}

void rsLadderResoShaped2::setDrive(double newDrive)
{
  drive   = newDrive;
  unDrive = pow(1/drive, driveComp);
}

void rsLadderResoShaped2::setDriveCompensation(double newCompensation)
{
  driveComp = newCompensation;
  unDrive   = pow(1/drive, driveComp);
}

void rsLadderResoShaped2::setSaturationAddConstant(double newOffset)
{
  addConst = newOffset;
}

void rsLadderResoShaped2::setSaturationAddInput(double newAmount)
{
  addIn = newAmount;
}

void rsLadderResoShaped2::setSaturationAddFiltered(double newAmount)
{
  addFlt = newAmount;
}

void rsLadderResoShaped2::setGateSensitivity(double newSensitivity)
{
  gate = newSensitivity;
}

void rsLadderResoShaped2::setGateMix(double newMix)
{
  gateMix = newMix;
}

void rsLadderResoShaped2::getSignalParts(double in, double *yf, double *yr)
{
  rsLadderResoShaped::getSignalParts(in, yf, yr);

  // apply saturation:
  double yd = *yr * drive;                              // driven resonance
  double s  = yd + addIn*in + addFlt*(*yf) + addConst;  // saturator sidechain signal
  *yr  = saturator.getSample(yd, s);                    // actual saturation  
  *yr *= unDrive;                                       // gain compensation

                                                        // apply gate (preliminary):
  if(gate > 0.0)
  {
    double gateLo = resonant.getLowerFeedbackLimit() / gate;
    double gateHi = resonant.getUpperFeedbackLimit() / gate;
    double tmp    = gateMix * (*yf) + (1-gateMix) * in;
    if(tmp < gateLo || tmp > gateHi)
      *yr = 0;  
    // This is a very crude hard-gate - later, a soft gate should be used and/or a filter applied
    // to avoid flutter with bandlimited waveforms. The code can also be optimized by precomputing
    // gateLo/gateHi. I think, it would be best, to wrap the gate-code into a class.
  }
}

//-------------------------------------------------------------------------------------------------

rsResoReplacer::rsResoReplacer()
{
  setResonanceWaveform(0);        // sine
  resoCutoffMultiplier = 10000.0; // maps to 20kHz at cutoff = 20Hz -> high enough to be neutral
  resoCutoffModulation = 0.0;

  resoFilter.setMode(rsOnePoleFilter::LOWPASS);  
    // maybe let the user choose, maybe use more interesting filter - like a ladder or biquad
}

void rsResoReplacer::setResonanceWaveform(int newWaveform)
{
  waveform = newWaveform;
    // later, we need to set up the wavetable oscillator here...
}

void rsResoReplacer::setResoCutoffMultiplier(double newMultiplier)
{
  resoCutoffMultiplier = newMultiplier;
}

void rsResoReplacer::setResoCutoffModulation(double newModulation)
{
  resoCutoffModulation = newModulation;
}

void rsResoReplacer::setSampleRate(double newSampleRate)
{
  rsLadderResoShaped::setSampleRate(newSampleRate);
  resoFilter.setSampleRate(newSampleRate);
}

void rsResoReplacer::getSignalParts(double in, double *yf, double *yr)
{
  double mag, phs;
  getNonresonantOutputAndResonanceParameters(in, yf, &mag, &phs);
  *yr = getWaveForm(mag, phs);
}

void rsResoReplacer::getNonresonantOutputAndResonanceParameters(double in, double *yf,
  double *mag, double *phs)
{
  // temporary variables:
  double r, q;

  *yf = nonResonant.getSample(in);
  r   = resonant.getSampleLinear(in) - *yf;   // pure resonance
  r   = attackSmoother.getSample(r);          // attack smoothing filter

  // compute magnitude and phase:
  q     = allpass.getSample(r);                // quadrature signal
  *mag  = sqrt(r*r + q*q);                    // extract magnitude
  *phs  = atan2(r, q);                        // exctract phase

  // manipulate phase (later add phase-modulation):
  *phs += phase;  // add static offset
  *phs += 2*PI;   // rsSqrWave, etc. expects phase in 0...2*PI - get rid when using a wavetable
}

double rsResoReplacer::getWaveForm(double a, double p)
{
  double r;

  // resonance waveform reconstruction (prelimiary - later use a mip mapped wavetable, where the
  // mip map level is chosen according the the phase difference between previous and current 
  // sample):
  switch(waveform)
  {
  case 1:  r =  a * rsTriWave(p); break;
  case 2:  r =  a * rsSqrWave(p); break;
  case 3:  r =  a * rsSawWave(p); break;
  case 4:  r = -a * rsSawWave(p); break;
  default: r =  a * sin(p);
  }

  // set up the resonance filter's cutoff frequency and filter the resonance waveform:
  double ctf = cutoff*(resoCutoffMultiplier + resoCutoffModulation*a);
  ctf = rsLimitToRange(ctf, 20.0, 20000.0);
  resoFilter.setCutoff(ctf);
  r = resoFilter.getSample(r);
  // todo: optimization: skip setCutoff if resoCutoffModulation == 0 - but we need to take care
  // that it's called in setResoCutoffMultiplier, setResoCutoffModulation, setCutoff

  //r /= resoFilter.getMagnitudeAt(cutoff); // gain correction - maybe not such a good idea

  return r;
}

//-------------------------------------------------------------------------------------------------

rsResoReplacerPhaseBumped::rsResoReplacerPhaseBumped()
{
  // todo: set up highpass and lowpass filters for the bumping signal

  inputDifferencer.setMode(rsOnePoleFilter::HIGHPASS);
  bumpSmoother1.setMode(rsOnePoleFilter::LOWPASS);
  bumpSmoother2.setMode(rsOnePoleFilter::LOWPASS);


  inputDifferencer.setCutoff(20000);

  // preliminary:
  bumpCutoff1 = 20;
  bumpCutoff2 = 100;
  bumpFactor  = 5.0;         // maybe this should be related to the 2 cutoffs
  bumpSmoother1.setCutoff(bumpCutoff1);
  bumpSmoother2.setCutoff(bumpCutoff2);

  averageSamples = 2000;
  bumpAverager.setLengthInSamples(averageSamples);
  bumpAverager.reset();


  setChaosParameter(0.0);

  setAmplitudeLimit(10.0);
  setInputRange(1.0);

  setSelfExitation(0.0);

  // move into reset function:
  bump     = 0.0;
  oldPhase = 0.0;
}

void rsResoReplacerPhaseBumped::setSampleRate(double newSampleRate)
{
  rsResoReplacer::setSampleRate(newSampleRate);
  inputDifferencer.setSampleRate(newSampleRate);
  bumpSmoother1.setSampleRate(newSampleRate);
  bumpSmoother2.setSampleRate(newSampleRate);
}

void rsResoReplacerPhaseBumped::setBumpFactor(double newFactor)
{
  bumpFactor = newFactor;
}

void rsResoReplacerPhaseBumped::setBumpCutoffs(double cutoff1, double cutoff2)
{
  bumpCutoff1 = cutoff1;
  bumpCutoff2 = cutoff2;
  bumpSmoother1.setCutoff(bumpCutoff1);
  bumpSmoother2.setCutoff(bumpCutoff2);
}

void rsResoReplacerPhaseBumped::setChaosParameter(double newParameter)
{
  chaosParam = newParameter;
}

void rsResoReplacerPhaseBumped::setAmplitudeLimit(double newLimit)
{
  rsAssert(newLimit > 0.0);
  ampLimit = newLimit;
  ampLimitInv = 1.0 / ampLimit;
}

void rsResoReplacerPhaseBumped::setInputRange(double newRange)
{
  inRange = newRange;
}

void rsResoReplacerPhaseBumped::setSelfExitation(double newExcitation)
{
  selfExcite = newExcitation;
}

void rsResoReplacerPhaseBumped::getSignalParts(double in, double *yf, double *yr)
{
  double ex = 0.0;      // exitation signal
  ex = getExcitation(); // comment to deactivate self-excitation code

  // copied from baseclass:
  double mag, phs;
  getNonresonantOutputAndResonanceParameters(in+ex, yf, &mag, &phs);
  mag  = limitMagnitude(mag);     // new - maybe optimize by inlining (2-level)
  *yr  = getWaveForm(mag, phs);  
  *yr *= inputScale(*yf);         // lowpass based limiting

  updateBump(in, *yf, *yr);            // also new
}

void rsResoReplacerPhaseBumped::updateBump(double in, double yf, double yr)
{
  double tmp;
  tmp  = inputDifferencer.getSample(in);
  tmp *= in + chaosParam * yr;   // should produce the chaos
  tmp *= bumpFactor;

  // apply moving average filter (turn spikes into plateaus):
  //tmp = averageSamples * bumpAverager.getSample(tmp);
  //tmp = bumpAverager.getSample(tmp);

  //tmp  = rsWrapToInterval(tmp, -PI, PI); // seems not good
  tmp  = bumpSmoother1.getSample(tmp);
  tmp  = bumpSmoother2.getSample(tmp);
  tmp *=  10000.0 / rsMin(bumpCutoff1, bumpCutoff2);// preliminary: use peak excursion formula, or better: attck/decay filter
   
  bump = tmp;

  oldPhase += bump;
}

double rsResoReplacerPhaseBumped::limitMagnitude(double mag)
{
  //return ampLimit * rsPositiveSigmoids::softClip(ampLimitInv*mag);

  return ampLimit * rsPositiveSigmoids::quartic(ampLimitInv*mag);
    // maybe try a different function
}

double rsResoReplacerPhaseBumped::inputScale(double in)
{
  // we use a kind of "Butterworth"-squared function for the multiplier

  double tmp; 

  tmp  = in / inRange;   // maybe include an offset later, use multiplication instead of division
  tmp *= tmp;            // ^2  -> 1st  order butterworth
  tmp *= tmp;            // ^4  -> 2nd  order   
  tmp *= tmp;            // ^8  -> 4th  order
  tmp *= tmp;            // ^16 -> 8th  order
  tmp *= tmp;            // ^32 -> 16th order
  tmp  = 1 / (1 + tmp);  

  return tmp;

  // todo: maybe use a piecewise polynomial function here that actually goes down to zero - we can
  // have a "transition width" parameter
}

double rsResoReplacerPhaseBumped::getExcitation()
{
  // experimental - supposed to work only when oldPhase is NOT updated in getSignalParts - 
  // excitation independent from input / feedback
  oldPhase += 2*PI*cutoff/sampleRate;
  return selfExcite * sin(oldPhase);
    // OK - this works and produces the beating / amplitude modulation effect - but no growl
    // maybe, we should use this and manipulate "oldPhase" in "updateBump"

  // maybe have also a self-excitation phase-offset (user parameter) and maybe we should add
  // "omega" to the oldPhase because the oldPhase is, well, old (lagging one omega behind)

}

//-------------------------------------------------------------------------------------------------

rsLadderMystran::rsLadderMystran()
{ 
  setMode(LP_24);
  reset(); 
}

void rsLadderMystran::reset()
{
  rsFillWithZeros(s, 4);
  zi = 0.0;
}
