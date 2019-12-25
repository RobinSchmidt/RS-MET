template<class TSig, class TPar>
rsLadderFilter2<TSig, TPar>::rsLadderFilter2()
{
  sampleRate = 44100.0;
  cutoff     = 1000.0;
  resonance  = 0.0;
  setMode(LP_24);
  calcCoeffs();
  reset();
}

template<class TSig, class TPar>
void rsLadderFilter2<TSig, TPar>::setCutoff(TPar newCutoff)
{
  cutoff = newCutoff;
  calcCoeffs();
}

template<class TSig, class TPar>
void rsLadderFilter2<TSig, TPar>::setSampleRate(TPar newSampleRate)
{
  sampleRate = newSampleRate;
  calcCoeffs();
}

template<class TSig, class TPar>
void rsLadderFilter2<TSig, TPar>::setResonance(TPar newResonance)
{
  resonance = newResonance;
  calcCoeffs();
}

template<class TSig, class TPar>
void rsLadderFilter2<TSig, TPar>::setMode(int newMode)
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

template<class TSig, class TPar>
void rsLadderFilter2<TSig, TPar>::getState(TSig *state)
{
  rsArrayTools::copy(y, state, 5);
}

template<class TSig, class TPar>
void rsLadderFilter2<TSig, TPar>::reset()
{
  rsArrayTools::fillWithZeros(y, 5);
}

template<class TSig, class TPar>
TPar rsLadderFilter2<TSig, TPar>::computeCompensationGain(TPar a, TPar b, TPar k)
{
  double b4 = b*b*b*b;     // b^4
  return ((((a+4) * a+6) * a+4) * a + k*b4 + 1) / b4;
  // todo: look at the limit when wc -> 0. then a -> 1, b -> 0
}

template<class TSig, class TPar>
TPar rsLadderFilter2<TSig, TPar>::computeFeedbackFactor(TPar fb, TPar cosWc, TPar a, TPar b)
{
  TPar g2 = b*b / (1 + a*a + 2*a*cosWc);
  return fb / (g2*g2);
}

template<class TSig, class TPar>
TPar rsLadderFilter2<TSig, TPar>::resonanceDecayToFeedbackGain(TPar decay, TPar cutoff)
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

template<class TSig, class TPar>
void rsLadderFilter2<TSig, TPar>::computeCoeffs(TPar wc, TPar fb, TPar *a, TPar *b, TPar *k,
  TPar *g)
{
  computeCoeffs(wc, fb, a, b, k);
  *g = computeCompensationGain(*a, *b, *k);
}

template<class TSig, class TPar>
void rsLadderFilter2<TSig, TPar>::computeCoeffs(TPar wc, TPar fb, TPar *a, TPar *b, TPar *k)
{
  TPar s, c, t;                     // sin(wc), cos(wc), tan((wc-PI)/4)
  rsSinCos(wc, &s, &c);
  t  = tan(0.25*(wc-PI));
  *a = t / (s-c*t);
  *b = 1.0 + *a;
  *k = computeFeedbackFactor(fb, c, *a, *b);
  // If the cutoff frequency goes to zero wc -> 0, then s -> 0, c -> 1, t -> -1, b -> 0, a -> -1.
  // The coefficient computation for the lowpass stages approaches this limit correctly, but the
  // formula for the feedback factor k runs into numerical problems when wc -> 0. However, we know
  // from the analog filter that the correct limit for the feedback factor is k = 4*fb, since the
  // analog limit corresponds to infinite samplerate. Also, the formula for the compensation gain
  // becomes invalid. Maybe we should look at the magnitude response at zero frequency of the
  // analog filter to get a gain formula from the feedback factor for this. Then, we could use a
  // lower threshold for wc, below which these limiting case formulas are used. For the feedback
  // computation, a lower limit of wc = 1.e-6 seems appropriate (for gain compensation, i did not
  // yet figure it out)
}

template<class TSig, class TPar>
void rsLadderFilter2<TSig, TPar>::calcCoeffs()
{
  TPar wc = 2 * PI * cutoff / sampleRate;
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

template<class TSig, class TPar>
rsLadderFilterZDF<TSig, TPar>::rsLadderFilterZDF()
{
  calcCoeffs();
}

template<class TSig, class TPar>
void rsLadderFilterZDF<TSig, TPar>::calcCoeffs()
{
  TPar wc, s, c, b4, q;

  wc = 2 * PI * this->cutoff / this->sampleRate;
  rsSinCos(wc, &s, &c);               // s = sin(wc), c = cos(wc)
  this->a = -1 / (s+c);
  this->b =  1 + this->a;
  this->k = this->computeFeedbackFactor(this->resonance, c, this->a, this->b);
  this->g = this->computeCompensationGain(this->a, this->b, this->k);

  // prediction coefficients:
  b4   =  this->b * this->b * this->b * this->b;   // b^4
  q    =  1 / (b4 * this->k + 1);                  // scaler
  p[4] = -q * this->a;
  p[3] =  this->b * p[4];
  p[2] =  this->b * p[3];
  p[1] =  this->b * p[2];
  p[0] =  q * b4;

  // for optimization, we may copy the code from computeFeedbackFactor, computeCompensationGain and
  // then avoid some redundant computations (b^4, for example) ... but i'm not sure if it's worth
  // it.
}

//-------------------------------------------------------------------------------------------------

template<class TSig, class TPar>
rsLadderFilterFeedbackSaturated<TSig, TPar>::rsLadderFilterFeedbackSaturated()
{
  drive   = 1.0;
  loLimit = -16777216;
  hiLimit = +16777216;
   // +-2^24 - large values - i think, it's good when they have an exact binary representation but
   // that may not be too critical

  sigmoid.setValueAt1(1.0);  // hardclipping mode by default
  mode    = 0;
  computeScaleAndShift();


  fbLpCutoff = 1.0;
  fbLpf.setMode(rsOnePoleFilter<TSig, TPar>::LOWPASS_IIT);
}

template<class TSig, class TPar>
void rsLadderFilterFeedbackSaturated<TSig, TPar>::setLowerFeedbackLimit(TPar newLimit)
{
  loLimit = newLimit;
  computeScaleAndShift();
}

template<class TSig, class TPar>
void rsLadderFilterFeedbackSaturated<TSig, TPar>::setUpperFeedbackLimit(TPar newLimit)
{
  hiLimit = newLimit;
  computeScaleAndShift();
}

template<class TSig, class TPar>
void rsLadderFilterFeedbackSaturated<TSig, TPar>::setFeedbackDrive(TPar newDrive)
{
  drive = newDrive;
}

template<class TSig, class TPar>
void rsLadderFilterFeedbackSaturated<TSig, TPar>::setSaturationGainAt1(TPar newGain)
{
  sigmoid.setValueAt1(newGain);
}

template<class TSig, class TPar>
void rsLadderFilterFeedbackSaturated<TSig, TPar>::setSaturationMode(int newMode)
{
  mode = newMode;
}

template<class TSig, class TPar>
void rsLadderFilterFeedbackSaturated<TSig, TPar>::computeScaleAndShift()
{
  rsRangeConversionCoefficients(loLimit, hiLimit,    -1.0,    +1.0, &scaleX, &shiftX);
  rsRangeConversionCoefficients(-1.0,       +1.0, loLimit, hiLimit, &scaleY, &shiftY);
}

template<class TSig, class TPar>
TPar rsLadderFilterFeedbackSaturated<TSig, TPar>::getCompensationGain()
{
  //double wc = 2*PI*cutoff/sampleRate;
  //double b0_1, a1_1, k_1;
  //computeCoeffs(wc, fb, &a1_1, &b0_1, &k_1, &g_1);

  TPar gain;
  gain = this->computeCompensationGain(this->a, this->b, this->k);  // maybe drive*k?
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

template<class TSig, class TPar>
rsLadderResoShaped<TSig, TPar>::rsLadderResoShaped()
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

  allpass.setMode(rsOnePoleFilter<TSig, TPar>::ALLPASS_BLT);
}

template<class TSig, class TPar>
void rsLadderResoShaped<TSig, TPar>::setSampleRate(TPar newSampleRate)
{
  sampleRate = newSampleRate;
  resonant.setSampleRate(sampleRate);
  nonResonant.setSampleRate(sampleRate);
  allpass.setSampleRate(sampleRate);
  setupAttackSmoother();
}

template<class TSig, class TPar>
void rsLadderResoShaped<TSig, TPar>::setInputGain(TPar newGain)
{
  inGain = newGain;
}

template<class TSig, class TPar>
void rsLadderResoShaped<TSig, TPar>::setInputLeak(TPar newLeak)
{
  leak = newLeak;
}

template<class TSig, class TPar>
void rsLadderResoShaped<TSig, TPar>::setCutoff(TPar newCutoff)
{
  cutoff = newCutoff;
  resonant.setCutoff(cutoff);
  nonResonant.setCutoff(cutoff);
  allpass.setCutoff(cutoff);
  setupScaledDecay();
}

template<class TSig, class TPar>
void rsLadderResoShaped<TSig, TPar>::setResonanceDecay(TPar newDecay)
{
  decay = newDecay;
  setupScaledDecay();
}

template<class TSig, class TPar>
void rsLadderResoShaped<TSig, TPar>::setDecayByFrequency(TPar newDecayByFreq)
{
  decayByFreq = newDecayByFreq;
  setupScaledDecay();
}

template<class TSig, class TPar>
void rsLadderResoShaped<TSig, TPar>::setResonanceAttack(TPar newAttack)
{
  attack = newAttack;
  setupAttackSmoother();
}

template<class TSig, class TPar>
void rsLadderResoShaped<TSig, TPar>::setResonanceGain(TPar newGain)
{
  resGain = newGain;
}

template<class TSig, class TPar>
void rsLadderResoShaped<TSig, TPar>::setResonancePhase(TPar newPhase)
{
  phase = newPhase;
}

template<class TSig, class TPar>
void rsLadderResoShaped<TSig, TPar>::setFeedbackDrive(TPar newDrive)
{
  resonant.setFeedbackDrive(newDrive);
    // todo: we may have to recalculate the DC gain compensation factor here
}

template<class TSig, class TPar>
void rsLadderResoShaped<TSig, TPar>::setFeedbackLowerLimit(TPar newLimit)
{
  resonant.setLowerFeedbackLimit(newLimit);
}

template<class TSig, class TPar>
void rsLadderResoShaped<TSig, TPar>::setFeedbackUpperLimit(TPar newLimit)
{
  resonant.setUpperFeedbackLimit(newLimit);
}

template<class TSig, class TPar>
void rsLadderResoShaped<TSig, TPar>::setFeedbackSaturationGainAt1(TPar newGain)
{
  resonant.setSaturationGainAt1(newGain);
}

template<class TSig, class TPar>
void rsLadderResoShaped<TSig, TPar>::setFeedbackSaturationPlace(int newPlace)
{
  resonant.setSaturationMode(newPlace);
}

template<class TSig, class TPar>
void rsLadderResoShaped<TSig, TPar>::getSignalParts(TSig in, TSig *yf, TSig *yr)
{
  double r; // temporary variable for resonance signal - maybe, get rid of it

  *yf = nonResonant.getSample(in);
  r   = resonant.getSample(in) - *yf;

  // todo: use getSampleNoGain function here, apply gain here...

  // apply attack smoothing filter:
  r = attackSmoother.getSample(r);

  // apply phase shift:
  TSig yq  = allpass.getSample(r);  // quadrature signal
  TSig gd, gq;                      // direct gain and quadrature gain
  rsSinCos(phase, &gq, &gd);          // make gd, gq members - move calculation to setResonancePhase
  r = gd*r + gq*yq;
  // apply phase-modulation by input, maybe the modulator can be pre-shaped as well...

  *yr = r;
}

template<class TSig, class TPar>
TSig rsLadderResoShaped<TSig, TPar>::getSample(TSig in)
{
  in *= inGain;
  TSig yf, yr;
  getSignalParts(in, &yf, &yr);
  return leak*in + yf + resGain*yr;
}

template<class TSig, class TPar>
void rsLadderResoShaped<TSig, TPar>::reset()
{
  resonant.reset();
  nonResonant.reset();
  attackSmoother.reset();
}

template<class TSig, class TPar>
void rsLadderResoShaped<TSig, TPar>::setupScaledDecay()
{
  scaledDecay = decay * pow(2.0, -decayByFreq*rsLog2(0.001*cutoff)); // neutral at cutoff = 1kHz
    // i think, this may be optimized so as to use only exp/log instead of pow/log2 (the latter
    // pair of which is more expensive to compute) by using some rules for exponentials and
    // logarithms

  TPar r = rsLadderFilter<TSig, TPar>::resonanceDecayToFeedbackGain(scaledDecay, cutoff);
  resonant.setResonance(r);

  // old:
  //resonant.setResonanceDecay(scaledDecay);

  setupAttackSmoother();
}

template<class TSig, class TPar>
void rsLadderResoShaped<TSig, TPar>::setupAttackSmoother()
{
  TPar w = 2*PI*cutoff/sampleRate;
  TPar d = attack*scaledDecay*sampleRate;
  attackSmoother.setFrequencyAndDecay(w, d);

  // enforce a unit-gain at the resonance frequency:
  attackSmoother.setOutputGain(1.0/attackSmoother.getMagnitudeAt(w));
}

//-------------------------------------------------------------------------------------------------

template<class TSig, class TPar>
rsLadderResoShaped2<TSig, TPar>::rsLadderResoShaped2()
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
  saturator.setMode(rsSidechainSaturator<TSig, TPar>::BYPASS);  // no waveshaping by default
}

template<class TSig, class TPar>
void rsLadderResoShaped2<TSig, TPar>::setSaturationMode(int newMode)
{
  saturator.setMode(newMode);
}

template<class TSig, class TPar>
void rsLadderResoShaped2<TSig, TPar>::setSaturationGainAt1(TPar newGain)
{
  saturator.setValueAt1(newGain);
}

template<class TSig, class TPar>
void rsLadderResoShaped2<TSig, TPar>::setDrive(TPar newDrive)
{
  drive   = newDrive;
  unDrive = pow(1/drive, driveComp);
}

template<class TSig, class TPar>
void rsLadderResoShaped2<TSig, TPar>::setDriveCompensation(TPar newCompensation)
{
  driveComp = newCompensation;
  unDrive   = pow(1/drive, driveComp);
}

template<class TSig, class TPar>
void rsLadderResoShaped2<TSig, TPar>::setSaturationAddConstant(TPar newOffset)
{
  addConst = newOffset;
}

template<class TSig, class TPar>
void rsLadderResoShaped2<TSig, TPar>::setSaturationAddInput(TPar newAmount)
{
  addIn = newAmount;
}

template<class TSig, class TPar>
void rsLadderResoShaped2<TSig, TPar>::setSaturationAddFiltered(TPar newAmount)
{
  addFlt = newAmount;
}

template<class TSig, class TPar>
void rsLadderResoShaped2<TSig, TPar>::setGateSensitivity(TPar newSensitivity)
{
  gate = newSensitivity;
}

template<class TSig, class TPar>
void rsLadderResoShaped2<TSig, TPar>::setGateMix(TPar newMix)
{
  gateMix = newMix;
}

template<class TSig, class TPar>
void rsLadderResoShaped2<TSig, TPar>::getSignalParts(TSig in, TSig *yf, TSig *yr)
{
  rsLadderResoShaped<TSig, TPar>::getSignalParts(in, yf, yr);

  // apply saturation:
  TSig yd = *yr * drive;                              // driven resonance
  TSig s  = yd + addIn*in + addFlt*(*yf) + addConst;  // saturator sidechain signal
  *yr  = saturator.getSample(yd, s);                    // actual saturation
  *yr *= unDrive;                                       // gain compensation


  // to be removed:
  //// apply gate (preliminary):
  //if(gate > 0.0)
  //{
  //  double gateLo = resonant.getLowerFeedbackLimit() / gate;
  //  double gateHi = resonant.getUpperFeedbackLimit() / gate;
  //  double tmp    = gateMix * (*yf) + (1-gateMix) * in;
  //  if(tmp < gateLo || tmp > gateHi)
  //    *yr = 0;
  //  // This is a very crude hard-gate - later, a soft gate should be used and/or a filter applied
  //  // to avoid flutter with bandlimited waveforms. The code can also be optimized by precomputing
  //  // gateLo/gateHi. I think, it would be best, to wrap the gate-code into a class.
  //}
}

//-------------------------------------------------------------------------------------------------

template<class TSig, class TPar>
rsResoReplacer<TSig, TPar>::rsResoReplacer()
{
  setResonanceWaveform(0);        // sine
  resoCutoffMultiplier = 10000.0; // maps to 20kHz at cutoff = 20Hz -> high enough to be neutral
  resoCutoffModulation = 0.0;

  resoFilter.setMode(rsOnePoleFilter<TSig, TPar>::LOWPASS_IIT);
    // maybe let the user choose, maybe use more interesting filter - like a ladder or biquad
}

template<class TSig, class TPar>
void rsResoReplacer<TSig, TPar>::setResonanceWaveform(int newWaveform)
{
  waveform = newWaveform;
    // later, we need to set up the wavetable oscillator here...
}

template<class TSig, class TPar>
void rsResoReplacer<TSig, TPar>::setResoCutoffMultiplier(TPar newMultiplier)
{
  resoCutoffMultiplier = newMultiplier;
}

template<class TSig, class TPar>
void rsResoReplacer<TSig, TPar>::setResoCutoffModulation(TPar newModulation)
{
  resoCutoffModulation = newModulation;
}

template<class TSig, class TPar>
void rsResoReplacer<TSig, TPar>::setSampleRate(TPar newSampleRate)
{
  rsLadderResoShaped<TSig, TPar>::setSampleRate(newSampleRate);
  resoFilter.setSampleRate(newSampleRate);
}

template<class TSig, class TPar>
void rsResoReplacer<TSig, TPar>::getSignalParts(TSig in, TSig *yf, TSig *yr)
{
  TSig mag, phs;
  getNonresonantOutputAndResonanceParameters(in, yf, &mag, &phs);
  *yr = getWaveForm(mag, phs);
}

template<class TSig, class TPar>
void rsResoReplacer<TSig, TPar>::getNonresonantOutputAndResonanceParameters(TSig in, TSig *yf,
  TSig *mag, TSig *phs)
{
  // temporary variables:
  double r, q;

  *yf = this->nonResonant.getSample(in);
  //r   = this->resonant.getSampleLinear(in) - *yf;  // pure resonance - old code (linear - no growl)
  r   = this->resonant.getSample(in) - *yf;          // pure resonance
  r   = this->attackSmoother.getSample(r);           // attack smoothing filter

  // compute magnitude and phase:
  q     = this->allpass.getSample(r);                // quadrature signal
  *mag  = sqrt(r*r + q*q);                    // extract magnitude
  *phs  = atan2(r, q);                        // exctract phase

  // manipulate phase (later add phase-modulation):
  *phs += this->phase;  // add static offset
  *phs += 2*PI;         // rsSqrWave, etc. expects phase in 0...2*PI - get rid when using a wavetable
}

template<class TSig, class TPar>
TSig rsResoReplacer<TSig, TPar>::getWaveForm(TSig a, TSig p)
{
  double r;

  // resonance waveform reconstruction (prelimiary - later use a mip mapped wavetable, where the
  // mip map level is chosen according the the phase difference between previous and current
  // sample):
  switch(waveform)
  {
  case 1:  r =  a * rsTriWave(p); break;
  case 2:  r =  a * rsSqrWave(p); break;
  case 3:  r = -a * rsSawWave(p); break;
  case 4:  r =  a * rsSawWave(p); break;
  default: r =  a * sin(p);
  }

  // set up the resonance filter's cutoff frequency and filter the resonance waveform:
  TPar ctf = this->cutoff*(resoCutoffMultiplier + resoCutoffModulation*a);
  ctf = rsClip(ctf, 20.0, 20000.0);
  resoFilter.setCutoff(ctf);
  r = resoFilter.getSample(r);
  // todo: optimization: skip setCutoff if resoCutoffModulation == 0 - but we need to take care
  // that it's called in setResoCutoffMultiplier, setResoCutoffModulation, setCutoff

  //r /= resoFilter.getMagnitudeAt(cutoff); // gain correction - maybe not such a good idea

  return r;
}

//-------------------------------------------------------------------------------------------------

template<class TSig, class TPar>
rsResoReplacerPhaseBumped<TSig, TPar>::rsResoReplacerPhaseBumped()
{
  // todo: set up highpass and lowpass filters for the bumping signal

  inputDifferencer.setMode(rsOnePoleFilter<TSig, TPar>::HIGHPASS_MZT);
  bumpSmoother1.setMode(rsOnePoleFilter<TSig, TPar>::LOWPASS_IIT);
  bumpSmoother2.setMode(rsOnePoleFilter<TSig, TPar>::LOWPASS_IIT);


  inputDifferencer.setCutoff(20000);

  // preliminary:
  bumpCutoff1 = 20;
  bumpCutoff2 = 100;
  bumpFactor  = 5.0;         // maybe this should be related to the 2 cutoffs
  bumpSmoother1.setCutoff(bumpCutoff1);
  bumpSmoother2.setCutoff(bumpCutoff2);

  averageSamples = 2000;
  bumpAverager.setLengthInSamples((int)averageSamples);
  bumpAverager.reset();


  setChaosParameter(0.0);

  setAmplitudeLimit(10.0);
  setInputRange(1.0);

  setSelfExitation(0.0);
  setAmplitudeOffset(0.0);

  // move into reset function:
  bump     = 0.0;
  oldPhase = 0.0;
}

template<class TSig, class TPar>
void rsResoReplacerPhaseBumped<TSig, TPar>::setSampleRate(TPar newSampleRate)
{
  rsResoReplacer<TSig, TPar>::setSampleRate(newSampleRate);
  inputDifferencer.setSampleRate(newSampleRate);
  bumpSmoother1.setSampleRate(newSampleRate);
  bumpSmoother2.setSampleRate(newSampleRate);
}

template<class TSig, class TPar>
void rsResoReplacerPhaseBumped<TSig, TPar>::setBumpFactor(TPar newFactor)
{
  bumpFactor = newFactor;
}

template<class TSig, class TPar>
void rsResoReplacerPhaseBumped<TSig, TPar>::setBumpCutoffs(TPar cutoff1, TPar cutoff2)
{
  bumpCutoff1 = cutoff1;
  bumpCutoff2 = cutoff2;
  bumpSmoother1.setCutoff(bumpCutoff1);
  bumpSmoother2.setCutoff(bumpCutoff2);
}

template<class TSig, class TPar>
void rsResoReplacerPhaseBumped<TSig, TPar>::setChaosParameter(TPar newParameter)
{
  chaosParam = newParameter;
}

template<class TSig, class TPar>
void rsResoReplacerPhaseBumped<TSig, TPar>::setAmplitudeLimit(TPar newLimit)
{
  rsAssert(newLimit > 0.0);
  ampLimit = newLimit;
  ampLimitInv = 1.0 / ampLimit;
}

template<class TSig, class TPar>
void rsResoReplacerPhaseBumped<TSig, TPar>::setInputRange(TPar newRange)
{
  inRange = newRange;
}

template<class TSig, class TPar>
void rsResoReplacerPhaseBumped<TSig, TPar>::setSelfExitation(TPar newExcitation)
{
  selfExcite = newExcitation;
}

template<class TSig, class TPar>
void rsResoReplacerPhaseBumped<TSig, TPar>::setAmplitudeOffset(TPar newOffset)
{
  ampOffset = newOffset;
}

template<class TSig, class TPar>
void rsResoReplacerPhaseBumped<TSig, TPar>::getSignalParts(TSig in, TSig *yf, TSig *yr)
{
  TSig ex = 0.0;      // exitation signal
  ex = getExcitation(in); // comment to deactivate self-excitation code

  // copied from baseclass:
  TSig mag, phs;
  this->getNonresonantOutputAndResonanceParameters(in+ex, yf, &mag, &phs);
  phs     += bump;                // new, compared to baseclass (maybe, we need to apply wrapToInterval?)
  mag     += ampOffset;
  mag      = limitMagnitude(mag);   // new - maybe optimize by inlining (2-level)
  *yr      = this->getWaveForm(mag, phs); // resonance reconstruction/replacement

  //*yr *= inputScale(in);        // input based limiting
    // creates too narrow spikes - use lowpassed version instead
  *yr *= inputScale(*yf);         // lowpass based limiting

  updateBump(in, *yf, *yr);            // also new
}

template<class TSig, class TPar>
void rsResoReplacerPhaseBumped<TSig, TPar>::updateBump(TSig in, TSig yf, TSig yr)
{
  TSig tmp;
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
}

template<class TSig, class TPar>
TSig rsResoReplacerPhaseBumped<TSig, TPar>::limitMagnitude(TSig mag)
{
  //return ampLimit * rsPositiveSigmoids::softClip(ampLimitInv*mag);

  return ampLimit * rsPositiveSigmoids<TSig>::quartic(ampLimitInv*mag);
    // maybe try a different function
}

template<class TSig, class TPar>
TSig rsResoReplacerPhaseBumped<TSig, TPar>::inputScale(TSig in)
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

template<class TSig, class TPar>
TSig rsResoReplacerPhaseBumped<TSig, TPar>::getExcitation(TSig in)
{
  // create self-excitation signal (very experimental - can be optimized):

  double omega = 2*PI*this->cutoff/this->sampleRate; // optimize: precompute in setCutoff

  oldPhase += omega;
  oldPhase  = fmod(oldPhase, 2*PI);

  //// experimental (phase sync/clamp):
  //if(in > inRange)
  //  oldPhase = 0.5*PI;
  //else if(in < -inRange)
  //  oldPhase = 1.5*PI;

  return selfExcite * sin(oldPhase);


  // todo: give the function a parameter: in (the input signal) such that we can apply input
  // dependent phase-manipulations here (such as sync, clamp, etc.)
  // clamp: when input is above high-threshold: phase is fixed to pi, below low-threshold: phase
  // fixed to -pi (should simulate the analog (supposed) phase-sticking effect)

  // maybe have also a self-excitation phase-offset (user parameter) and maybe we should add
  // "omega" to the oldPhase because the oldPhase is, well, old (lagging one omega behind)
}

//-------------------------------------------------------------------------------------------------

template<class TSig, class TPar>
rsLadderMystran<TSig, TPar>::rsLadderMystran()
{
  this->setMode(this->LP_24);
  reset();
}

template<class TSig, class TPar>
void rsLadderMystran<TSig, TPar>::reset()
{
  rsArrayTools::fillWithZeros(s, 4);
  zi = 0.0;
}
