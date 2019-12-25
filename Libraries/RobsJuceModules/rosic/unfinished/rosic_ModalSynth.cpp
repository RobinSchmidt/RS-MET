using namespace rosic;

void rsModalFrequencyGenerator::allHarmonics(double* r, int N)
{
  for(int i = 0; i < N; i++)
    r[i] = double(i+1);
}

void rsModalFrequencyGenerator::oddHarmonics(double* r, int N)
{
  for(int i = 0; i < N; i++)
    r[i] = double(2*i+1);
}

void rsModalFrequencyGenerator::twelveTone(double* r, int N)
{
  twelveTone21(r, N);
  long double s = pow(2.0, 1.0/12.0); // basis
  for(int n = 21; n < N; n++)
    r[n] = pow(s, 54+n-21);

  // maybe use a formula:
  // b  = pow(2.0, 1.0/12.0); // basis - can be generalized
  // fn = b^k where k = round(logB(n*f0, b)...or maybe use a special rounding:

  // roundedExponent(target, basis)
  // a  = logB(target, basis);
  // af = floor(a);
  // ac = ceil(a);
  // yf = pow(basis, yf);
  // yc = pow(basis, yc);
  // if(fabs(target-yc) < fabs(target-yf))
  //   return ac;
  // else
  //   return af;
}

void rsModalFrequencyGenerator::twelveTonePseudoHarmonic(double* r, int N)
{
  twelveTone21(r, N);
  for(int i = 21; i < N; i++)
    r[i] = double(i+1);
}

void rsModalFrequencyGenerator::rodFreeFree(double* r, int N)
{
  RAPT::rsAssert(N >= 6);
  r[0] = 4.7300407448627040260240481;
  r[1] = 7.8532046240958375564770667;
  r[2] = 10.9956078380016709066690325;
  r[3] = 14.1371654912574641771059179;
  r[4] = 17.2787596573994814380910740;
  r[5] = 20.4203522456260610909364112;
  double m, c;
  int i;
  for(i = 6; i < N; i++) {
    m = i+1;
    c = (m+0.5)*PI;
    r[i] = c - pow(-1, m) * 2 * exp(-c) - 4*exp(-2*c);
  }
  c = r[0]*r[0];
  for(i = 0; i < N; i++)
    r[i] = r[i]*r[i] / c;
  // ...verify and optimize....
  // Reference: Music: A Mathematical Offering (Dave Benson), page 119 
}

void rsModalFrequencyGenerator::rodFreeClamped(double* r, int N)
{
  RAPT::rsAssert(N >= 6);
  r[0] = 1.8751040687119611664453082;
  r[1] = 4.6940911329741745764363918;
  r[2] = 7.8547574382376125648610086;
  r[3] = 10.9955407348754669906673491;
  r[4] = 14.1371683910464705809170468;
  r[5] = 17.2787595320882363335439284;
  double m, c;
  int i;
  for(i = 6; i < N; i++) {
    m = i+1;
    c = (m-0.5)*PI;
    r[i] = c - pow(-1, m) * 2 * exp(-c) - 4 * exp(-2*c); 
  }
  c = r[0]*r[0];
  for(i = 0; i < N; i++)
    r[i] = r[i]*r[i] / c;
  // ....verify and optimize...
  // Reference: Music: A Mathematical Offering (Dave Benson), page 123 
}

/*
void rsModalFrequencyGenerator::idealBar(double* r, int N)
{

}
*/

void rsModalFrequencyGenerator::stiffString(double* r, int N, double B)
{
  // todo: use RAPT::rsStiffStringFreqRatio
  for(int i = 0; i < N; i++) {
    double n = double(i+1);
    r[i] = n*sqrt(1+B*n*n);
  }
  // see:
  // http://www.simonhendry.co.uk/wp/wp-content/uploads/2012/08/inharmonicity.pdf  Eq.10
  // http://www.jbsand.dk/div/StivStreng.pdf
}

void rsModalFrequencyGenerator::twelveTone21(double* r, int N)
{
  double tmp[21];
  long double s = pow(2.0, 1.0/12.0); // basis

                         //  #    ratio
  tmp[0]  = pow(s,  0);  //  1    1.0
  tmp[1]  = pow(s, 12);  //  2    2.0
  tmp[2]  = pow(s, 19);  //  3    2.9966141537533639
  tmp[3]  = 4.0;         //  4    4.0
  tmp[4]  = pow(s, 28);  //  5    5.0396841995794937
  tmp[5]  = pow(s, 31);  //  6    5.9932283075067279
  tmp[6]  = pow(s, 34);  //  7    7.1271897451227169
  tmp[7]  = 8.0;         //  8    8.0
  tmp[8]  = pow(s, 38);  //  9    8.9796963864749877
  tmp[9]  = pow(s, 40);  // 10   10.079368399158989
  tmp[10] = pow(s, 42);  // 11   11.313708498984765
  tmp[11] = pow(s, 43);  // 12   11.986456615013459
  tmp[12] = pow(s, 44);  // 13   12.699208415745600
  tmp[13] = pow(s, 46);  // 14   14.254379490245437
  tmp[14] = pow(s, 47);  // 15   15.101989002907104
  tmp[15] = 16.0;        // 16   16.0  
  tmp[16] = pow(s, 49);  // 17   16.951409509748732
  tmp[17] = pow(s, 50);  // 18   17.959392772949979
  tmp[18] = pow(s, 51);  // 19   19.027313840043551
  tmp[19] = pow(s, 52);  // 20   20.158736798317982
  tmp[20] = pow(s, 53);  // 21   21.357437666720561
  // optimize - get rid of calling pow - use direct constants (or can the compiler do that?)

  for(int i = 0; i < RAPT::rsMin(N, 21); i++)
    r[i] = tmp[i];
}

//=================================================================================================






/*
void rsModalAlgoParameters::setFromUserParameters(const rsModalUserParameters& usrPars, double fs)
{
  w      = 2 * PI * usrPars.frequency / fs;
  A      = usrPars.amplitude;
  p      = RAPT::rsDegreeToRadiant(usrPars.phase);
  att    = 0.001 * usrPars.attack * fs;
  dec    = 0.001 * usrPars.decay  * fs;
  dw     = 2 * PI * usrPars.freqSpread / fs;
  dp     = RAPT::rsDegreeToRadiant(usrPars.phaseDelta);
  b      = usrPars.blend;
  attScl = usrPars.attackScale;
  decScl = usrPars.decayScale;
}
*/

//=================================================================================================

rsModalSynth::rsModalSynth()
{
  setFreqRatioProfileTopLeft(ALL_HARMONICS);
  setFreqRatioProfileTopRight(ALL_HARMONICS);
  setFreqRatioProfileBottomLeft(ALL_HARMONICS);
  setFreqRatioProfileBottomRight(ALL_HARMONICS);
}

void rsModalSynth::setFreqRatioProfileTopLeft(int newProfile)
{
  freqRatioProfileTopLeft = newProfile;
  fillFreqRatios(freqRatiosTopLeft, freqRatiosLogTopLeft, freqRatioProfileTopLeft);
}

void rsModalSynth::setFreqRatioProfileTopRight(int newProfile)
{
  freqRatioProfileTopRight = newProfile;
  fillFreqRatios(freqRatiosTopRight, freqRatiosLogTopRight, freqRatioProfileTopRight);
}

void rsModalSynth::setFreqRatioProfileBottomLeft(int newProfile)
{
  freqRatioProfileBottomLeft = newProfile;
  fillFreqRatios(freqRatiosBottomLeft, freqRatiosLogBottomLeft, freqRatioProfileBottomLeft);
}

void rsModalSynth::setFreqRatioProfileBottomRight(int newProfile)
{
  freqRatioProfileBottomRight = newProfile;
  fillFreqRatios(freqRatiosBottomRight, freqRatiosLogBottomRight, freqRatioProfileBottomRight);
}

void rsModalSynth::setFreqRatioMixX(double newMix)
{
  freqRatioMixX = newMix;
  updateFreqRatios();
}

void rsModalSynth::setFreqRatioMixY(double newMix)
{
  freqRatioMixY = newMix;
  updateFreqRatios();
}

//-------------------------------------------------------------------------------------------------

double rsModalSynth::getModeFreqRatio(int i, double k, double v) const
{
  return freqRatios[i];
}

double rsModalSynth::getModeFrequency(int i, double k, double v) const
{
  return RAPT::rsPitchToFreq(k) * getModeFreqRatio(i, k, v);
}

double rsModalSynth::getModeRelativePitch(int i, double k, double v) const
{
  double r = getModeFreqRatio(i, k, v);
  return FREQFACTOR2PITCHOFFSET(r); // move as function to rapt
  //return RAPT::rsFreqFactorToPitchOffset(r);
  // todo: (maybe) store mode-pitches instead of the natural log - we'll see which one is more
  // convenient in terms of writing optimized code - i think, we could return here
  // something * freqRatiosLog[i]
}

double rsModalSynth::getModeAmplitude(int i, double k, double v) const
{
  double dB = getModeLevel(i, k, v);
  return RAPT::rsDbToAmp(dB);
}

double rsModalSynth::getModeLevel(int i, double k, double v) const
{
  double dB   = level + levelByKey*(k-64)/63. + levelByVel*(v-64)/63.;
  double octs = getModeRelativePitch(i, k, v) / 12.0; // how many octaves is the mode above f0
  double kk = ampSlopeByKey*(k-64)/63.;
  double kv = ampSlopeByVel*(v-64)/63.;
  dB += octs * (ampSlope + kk + kv);
  return dB;
}

double rsModalSynth::getModeAttack(int i, double k, double v) const
{
  double cr  = getTimeCoeffByIndex(i);
  double ck  = getTimeCoeffByKey(k);
  double cv  = getTimeCoeffByVel(v);
  double att = attack * exp(cr*attackByRatio  + ck*attackByKey  + cv*attackByVel);
  return att;
}

double rsModalSynth::getModeDecay(int i, double k, double v) const
{
  double cr  = getTimeCoeffByIndex(i);
  double ck  = getTimeCoeffByKey(k);
  double cv  = getTimeCoeffByVel(v);
  double dec = decay * exp(cr*decayByRatio  + ck*decayByKey  + cv*decayByVel);
  return dec;
}

double rsModalSynth::getTimeCoeffByIndex(int i) const
{
  double cr = freqRatiosLog[i];
  return cr;
}

double rsModalSynth::getTimeCoeffByKey(double key) const
{
  double ck = (LN2/12.0) * (key-64.0); // 64 is the neutral reference key
  return ck;
}

double rsModalSynth::getTimeCoeffByVel(double vel) const
{
  double cv = (LN2/63.0) * (vel-64.0); // 64 is the neutral reference velocity
  return cv;
}
// i think, for scaling the attack/decay time-values by key, we should use the scaler
// 2^( decByKey * (key - refKey) / 12 )
// if decByKey = 1, this would scale the time up by factor 2 when key is and octave above
// the reference key (which should probably be 64. this translates to:
// exp( log(2) * decByKey * (key - refKey) / 12 ) = exp(ck*decByKey)
// i think., the cv formula would make times twice as long at full velocity 127 and half
// as long at minimum nonzero velocity 1 ...check this

//-------------------------------------------------------------------------------------------------

void rsModalSynth::fillFreqRatios(double* ratios, double *logRatios, int profile)
{
  typedef rsModalFrequencyGenerator MFG;
  int N = maxNumModes;
  switch(profile)
  {
  case ALL_HARMONICS:               MFG::allHarmonics(            ratios, N); break;
  case ODD_HARMONICS:               MFG::oddHarmonics(            ratios, N); break;
  case TWELVE_TONE_PSEUDO_HARMONIC: MFG::twelveTonePseudoHarmonic(ratios, N); break;
  case TWELVE_TONE_EQUAL:           MFG::twelveTone(              ratios, N); break;
  case ROD_FREE_FREE:               MFG::rodFreeFree(             ratios, N); break;
  case ROD_FREE_CLAMPED:            MFG::rodFreeClamped(          ratios, N); break;

  //case STIFF_STRING: MFG::stiffString(ratios, N, inharmonicity); break;

  default: MFG::allHarmonics(ratios, N);
  };
  RAPT::rsArrayTools::applyFunction(ratios, logRatios, maxNumModes, &log);
  updateFreqRatios();
}

void rsModalSynth::noteOn(double key, double vel)
{
  if(vel == 0)
    return; // note off

  // todo: update the freqRatios according to key/vel - these should depend on key/vel

  double ck = getTimeCoeffByKey(key);  // maybe rename to ckt - or can the same coeffs be used for
  double cv = getTimeCoeffByVel(vel);  // other things too (amplitudes)?

  phaseGenerator.setSeed(phaseRandomSeed);
  phaseGenerator.setRange(-phaseRandomness*PI, phaseRandomness*PI);

  double f0 = RAPT::rsPitchToFreq(key);
  double r, f;
  double w, A, p, att, dec;
  int m = 0;

  for(m = lowestMode-1; m < highestMode; m++) {

    r = freqRatios[m];
    f = f0 * r;
    if( f > 0.5*sampleRate )
      break;

    // compute modal radian frequency and phase:
    w = 2 * PI * f / sampleRate;              // optimize - precompute 2*PI/fs
    p = phaseGenerator.getSample();

    // compute modal amplitude:
    A = getModeAmplitude(m, key, vel);
    // calling this function here is HORRIBLY inefficient - some things can be precomputed outside
    // the loop and a lot of back-and-forth between linear and logarithmic representations of
    // frequency/pitch and amplitude/decibels should be consolidated into a single call to exp 
    // (if possible - or at least, less back-and forth), similarly to the envelope time 
    // calcluations below - optimize!

    // compute modal envelope parameters:
    double cr = getTimeCoeffByIndex(m);  // optimize later ...or inline computation
    double ts = 0.001 * sampleRate;      // time-scaler - take out of the loop
    att = ts * attack * exp(cr*attackByRatio + ck*attackByKey + cv*attackByVel);
    dec = ts * decay  * exp(cr*decayByRatio  + ck*decayByKey  + cv*decayByVel);


    modalBank.setModalFilterParameters(m, w, A, p, att, dec
      /*,deltaOmega, phaseDelta, blend, attackScale, decayScale*/);
  }

  modalBank.setLowModeIndex( lowestMode-1); // modalBank counts modes from 0, rsModalSynth from 1
  //modalBank.setHighModeIndex(highestMode-1);
  //modalBank.setHighModeIndex(RAPT::rsMin(m, highestMode-1)); // check this - correct or off by 1?
  modalBank.setHighModeIndex(RAPT::rsMin(m, highestMode)-1);  // ...this looks better
  noteAge = 0;

  modalBank.reset();        // maybe allow a partial reset
  // i think, i in this note more modes become (re)activated than in the previous node, then 
  // without reset, the new modes will start from an old state (when they were active the last 
  // time) - this is not good - maybe we should reset modes when they are turned off, i.e.
  // loop from highestMode to maxNumModes and call reset for each - or something
}

void rsModalSynth::updateFreqRatios()
{
  freqRatiosAreReady = false;

  // obtain frequency ratios by bilinear interpolation of the associated pitches:
  double rt, rb, r;  // ratio for top and bottom and final result


  double mixX = 0.5 * (freqRatioMixX + 1.0);
  double mixY = 0.5 * (freqRatioMixY + 1.0);

  for(int i = 0; i < maxNumModes; i++) {
    rt = (1-mixX)*freqRatiosLogTopLeft[i]    + mixX*freqRatiosLogTopRight[i];
    rb = (1-mixX)*freqRatiosLogBottomLeft[i] + mixX*freqRatiosLogBottomRight[i];
    //r  = (1-mixY)*rt + mixY*rb;
    r  = (1-mixY)*rb + mixY*rt;
    freqRatiosLog[i] = r;
    freqRatios[i] = exp(r);
  }

  freqRatiosAreReady = true;
}


/*
ToDo:
-implement ratio/key/vel tracking for amp/attack/decay
-implement attackScale / decayScale / Blend (blend also by ratio/key/vel)
-allow negative mixing, i.e. subtraction of the 2nd pair of filters -> let the mixing coeff for
 the two pairs be outside the 0..1 range
-to implement different mode ratios at different keys, render two sets of mode-ratios, one for key 
 0 (or 1) and one for key 127 and interpolate according to the played key - at key=64 we should 
 get the unmodified mode ratios
 -it seems impractical to recompute all mode-ratios from the underlying formulas each time a new 
  key is struck - this way, we can efficiently implement keyscaling of the mode-ratios
 -maybe do a similar thing for velocity scaling (but key-scaling first - it's more important)
 -this would mean that in the end, we have a twofold bilinear interpolation. the first according
  to the vector pad, the second according to key/vel

Ideas:
-actually, it would be a good opportunity to insert tah as module into the NewSynth
-let the user insert different types of mode filters (simple decaying sines, attack/decay-sines, 
 4-env-sines, nonlinear modes (perhaps with amplitude dependent frequency - they should have a 
 second "sidechain" input where we feed back the total summed output - so the nonlinear effects may
 depend on the total output value

-we need something like a highpass with keytracking or some sort of additional highpassish 
 amplitude weighting function

 a lot of code can be dragged in from ModalExample.h/cpp in rosic_tests/PortedFromRSLib/Examples

-look at the straightliner osc fo inspiration for more parameters for ampltude and pahse spectrum
 for example, use LowestHarmonic/HighestHarmonic



Observations:
-with Attack = 50ms and Decay = 1000ms with -100% byRatio, there are strange comb-like effects
 -maybe it's because at the higher frequencies the attack (which has no ratio-tracking) actually
  becomes longer than the Decay such that attck/decay sort of reverse their roles?


see here for how to figure out mode decay times from samples:
  https://ccrma.stanford.edu/~jos/pasp/Sinusoidal_modeling.html
*/

