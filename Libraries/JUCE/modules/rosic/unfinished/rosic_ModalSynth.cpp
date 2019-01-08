using namespace rosic;

void rsModalFrequencyGenerator::harmonic(double* r, int N)
{
  for(int i = 0; i < N; i++)
    r[i] = double(i+1);
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
}

/*
void rsModalFrequencyGenerator::idealBar(double* r, int N)
{

}
*/

void rsModalFrequencyGenerator::stiffString(double* r, int N, double B)
{
  for(int i = 0; i < N; i++) {
    double n = double(i+1);
    r[i] = n*sqrt(1+B*n*n);
  }
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

void rsModalSynth::setFreqRatioProfile1(int newProfile)
{
  freqRatioProfile1 = newProfile;
  fillFreqRatios(freqRatios1, freqRatiosLog1, freqRatioProfile1);
}

void rsModalSynth::setFreqRatioProfile2(int newProfile)
{
  freqRatioProfile2 = newProfile;
  fillFreqRatios(freqRatios2, freqRatiosLog2, freqRatioProfile2);
}

void rsModalSynth::setFreqRatioProfile3(int newProfile)
{
  freqRatioProfile3 = newProfile;
  fillFreqRatios(freqRatios3, freqRatiosLog3, freqRatioProfile3);
}

void rsModalSynth::setFreqRatioProfile4(int newProfile)
{
  freqRatioProfile4 = newProfile;
  fillFreqRatios(freqRatios4, freqRatiosLog4, freqRatioProfile4);
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


void rsModalSynth::fillFreqRatios(double* ratios, double *logRatios, int profile)
{
  typedef rsModalFrequencyGenerator MFG;
  int N = maxNumModes;
  switch(profile)
  {
  case HARMONIC:     MFG::harmonic(   ratios, N);                break;

  //case STIFF_STRING: MFG::stiffString(ratios, N, inharmonicity); break;
  default: MFG::harmonic(ratios, N);
  };
  RAPT::rsArray::applyFunction(ratios, logRatios, maxNumModes, &log);
  updateFreqRatios();
}

void rsModalSynth::noteOn(int key, int vel)
{
  // todo: update the freqRatios according to key/vel - these should depend on key/vel


  phaseGenerator.setSeed(phaseRandomSeed);
  phaseGenerator.setRange(-phaseRandomness*PI, phaseRandomness*PI);

  double f0 = RAPT::rsPitchToFreq(key);
  double r, f;
  double w, A, p, att, dec;
  double ampSlopeExponent = -1;  // preliminary, -1 should result from a setting of -6.02 dB/oct
  int m = 0;
  for(m = 0; m < numPartialsLimit; m++) {
    r = freqRatios[m];
    f = f0 * r;
    if( f > 0.5*sampleRate )
      break;

    // compute modal filter algorithm parameters (preliminary) and set up the filter:
    w = 2 * PI * f / sampleRate;              // optimize - precompute 2*PI/fs
    A = amplitude * pow(r, ampSlopeExponent); // maybe this can be expressed via exp?
    p = phaseGenerator.getSample();
    att = 0.001 * attack * sampleRate;
    dec = 0.001 * decay  * sampleRate;
    modalBank.setModalFilterParameters(m, w, A, p, att, dec
      /*,deltaOmega, phaseDelta, blend, attackScale, decayScale*/);
  }

  modalBank.setNumModes(m); // is this correct? or off by one? check this!
  modalBank.reset();        // maybe allow a partial reset
  noteAge = 0;
}

/*
void rsModalSynth::fillFreqRatiosHarmonic(double* ratios)
{
  for(int i = 0; i < maxNumModes; i++)
    ratios[i] = double(i+1);
}

void rsModalSynth::fillFreqRatiosStiffString(double* ratios, double B)
{
  for(int i = 0; i < maxNumModes; i++) {
    double n = double(i+1);
    ratios[i] = n*sqrt(1+B*n*n);
  }
}
*/

void rsModalSynth::updateFreqRatios()
{
  freqRatiosAreReady = false;

  // obtain frequency ratios by bilinear interpolation either of the frequencies themselves or
  // the associated pitches:
  double rt, rb, r;  // ratio for top and bottom and final result
  if(logLinearFreqInterpolation) // maybe rename to pitchInterpolation
    for(int i = 0; i < maxNumModes; i++) {
      rt = (1-freqRatioMixX)*freqRatiosLog1[i] + freqRatioMixX*freqRatiosLog2[i];
      rb = (1-freqRatioMixX)*freqRatiosLog3[i] + freqRatioMixX*freqRatiosLog4[i];
      r  = (1-freqRatioMixY)*rt + freqRatioMixY*rb;
      freqRatios[i] = exp(r);
    }
  else
    for(int i = 0; i < maxNumModes; i++) {
      rt = (1-freqRatioMixX)*freqRatios1[i] + freqRatioMixX*freqRatios2[i];
      rb = (1-freqRatioMixX)*freqRatios3[i] + freqRatioMixX*freqRatios4[i];
      r  = (1-freqRatioMixY)*rt + freqRatioMixY*rb;
      freqRatios[i] = r;
    }

  freqRatiosAreReady = true;
}


/*
Ideas:
-let the user insert different types of mode filters (simple decaying sines, attack/decay-sines, 
 4-env-sines, nonlinear modes (perhaps with amplitude dependent frequency - they should have a 
 second "sidechain" input where we feed back the total summed output - so the nonlinear effects may
 depend on the total output value


 a lot of code can be dragged in from ModalExample.h/cpp in rosic_tests/PortedFromRSLib/Examples

*/

