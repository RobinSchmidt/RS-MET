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
  RAPT::rsArray::applyFunction(ratios, logRatios, maxNumModes, &log);
  updateFreqRatios();
}

void rsModalSynth::noteOn(int key, int vel)
{
  if(vel == 0)
    return; // note off

  // todo: update the freqRatios according to key/vel - these should depend on key/vel

  // i think, for scaling the attack/decay time-values by key, we should use the scaler
  // 2^( decByKey * (key - refKey) / 12 )
  // if decByKey = 1, this would scale the time up by factor 2 when key is and octave above
  // the reference key (which should probably be 64. this translates to:
  // exp( log(2) * decByKey * (key - refKey) / 12 ) = exp(ck*decByKey)
  // with ck = log(2)*(key-refKey)/12
  //int refKey = 64;
  double ck = (LN2/12.0) * double(key-64);  // 64 is the neutral reference key
  double cv = (LN2/63.0) * double(vel-64);  // 64 is the neutral reference velocity
  //cv = 0; // preliminary
  // i think., the cv formula would make time stwice as long at full velocity 127 and half
  // as long at minimum nonzero velocity 1 ...check this



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

    double rLog = freqRatiosLog[m]; // maybe rename to cr for consistency with ck, cv

    // compute modal filter algorithm parameters (preliminary) and set up the filter:
    w = 2 * PI * f / sampleRate;              // optimize - precompute 2*PI/fs
    p = phaseGenerator.getSample();

    //A = amplitude * pow(r, ampSlopeExponent); // maybe this can be expressed via exp?
    //A = amp * pow(r, ampByRatio); 
    A = amp * exp(rLog * ampByRatio);  // should be the same - yes - works

    // actually, it doesn't semm to be enough to scale amplitude by ratio, key and vel - we
    // also need to scale the ampByRatio itself by key and vel - which is not the same thing

    // i think, we should have Level with Key/Vel and Slope also with Key/Vel

    att = 0.001 * attack * sampleRate;
    dec = 0.001 * decay  * sampleRate;

    // this is more efficient and extensible at negligible cost by adding more scalers later:
    att *= exp(rLog*attackByRatio + ck*attackByKey + cv*attackByVel);
    dec *= exp(rLog*decayByRatio  + ck*decayByKey  + cv*decayByVel);

    // scale attack/decay according to key, vel and - importantly - r

    //att *= pow(r, attackByRatio);
    //dec *= pow(r, decayByRatio);

    // to optimize this, use pow(r, x) = exp(x*log(r)) - log r is computed anyway in 
    // updateFreqRatios - we should store it in a member array there and should get completely rid
    // of linear interpolation of mode frequency (it doesn't seem to be useful anyway) - do it 
    // always in logarithmic (pitch) domain as we need the log of the ratio here anyway
    // at the end, all the various scalings (by key, vel, ratio) should be incorporated into a 
    // weigthed sum and a single call to exp should be done, like this:
    // dec *= exp(c1*decByRatio*log(r) + c2*decByKey*(key-refKey) + c3*decByVel*(vel-refVel))
    // here c1,c2,c3 are appropriate fixed constants (maybe log(2)? - we'll see - check how 
    // straightliner handles key and vel scaling and use the same rule here)
    // OK - the optimization is implemented below - todo add the pther scalers...

    // test:
    //double s1, s2;
    //s1 = pow(r, decayByRatio);
    //s2 = exp(rLog*decayByRatio);  // should be the same - yep, works


    modalBank.setModalFilterParameters(m, w, A, p, att, dec
      /*,deltaOmega, phaseDelta, blend, attackScale, decayScale*/);
  }

  modalBank.setNumModes(m); // is this correct? or off by one? check this!
  //modalBank.reset();        // maybe allow a partial reset
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

Ideas:
-actually, it would be a good opportunity to insert tah as module into the NewSynth
-let the user insert different types of mode filters (simple decaying sines, attack/decay-sines, 
 4-env-sines, nonlinear modes (perhaps with amplitude dependent frequency - they should have a 
 second "sidechain" input where we feed back the total summed output - so the nonlinear effects may
 depend on the total output value

-we need something like a highpass with keytracking or some sort of additional highpassish 
 amplitude weighting function

 a lot of code can be dragged in from ModalExample.h/cpp in rosic_tests/PortedFromRSLib/Examples


Observations:
-with Attack = 50ms and Decay = 1000ms with -100% byRatio, there are strange comb-like effects
 -maybe it's because at the higher frequencies the attack (which has no ratio-tracking) actually
  becomes longer than the Decay such that attck/decay sort of reverse their roles?
*/

