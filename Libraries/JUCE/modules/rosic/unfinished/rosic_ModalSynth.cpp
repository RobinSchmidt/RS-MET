using namespace rosic;


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
  switch(profile)
  {
  case HARMONIC:     fillFreqRatiosHarmonic(ratios);                   break;
  case STIFF_STRING: fillFreqRatiosStiffString(ratios, inharmonicity); break;
  default: fillFreqRatiosHarmonic(ratios);
  };
  RAPT::rsArray::applyFunction(ratios, logRatios, maxNumModes, &log);
  updateFreqRatios();
}

void rsModalSynth::noteOn(int key, int vel)
{
  // todo: update the freqRatios according to key/vel

  // todo: update the modal filter bank - recompute filter coeffs, reset states,...

  double f0 = rsPitchToFreq(key);
  double r, f;
  for(int m = 0; m < numPartialsLimit; m++) {
    r = freqRatios[m];
    f = f0 * r;


  }

  noteAge = 0;
}

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