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


/*
Ideas:
-let the user insert different types of mode filters (simple decaying sines, attack/decay-sines, 
 4-env-sines, nonlinear modes (perhaps with amplitude dependent frequency - they should have a 
 second "sidechain" input where we feed back the total summed output - so the nonlinear effects may
 depend on the total output value


*/

