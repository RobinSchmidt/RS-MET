using namespace rosic;

void rsModalSynth::setFreqRatioProfile1(int newProfile)
{
  freqRatioProfile1 = newProfile;
  fillFreqRatios(freqRatios1, freqRatioProfile1);
}

void rsModalSynth::setFreqRatioProfile2(int newProfile)
{
  freqRatioProfile2 = newProfile;
  fillFreqRatios(freqRatios2, freqRatioProfile2);
}

void rsModalSynth::setFreqRatioMix(double newMix)
{
  freqRatioMix = newMix;
  updateFreqRatios();
}

void rsModalSynth::fillFreqRatios(double* ratios, int profile)
{
  switch(profile)
  {
  case HARMONIC:     fillFreqRatiosHarmonic(ratios);                   break;
  case STIFF_STRING: fillFreqRatiosStiffString(ratios, inharmonicity); break;
  default: fillFreqRatiosHarmonic(ratios);
  };
  updateFreqRatios();
}

void rsModalSynth::noteOn(int key, int vel)
{
  // todo: update the modal filter bank - recompute filter coeffs, reset states,...

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

  if(logLinearFreqInterpolation) // maybe rename to pitchInterpolation
    for(int i = 0; i < maxNumModes; i++)
      freqRatios[i] = exp((1-freqRatioMix)*log(freqRatios1[i]) + freqRatioMix*log(freqRatios2[i]));
  else
    for(int i = 0; i < maxNumModes; i++)
      freqRatios[i] = (1-freqRatioMix)*freqRatios1[i] + freqRatioMix*freqRatios2[i];

  freqRatiosAreReady = true;
}