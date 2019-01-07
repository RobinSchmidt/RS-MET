using namespace rosic;

void rsModalSynth::setFreqRatioProfile(int newProfile)
{

}

void rsModalSynth::setFreqRatioAdjustment(double newAdjustmentAmount)
{

}

void rsModalSynth::updateFreqRatios()
{
  freqRatiosAreReady = false;

  // todo: switch between freqRatioProfiles

  for(int i = 0; i < rsModalBankFloatSSE2::maxNumModes; i++)
  {
    int n = i+1;
    freqRatios[i] = double(i+1);
  }

  freqRatiosAreReady = true;
}