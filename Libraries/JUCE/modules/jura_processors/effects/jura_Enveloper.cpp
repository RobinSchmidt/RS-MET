
Enveloper::Enveloper(CriticalSection *lockToUse) 
  : AudioModuleWithMidiIn(lockToUse)
  , envGenWrapper(lockToUse, &envGen)
{
  ScopedLock scopedLock(*plugInLock);
  moduleName = "Enveloper";
  setActiveDirectory(getApplicationDirectory() + "/EnveloperPresets");
}

AudioModuleEditor* Enveloper::createEditor()
{
  jura::BreakpointModulatorEditor* editor = 
    new jura::BreakpointModulatorEditor(plugInLock, &envGenWrapper);
  //editor->setLayout(1);
  editor->setSize(500, 260);
  return editor;
  // somehow, the widgets below the plot are messed up ..or, there is something in the background
  // also, the arrangement of the subsections is not yet optimal - we need to get rid of some 
  // margins
}

void Enveloper::processBlock(double **inOutBuffer, int numChannels, int numSamples)
{
  jassert(numChannels == 2);
  for(int n = 0; n < numSamples; n++)
  {
    double env = envGen.getSample();
    inOutBuffer[0][n] *= env;
    inOutBuffer[1][n] *= env;
  }
}

void Enveloper::setSampleRate(double newSampleRate)
{
  envGen.setSampleRate(newSampleRate); 
}

void Enveloper::noteOn(int noteNumber, int velocity)
{
  envGen.noteOn(true, noteNumber, velocity);
}

void Enveloper::noteOff(int noteNumber)
{
  envGen.noteOff(true);
}

//void Enveloper::reset()
//{
// 
//}

