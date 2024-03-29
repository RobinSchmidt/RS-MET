
Enveloper::Enveloper(CriticalSection *lockToUse) 
  : AudioModuleWithMidiIn(lockToUse) /* , envGenWrapper(lockToUse, &envGen)*/
{
  ScopedLock scopedLock(*lock);
  setModuleTypeName("Enveloper");
  envGenWrapper = new BreakpointModulatorAudioModule(lockToUse, &envGen);
  addChildAudioModule(envGenWrapper);
}

AudioModuleEditor* Enveloper::createEditor(int type)
{
  //jura::BreakpointModulatorEditor* editor = new jura::BreakpointModulatorEditor(lock, &envGenWrapper);
  jura::BreakpointModulatorEditor* editor = new jura::BreakpointModulatorEditor(lock, envGenWrapper);

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

XmlElement* Enveloper::getStateAsXml(const juce::String& stateName, bool markAsClean)
{
  return envGenWrapper->getStateAsXml(stateName, markAsClean);
}

void Enveloper::setStateFromXml(const XmlElement& xmlState, const juce::String& stateName,
  bool markAsClean)
{
  envGenWrapper->setStateFromXml(xmlState, stateName, markAsClean);
}

//void Enveloper::reset()
//{
// 
//}

