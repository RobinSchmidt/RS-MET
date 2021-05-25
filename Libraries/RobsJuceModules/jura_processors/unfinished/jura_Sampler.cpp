
SamplerModule::SamplerModule(CriticalSection *lockToUse) : AudioModule(lockToUse)
{
  ScopedLock scopedLock(*lock);
  setModuleTypeName("Sampler");
  setModuleName("Sampler");
}

void SamplerModule::createParameters()
{
  ScopedLock scopedLock(*lock);

  typedef SamplerModule SM;
  typedef Parameter Param;
  Param* p;

  p = new Param("Gain", -48.0, 12.0, 0.1, Parameter::LINEAR); 
  addObservedParameter(p);
  p->setValueChangeCallback<SM>(this, &SM::setGain);
}

AudioModuleEditor* SamplerModule::createEditor(int type)
{
  return new jura::SamplerEditor(this);
}

void SamplerModule::setSampleRate(double newSampleRate)
{ 
  engine.setSampleRate(newSampleRate); 
}

void SamplerModule::setGain(double newGain)
{
  //engine.setGlobalGain(newGain);
}

void SamplerModule::setStateFromXml(const XmlElement& xmlState, const juce::String& stateName,
  bool markAsClean)
{
  // Recall the global playback parameters such as gain, max num layers, max polyphony, resampling 
  // quality, etc.:
  AudioModule::setStateFromXml(xmlState, stateName, markAsClean);

  // The actual instrument definition is loaded from an sfz file that is defined in the xml (just 
  // like sample-files are defined in the xml for the wavetable oscillator):
  juce::String jSfzPath = xmlState.getStringAttribute("InstrumentFile", juce::String());
  std::string  sSfzPath = jSfzPath.toStdString();
  engine.loadFromSFZ(sSfzPath.c_str());
  // doesn't work yet - probably because the file is not in the directory where the engine expects 
  // it to be. ToDo: We need to set up the sfz directory in the engine and then the engine must 
  // make use of it

  // ToDo:
  // -obtain sfz file that should be loaded from the xml
  // -retrieve the sample folder ...either from the xml or from sfz...not sure yet, what's best
  // -set up the sample folder in the core
  // -instruct the engine to load the instrument definition from the sfz file
}

XmlElement* SamplerModule::getStateAsXml(const juce::String& stateName, bool markAsClean)
{
  XmlElement *xmlState = AudioModule::getStateAsXml(stateName, markAsClean);
  juce::String sfzPath = sfzFile.getRelativePathFrom(getPresetDirectory());
  xmlState->setAttribute("InstrumentFile", sfzPath);
  //xmlState->setAttribute("SampleDirectory", sampleDir);
  return xmlState;
}

void SamplerModule::processBlock(double **inOutBuffer, int numChannels, int numSamples)
{
  jassert(numChannels == 2);
  for(int n = 0; n < numSamples; n++)
    engine.processFrame(&inOutBuffer[0][n], &inOutBuffer[1][n]);
}

void SamplerModule::processStereoFrame(double *dL, double *dR)
{
  engine.processFrame(dL, dR);
}

void SamplerModule::reset()
{ 
  engine.reset();
}

//=================================================================================================

SamplerEditor::SamplerEditor(SamplerModule* samplerToEdit)
  //: AudioModuleEditor(samplerToEdit->lock) 
  : AudioModuleEditor(samplerToEdit)
{
  ScopedLock scopedLock(*lock);

  // ...

  setSize(400, 200);
}

void SamplerEditor::resized()
{
  ScopedLock scopedLock(*lock);
  AudioModuleEditor::resized();
  int x = 0;
  int y = getPresetSectionBottom()+4;
  int w = getWidth();
  int h = getHeight() - y;

  // ...
}

/*
bugs: 


Maybe the xml presets for the SamplerModule should contain the filename for an sfz file that sould 
be loaded plus some global settings such as the maximum number of layers, the resampling algo, etc. 
Don't use the xml preset to store the actual sfz opcodes. Preset loading would become a nested 
process the outer level is the xml and the inner level is the sfz. The GUI should provide an 
additional set of load/save widgets for the sfz.
*/