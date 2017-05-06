
//-------------------------------------------------------------------------------------------------
// construction/destruction:

StraightlinerAudioModule::StraightlinerAudioModule(CriticalSection *newPlugInLock, rosic::Straightliner *straightlinerToWrap)
: PolyphonicInstrumentAudioModule(newPlugInLock, straightlinerToWrap)
{
  jassert(straightlinerToWrap != NULL); // you must pass a valid rosic-object to the constructor
  wrappedStraightliner      = straightlinerToWrap;
  underlyingRosicInstrument = straightlinerToWrap;
  setModuleName(juce::String("Straightliner"));
  setActiveDirectory(getApplicationDirectory() + juce::String("/StraightlinerPresets") );
  //oscSectionEditor->setActiveDirectory(pluginDir + juce::String(T("/StraightlinerPresets/OscSectionPresets")) );

  oscSectionModule = new FourOscSectionAudioModule(lock, &wrappedStraightliner->voiceArray[0].oscSection);
  oscSectionModule->setModuleName(juce::String("OscSection"));
  addChildAudioModule(oscSectionModule);

  filterModule = new MultiModeFilterAudioModule(lock, &wrappedStraightliner->voiceArray[0].filter);
  filterModule->setModuleName(juce::String("Filter"));
  addChildAudioModule(filterModule);

  pitchEnvModule = new BreakpointModulatorAudioModule(lock, &wrappedStraightliner->voiceArray[0].pitchEnv);
  pitchEnvModule->setModuleName(juce::String("PitchEnvelope"));
  addChildAudioModule(pitchEnvModule);

  filterEnvModule = new BreakpointModulatorAudioModule(lock, &wrappedStraightliner->voiceArray[0].filterEnv);
  filterEnvModule->setModuleName(juce::String("FilterEnvelope"));
  addChildAudioModule(filterEnvModule);

  ampEnvModule = new BreakpointModulatorAudioModule(lock, &wrappedStraightliner->voiceArray[0].ampEnv);
  ampEnvModule->setModuleName(juce::String("AmpEnvelope"));
  addChildAudioModule(ampEnvModule);
}

void StraightlinerAudioModule::setStateFromXml(const XmlElement &xmlState, const juce::String &stateName, bool markAsClean)
{
  // retrieve the patch format of the xml-file to enable different interpretations of the patch for 
  // backwards compatibility:
  int xmlPatchFormat = xmlState.getIntAttribute("PatchFormat", 0);

  // this override is specific to straightliner - in other plugins, we may hopefully rely on the StateManager 
  // infrastructure...
  if( xmlPatchFormat == 0 ) // this is an old preset
  {
    // in the old versions, we has the 4 oscillators each as child AudioModule, whereas in newer 
    // versions we have the whole oscillatro-section as one child module ....
    XmlElement* oscState = xmlState.getChildByName(T("Osc1"));
    if( oscState != NULL )
      oscSectionModule->osc1Module->setStateFromXml(*oscState, juce::String::empty, markAsClean);

    oscState = xmlState.getChildByName(T("Osc2"));
    if( oscState != NULL )
      oscSectionModule->osc2Module->setStateFromXml(*oscState, juce::String::empty, markAsClean);

    oscState = xmlState.getChildByName(T("Osc3"));
    if( oscState != NULL )
      oscSectionModule->osc3Module->setStateFromXml(*oscState, juce::String::empty, markAsClean);

    oscState = xmlState.getChildByName(T("Osc4"));
    if( oscState != NULL )
      oscSectionModule->osc4Module->setStateFromXml(*oscState, juce::String::empty, markAsClean);

    oscSectionModule->setStateName(juce::String::empty, true);

    // for restoring the other modules, we may invoke the baseclass' method:
    PolyphonicInstrumentAudioModule::setStateFromXml(xmlState, stateName, markAsClean);
  }
  else
    PolyphonicInstrumentAudioModule::setStateFromXml(xmlState, stateName, markAsClean);
}

bool StraightlinerAudioModule::checkForCrack()
{
  /*
  // when uncommented, this code will result in displaying Straightliner - cracked by ViP 
  // in the headline:
  juce::String directoryjuce::String = getApplicationDirectory();
  File   crackedKeyFile  = File(directoryjuce::String + File::separatorString + juce::String(T("Straightliner.xml")) );
  if( crackedKeyFile.existsAsFile() )
  {
    XmlDocument xmlDoc(crackedKeyFile);
    XmlElement* xmlKey = xmlDoc.getDocumentElement();
    if( xmlKey != NULL )
    {
      bool result = false;
      //if( !xmlKey->hasTagName(juce::String(T("Keyfile")))  )
      if( xmlKey->hasTagName(juce::String(T("Registration")))  )
      {
        moduleNameAppendix = juce::String(T(" - cracked by ViP"));
        result = true;
      }
      delete xmlKey;
      return result;
    }
  }
  \todo: keep strings that deal with copy-protection only in encoded form in memory (maybe not even temporarily) - decode them on-the-fly
  // while drawing
  */
  return false;
}