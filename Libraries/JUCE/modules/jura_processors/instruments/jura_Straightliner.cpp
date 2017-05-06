
//-------------------------------------------------------------------------------------------------
// construction/destruction:

StraightlinerAudioModule::StraightlinerAudioModule(CriticalSection *newPlugInLock, 
  rosic::Straightliner *straightlinerToWrap)
: PolyphonicInstrumentAudioModule(newPlugInLock, straightlinerToWrap)
{
  jassert(straightlinerToWrap != NULL); // you must pass a valid rosic-object to the constructor
  wrappedStraightliner      = straightlinerToWrap;
  underlyingRosicInstrument = straightlinerToWrap;
  setModuleName(juce::String("Straightliner"));
  setActiveDirectory(getApplicationDirectory() + juce::String("/StraightlinerPresets") );
  //oscSectionEditor->setActiveDirectory(pluginDir + juce::String(T("/StraightlinerPresets/OscSectionPresets")) );

  oscSectionModule = new FourOscSectionAudioModule(lock, 
    &wrappedStraightliner->voiceArray[0].oscSection);
  oscSectionModule->setModuleName(juce::String("OscSection"));
  addChildAudioModule(oscSectionModule);

  filterModule = new MultiModeFilterAudioModule(lock, 
    &wrappedStraightliner->voiceArray[0].filter);
  filterModule->setModuleName(juce::String("Filter"));
  addChildAudioModule(filterModule);

  pitchEnvModule = new BreakpointModulatorAudioModule(lock, 
    &wrappedStraightliner->voiceArray[0].pitchEnv);
  pitchEnvModule->setModuleName(juce::String("PitchEnvelope"));
  addChildAudioModule(pitchEnvModule);

  filterEnvModule = new BreakpointModulatorAudioModule(lock, 
    &wrappedStraightliner->voiceArray[0].filterEnv);
  filterEnvModule->setModuleName(juce::String("FilterEnvelope"));
  addChildAudioModule(filterEnvModule);

  ampEnvModule = new BreakpointModulatorAudioModule(lock, 
    &wrappedStraightliner->voiceArray[0].ampEnv);
  ampEnvModule->setModuleName(juce::String("AmpEnvelope"));
  addChildAudioModule(ampEnvModule);
}

void StraightlinerAudioModule::setStateFromXml(const XmlElement &xmlState, 
  const juce::String &stateName, bool markAsClean)
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
    XmlElement* oscState = xmlState.getChildByName("Osc1");
    if( oscState != NULL )
      oscSectionModule->osc1Module->setStateFromXml(*oscState, juce::String::empty, markAsClean);

    oscState = xmlState.getChildByName("Osc2");
    if( oscState != NULL )
      oscSectionModule->osc2Module->setStateFromXml(*oscState, juce::String::empty, markAsClean);

    oscState = xmlState.getChildByName("Osc3");
    if( oscState != NULL )
      oscSectionModule->osc3Module->setStateFromXml(*oscState, juce::String::empty, markAsClean);

    oscState = xmlState.getChildByName("Osc4");
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

//=================================================================================================

StraightlinerModuleEditor::StraightlinerModuleEditor(CriticalSection *newPlugInLock, 
  StraightlinerAudioModule* newStraightlinerAudioModule) 
  : PolyphonicInstrumentEditor(newPlugInLock, newStraightlinerAudioModule)
{
  setHeadlineStyle(MAIN_HEADLINE);

  jassert(newStraightlinerAudioModule != NULL ); // you must pass a valid module here
  straightlinerAudioModule = newStraightlinerAudioModule;
  stateWidgetSet->setLayout(StateLoadSaveWidgetSet::LABEL_AND_BUTTONS_ABOVE);

  oscSectionEditor = new FourOscSectionModuleEditor(lock, straightlinerAudioModule->oscSectionModule);
  addChildEditor( oscSectionEditor );
  oscSectionEditor->setDescriptionField(infoField, true);

  filterEditor = new MultiModeFilterModuleEditor(lock, straightlinerAudioModule->filterModule);
  addChildEditor( filterEditor );
  filterEditor->setDescriptionField(infoField, true);

  envelopeEditor = new BreakpointModulatorEditorMulti(lock, straightlinerAudioModule->pitchEnvModule);
  envelopeEditor->addModulatorToEdit(straightlinerAudioModule->filterEnvModule);
  envelopeEditor->addModulatorToEdit(straightlinerAudioModule->ampEnvModule);
  envelopeEditor->setModulatorLabel(0, juce::String("Pitch"));
  envelopeEditor->setModulatorLabel(1, juce::String("Filter"));
  envelopeEditor->setModulatorLabel(2, juce::String("Amplitude"));
  envelopeEditor->selectModulatorToEdit(2); // select amp envelope for editing
  addChildEditor( envelopeEditor );
  envelopeEditor->setHeadlineText(juce::String("Envelopes"));
  envelopeEditor->setDescription(juce::String("This is the editor for the 3 modulation generators"));
  envelopeEditor->addChangeListener(this);
  envelopeEditor->setDescriptionField(infoField, true);

  numHueOffsets = 2; // for osc- and filter-section

                     //loadColourScheme();
  updateSubEditorColourSchemes();
  updateWidgetsAccordingToState();
}

//-----------------------------------------------------------------------------------------------------------------------------------------
// callbacks:

void StraightlinerModuleEditor::copyColourSettingsFrom(const ColourSchemeComponent *componentToCopyFrom)
{
  AudioModuleEditor::copyColourSettingsFrom(componentToCopyFrom);
  oscSectionEditor->setCentralHue(editorColourScheme.getCentralHue() + editorColourScheme.getHueOffset(0));
  filterEditor->setCentralHue(    editorColourScheme.getCentralHue() + editorColourScheme.getHueOffset(1));
}

void StraightlinerModuleEditor::updateWidgetsAccordingToState()
{
  if( straightlinerAudioModule == NULL )
    return;
  if( straightlinerAudioModule->wrappedStraightliner == NULL )
    return;

  // remember if the preset was clean or dirty before making a few calls that may lead to a 
  // dirtification of the preset-state:
  bool presetIsDirty = straightlinerAudioModule->isStateDirty();
  const MessageManagerLock mmLock;     
  // the event loop will now be locked so it's safe to make a few calls..

  // update the global widgets and automatable sliders:
  PolyphonicInstrumentEditor::updateWidgetsAccordingToState();

  // update the sub-editors:
  oscSectionEditor->updateWidgetsAccordingToState();
  filterEditor->updateWidgetsAccordingToState();
  envelopeEditor->updateWidgetsAccordingToState();

  // preserve the clean/dirty state of the preset regardless of any parameter changes that may take
  // place that may take place - note that not sending a change notification from the widgets is 
  // not enough to make that sure because some of the have AutomatableParameters associated with 
  // them which themselves may dirtify the preset:
  if( presetIsDirty )
    straightlinerAudioModule->markStateAsDirty();
  else
    straightlinerAudioModule->markStateAsClean();
}

void StraightlinerModuleEditor::paint(Graphics &g)
{
  if( drawGradientsBasedOnOutlines == false ) 
  {
    fillRectWithBilinearGradient(g, 0, 0, getWidth(), getHeight(), editorColourScheme.topLeft, editorColourScheme.topRight, 
      editorColourScheme.bottomLeft, editorColourScheme.bottomRight);  
    Editor::drawHeadline(g);
    return;
  }

  int x1, x2, y1, y2, w, w2, w3, h, h2, h3;

  Colour gradientColourOsc = oscSectionEditor->getTopLeftColour();
  Colour gradientColourFlt = filterEditor->getTopLeftColour();
  Colour gradientColourEnv = envelopeEditor->getTopLeftColour();

  // gradient line between osc-section and filter:
  x1 = oscSectionEditor->getRight();
  x2 = filterEditor->getX();
  y1 = oscSectionEditor->getY();
  y2 = oscSectionEditor->getBottom();
  w  = x2-x1;
  h  = y2-y1;
  h2 = h/2;
  h3 = h-h2;

  fillRectWithBilinearGradient(g, x1, y1+h2, w, h3, 
    gradientMidColour, gradientMidColour, 
    gradientColourOsc, gradientColourFlt );

  Colour mixColour = getMixedColour(gradientColourOsc, gradientColourFlt);

  fillRectWithBilinearGradient(g, x1, y1, w, h2, 
    mixColour,         mixColour, 
    gradientMidColour, gradientMidColour );

  // the square in the center:
  y1 = oscSectionEditor->getBottom();
  y2 = envelopeEditor->getY();
  h  = y2-y1;
  fillRectWithBilinearGradient(g, x1, y1, w, h, 
    gradientColourOsc, gradientColourFlt, 
    gradientColourEnv, gradientColourEnv );

  // the line between amp-env and filter-env:
  y1 = envelopeEditor->getY();
  y2 = envelopeEditor->getBottom();
  h  = y2-y1;
  h2 = h/2;
  h3 = h-h2;
  fillRectWithBilinearGradient(g, x1, y1, w, h2, 
    gradientColourEnv, gradientColourEnv, 
    gradientMidColour, gradientMidColour );

  mixColour = getMixedColour(gradientColourEnv, gradientColourEnv);

  fillRectWithBilinearGradient(g, x1, y1+h2, w, h3, 
    gradientMidColour, gradientMidColour,
    mixColour,         mixColour);

  // the line between osc-section and amp-env:
  x1 = oscSectionEditor->getX();
  x2 = oscSectionEditor->getRight();
  y1 = oscSectionEditor->getBottom();
  y2 = envelopeEditor->getY();
  w  = x2-x1;
  w2 = w/2;
  w3 = w-w2;
  h  = y2-y1;
  fillRectWithBilinearGradient(g, x1+w2, y1, w3, h, gradientMidColour, gradientColourOsc, gradientMidColour, gradientColourEnv );

  mixColour = getMixedColour(gradientColourOsc, gradientColourEnv);

  fillRectWithBilinearGradient(g, x1, y1, w2, h, mixColour, gradientMidColour, mixColour, gradientMidColour );

  // the line between filter and filter-env:
  x1 = filterEditor->getX();
  x2 = filterEditor->getRight();
  w  = x2-x1;
  w2 = w/2;
  w3 = w-w2;
  fillRectWithBilinearGradient(g, x1, y1, w2, h, gradientColourFlt, gradientMidColour, gradientColourEnv, gradientMidColour);

  mixColour = getMixedColour(gradientColourFlt, gradientColourEnv);

  fillRectWithBilinearGradient(g, x1+w2, y1, w3, h, gradientMidColour, mixColour, gradientMidColour, mixColour );

  // the line left to osc-section and amp-env:
  x1 = 0;
  x2 = oscSectionEditor->getX();
  y1 = oscSectionEditor->getY();
  y2 = envelopeEditor->getBottom();
  w  = x2-x1;
  h  = y2-y1;
  fillRectWithBilinearGradient(g, x1, y1, w, h, gradientColourOsc, gradientColourOsc, gradientColourEnv, gradientColourEnv  );

  // the line right to filter and filter env:
  x1 = filterEditor->getRight();
  x2 = getWidth();
  w  = x2-x1;
  fillRectWithBilinearGradient(g, x1, y1, w, h, gradientColourFlt, gradientColourFlt, gradientColourEnv, gradientColourEnv  );

  // the global section at the top:
  x1 = 0;
  x2 = getRight();
  y2 = oscSectionEditor->getY();
  y1 = 0;
  w  = x2-x1;
  h  = y2-y1;
  fillRectWithBilinearGradient(g, x1, y1, w, h, editorColourScheme.topLeft, editorColourScheme.topRight, 
    gradientColourOsc, gradientColourFlt);  

  // the info-line:
  y1 = envelopeEditor->getBottom();
  y2 = getBottom();
  h  = y2-y1;
  fillRectWithBilinearGradient(g, x1, y1, w, h, gradientColourEnv, gradientColourEnv, 
    editorColourScheme.bottomLeft, editorColourScheme.bottomRight);  

  Editor::drawHeadline(g);
}

void StraightlinerModuleEditor::resized()
{
  AudioModuleEditor::resized();
  int x = 0;
  int y = 0;
  int w = getWidth();
  int h = getHeight();

  x = 0;
  y = getHeadlineBottom();
  w = getWidth()/5;

  stateWidgetSet->setBounds(x+4, y+4, w-8, 40);

  x = w;
  y = getHeadlineBottom()+4;
  levelSlider->setBounds(       x+4,     y, w-8,   20);
  y += 18;
  levelByKeySlider->setBounds(  x+4,     y, w/2-8, 16);
  levelByVelSlider->setBounds(  x+w/2+4, y, w/2-8, 16);

  x = 2*w;
  y = getHeadlineBottom()+4;
  midSideRatioSlider->setBounds(x+4, y, w-8, 14);
  y += 16;
  numVoicesSlider->setBounds(   x+4, y, w-8, 14);
  y += 16;
  compSlider->setBounds(        x+4, y, w-8, 14);

  x = 3*w;
  y = getHeadlineBottom();
  glideButton->setBounds(x+4, y+4, 40, 20);
  glideTimeSlider->setBounds(glideButton->getRight()+4, y+4, w-glideButton->getWidth()-12, 20);
  y += 24;
  masterTuneSlider->setBounds(x+4, y+4, w/2-8, 16);
  x += w/2;
  wheelRangeSlider->setBounds(x+4, y+4, w/2-8, 16);

  x = 4*w;
  y = getHeadlineBottom();

  tuningLabel->setBounds(x+4, y+4, 92, 20);
  tuningPlusButton->setBounds(x+w-4-20, y+4, 20, 20);
  tuningMinusButton->setBounds(tuningPlusButton->getX()-18, y+4, 20, 20);
  tuningLoadButton->setBounds(tuningMinusButton->getX()-40+2, y+4, 40, 20);

  y = tuningLabel->getBottom()-2;
  tuningFileNameLabel->setBounds(tuningLabel->getX(), y, w-8, 20);

  x = 0;
  y = stateWidgetSet->getBottom()+6;
  w = getWidth()/2;
  h = 304;

  oscSectionEditor->setBounds(x+4, y+8, w-8, h);

  y = oscSectionEditor->getY();
  x = oscSectionEditor->getRight()+8;
  h = oscSectionEditor->getBottom()-y;
  w = oscSectionEditor->getWidth();

  filterEditor->setBounds(x, y, w, h);

  x = oscSectionEditor->getX();
  y = oscSectionEditor->getBottom()+8;
  w = oscSectionEditor->getWidth();
  h = 256;

  x = oscSectionEditor->getX();
  y = oscSectionEditor->getBottom()+8;
  w = filterEditor->getRight()-x;
  //h = 256;
  h = 260;
  envelopeEditor->setBounds(x, y, w, h);
}

//-----------------------------------------------------------------------------------------------------------------------------------------
// others:

/*
void StraightlinerModuleEditor::loadPreferencesFromFile()
{
AudioModuleEditor::loadPreferencesFromFile();
envelopeEditor->copyColourSettingsFrom(this);
envelopeEditor->repaint();
}
*/

/*
void StraightlinerModuleEditor::setColourSchemeFromXml(const XmlElement& xmlColorScheme)
{
PolyphonicInstrumentEditor::setColourSchemeFromXml(xmlColorScheme);
envelopeEditor->copyColourSettingsFrom(this);
}
*/

void StraightlinerModuleEditor::updateSubEditorColourSchemes()
{
  envelopeEditor->copyColourSettingsFrom(this);
  oscSectionEditor->copyColourSettingsFrom(this);
  filterEditor->copyColourSettingsFrom(this);

  oscSectionEditor->setCentralHue(editorColourScheme.getCentralHue() + editorColourScheme.getHueOffset(0));
  filterEditor->setCentralHue(    editorColourScheme.getCentralHue() + editorColourScheme.getHueOffset(1));
}