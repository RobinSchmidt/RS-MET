
VectorSamplePlayerAudioModule::VectorSamplePlayerAudioModule(CriticalSection *newPlugInLock,      
  rosic::VectorSamplePlayer *vectorSamplePlayerToWrap)  : AudioModule(newPlugInLock)
{
  jassert(vectorSamplePlayerToWrap != NULL); // you must pass a valid rosic-object to the constructor
  wrappedVectorSamplePlayer = vectorSamplePlayerToWrap;
  moduleName = juce::String("VectorSamplePlayer");

  // initialize the current directory for preset loading and saving:
  setActiveDirectory(getApplicationDirectory() + juce::String("/VectorSamplePlayerPresets"));

  vectorMixerModule = new VectorMixerAudioModule(lock, &wrappedVectorSamplePlayer->vectorMixer);
  vectorMixerModule->setModuleName(juce::String("VectorMixer"));
  addChildAudioModule(vectorMixerModule);

  samplePlayerTopLeftModule = new SamplePlayerAudioModule(lock, 
    &wrappedVectorSamplePlayer->samplePlayerTopLeft);
  samplePlayerTopLeftModule->setModuleName(juce::String("SamplePlayerTopLeft"));
  addChildAudioModule(samplePlayerTopLeftModule);

  samplePlayerTopRightModule = new SamplePlayerAudioModule(lock, 
    &wrappedVectorSamplePlayer->samplePlayerTopRight);
  samplePlayerTopRightModule->setModuleName(juce::String("SamplePlayerTopRight"));
  addChildAudioModule(samplePlayerTopRightModule);

  samplePlayerBottomLeftModule = new SamplePlayerAudioModule(lock, 
    &wrappedVectorSamplePlayer->samplePlayerBottomLeft);
  samplePlayerBottomLeftModule->setModuleName(juce::String("SamplePlayerBottomLeft"));
  addChildAudioModule(samplePlayerBottomLeftModule);

  samplePlayerBottomRightModule = new SamplePlayerAudioModule(lock, 
    &wrappedVectorSamplePlayer->samplePlayerBottomRight);
  samplePlayerBottomRightModule->setModuleName(juce::String("SamplePlayerBottomRight"));
  addChildAudioModule(samplePlayerBottomRightModule);

  xLfoModule = new LowFrequencyOscillatorAudioModule(lock, &wrappedVectorSamplePlayer->xLfo);
  xLfoModule->setModuleName(juce::String("X-LFO"));
  addChildAudioModule(xLfoModule);

  yLfoModule = new LowFrequencyOscillatorAudioModule(lock, &wrappedVectorSamplePlayer->yLfo);
  yLfoModule->setModuleName(juce::String("Y-LFO"));
  addChildAudioModule(yLfoModule);
}

//=================================================================================================

VectorSamplePlayerSampleEditor::VectorSamplePlayerSampleEditor(CriticalSection *newPlugInLock,                                                              
  SamplePlayerAudioModule* newSamplePlayerAudioModule)
  : SamplePlayerModuleEditor(newPlugInLock, newSamplePlayerAudioModule)
  , AudioModuleEditor( newSamplePlayerAudioModule)
{ 
  setLinkPosition(INVISIBLE);
  setPresetSectionPosition(INVISIBLE);
  sampleDisplayZoomer->setVisible(false);
  layout        = 0;
  headlineWidth = 48;
}

void VectorSamplePlayerSampleEditor::setLayout(int newLayout) 
{ 
  layout = newLayout; 
  switch( layout )
  {
  case TOP_LEFT:     Editor::setHeadlinePosition(Editor::BOTTOM_LEFT);  break;
  case TOP_RIGHT:    Editor::setHeadlinePosition(Editor::BOTTOM_RIGHT); break;
  case BOTTOM_LEFT:  Editor::setHeadlinePosition(Editor::TOP_LEFT);     break;
  case BOTTOM_RIGHT: Editor::setHeadlinePosition(Editor::TOP_RIGHT);    break;
  }
  resized(); 
}

void VectorSamplePlayerSampleEditor::resized()
{
  Editor::resized();
  int x = 0;
  int y = 0;
  int w = getWidth();
  int h = getHeight();

  w = 150;
  h = 16;
  if( layout == TOP_LEFT )
  {
    x = 0;
    y = 0;
    w = getWidth()-horizontalIntrusion;
    h = getHeight()-y-verticalIntrusion;
    sampleDisplay->setBounds(x, y, w, h);

    y = sampleDisplay->getBottom();
    int labelWidth = w-72;
    x = w-labelWidth;
    sampleFileNameLabel->setBounds(x, y, labelWidth, 16);
    x = sampleFileNameLabel->getX()-40-26;
    sampleLoadButton->setBounds(x, y, 40, 16);
    sampleMinusButton->setBounds(sampleLoadButton->getRight()-2,  y, 16, 16);
    samplePlusButton->setBounds( sampleMinusButton->getRight()-2, y, 16, 16);

    x = getWidth()-horizontalIntrusion-2;
    y = 4;
    muteButton->setBounds(          x   +4, y   +4, 40, 16);
    soloButton->setBounds(          x+38+4, y   +4, 40, 16);
    loopButton->setBounds(          x   +4, y+14+4, 40, 16);
    phaseRandomizeButton->setBounds(x   +4, y+28+4, 40, 16);
    y = loopButton->getY()+8;
    x = loopButton->getRight()+4;
    moreButton->setBounds(x, y, 40, 16);

    y = getHeight()-16;
    x = headlineWidth;
    w = (getWidth()-horizontalIntrusion-x)/2;
    tuneSlider->setBounds( x+4,   y-4, w-8, 16);
    levelSlider->setBounds(x+w+4, y-4, w-8, 16);
  }
  else if( layout == TOP_RIGHT )
  {
    x = horizontalIntrusion;
    y = 0;
    w = getWidth()-x;
    h = getHeight()-y-verticalIntrusion;
    sampleDisplay->setBounds(x, y, w, h);

    y = sampleDisplay->getBottom();
    sampleFileNameLabel->setBounds(x, y, w-72, 16);
    x = sampleFileNameLabel->getRight();
    sampleLoadButton->setBounds(x-2, y, 40, 16);
    sampleMinusButton->setBounds(sampleLoadButton->getRight()-2,  y, 16, 16);
    samplePlusButton->setBounds( sampleMinusButton->getRight()-2, y, 16, 16);

    x = horizontalIntrusion-40-6;
    y = 4;
    muteButton->setBounds(          x   +4, y   +4, 40, 16);
    soloButton->setBounds(          x-38+4, y   +4, 40, 16);
    loopButton->setBounds(          x   +4, y+14+4, 40, 16);
    phaseRandomizeButton->setBounds(x   +4, y+28+4, 40, 16);
    y = loopButton->getY()+8;
    x = loopButton->getX()-40-4;
    moreButton->setBounds(x, y, 40, 16);

    x = horizontalIntrusion;
    y = getHeight()-16;
    w = (headlineX-x)/2;
    levelSlider->setBounds( x+4,   y-4, w-8, 16);
    tuneSlider->setBounds(x+w+4, y-4, w-8, 16);
  }
  else if( layout == BOTTOM_LEFT )
  {
    x = 0;
    y = verticalIntrusion;
    w = getWidth()-horizontalIntrusion;
    h = getHeight()-y;
    sampleDisplay->setBounds(x, y, w, h);

    int labelWidth = w-72;
    x  = w-labelWidth;
    y -= 16;
    sampleFileNameLabel->setBounds(x, y, labelWidth, 16);
    x = sampleFileNameLabel->getX()-40-26;
    sampleLoadButton->setBounds(x, y, 40, 16);
    sampleMinusButton->setBounds(sampleLoadButton->getRight()-2,  y, 16, 16);
    samplePlusButton->setBounds( sampleMinusButton->getRight()-2, y, 16, 16);

    x = getWidth()-horizontalIntrusion-2;
    y = getHeight()-10-16;
    muteButton->setBounds(          x   +4, y   +4, 40, 16);
    soloButton->setBounds(          x+38+4, y   +4, 40, 16);
    loopButton->setBounds(          x   +4, y-14+4, 40, 16);
    phaseRandomizeButton->setBounds(x   +4, y-28+4, 40, 16);
    y = loopButton->getY()-8;
    x = loopButton->getRight()+4;
    moreButton->setBounds(x, y, 40, 16);

    x = headlineWidth;
    y = 0;
    w = (getWidth()-horizontalIntrusion-x)/2;
    tuneSlider->setBounds( x+4,   y+4, w-8, 16);
    levelSlider->setBounds(x+w+4, y+4, w-8, 16);
  }
  else if( layout == BOTTOM_RIGHT )
  {
    x = horizontalIntrusion;
    y = verticalIntrusion;
    w = getWidth()-x;
    h = getHeight()-y;
    sampleDisplay->setBounds(x, y, w, h);
    y -= 16;
    sampleFileNameLabel->setBounds(x, y, w-72, 16);
    x = sampleFileNameLabel->getRight();
    sampleLoadButton->setBounds(x-2, y, 40, 16);
    sampleMinusButton->setBounds(sampleLoadButton->getRight()-2,  y, 16, 16);
    samplePlusButton->setBounds( sampleMinusButton->getRight()-2, y, 16, 16);

    x = horizontalIntrusion-40-6;
    y = getHeight()-10-16;
    muteButton->setBounds(          x   +4, y   +4, 40, 16);
    soloButton->setBounds(          x-38+4, y   +4, 40, 16);
    loopButton->setBounds(          x   +4, y-14+4, 40, 16);
    phaseRandomizeButton->setBounds(x   +4, y-28+4, 40, 16);
    y = loopButton->getY()-8;
    x = loopButton->getX()-40-4;
    moreButton->setBounds(x, y, 40, 16);

    x = horizontalIntrusion;
    y = 0;
    w = (headlineX-x)/2;
    levelSlider->setBounds( x+4,   y+4, w-8, 16);
    tuneSlider->setBounds(x+w+4, y+4, w-8, 16);
  }
}

//=================================================================================================
// class VectorSamplePlayerLfoEditor:

VectorSamplePlayerLfoEditor::VectorSamplePlayerLfoEditor(CriticalSection *newPlugInLock, 
  LowFrequencyOscillatorAudioModule* newLowFrequencyOscillatorAudioModule)
  : LowFrequencyOscillatorEditor(newPlugInLock, newLowFrequencyOscillatorAudioModule)
{ 
  setLinkPosition(INVISIBLE);
  setPresetSectionPosition(INVISIBLE);
  layout        = 0;
  headlineWidth = 60;

  cycleLengthSlider->setStringConversionFunction(valueToString3);
}

void VectorSamplePlayerLfoEditor::setLayout(int newLayout) 
{ 
  layout = newLayout; 
}

void VectorSamplePlayerLfoEditor::paint(Graphics &g)
{
  LowFrequencyOscillatorEditor::paint(g);

  /*
  float x1, x2, x3, y1, y2, y3;
  if( layout == LEFT )
  {
  x1 = (float) waveformDisplay->getRight();
  y1 = (float) waveformDisplay->getY();
  x2 = (float) waveformDisplay->getRight();
  y2 = (float) waveformDisplay->getBottom();
  x3 = (float) getWidth();
  y3 = (float) (y1+y2)/2;
  }
  else
  {
  x1 = (float) waveformDisplay->getX();
  y1 = (float) waveformDisplay->getY();
  x2 = (float) waveformDisplay->getX();
  y2 = (float) waveformDisplay->getBottom();
  x3 = (float) 0;
  y3 = (float) (y1+y2)/2;
  }  

  g.setColour(Colours::white); 
  rojue::drawTriangle(g, x1, y1, x2, y2, x3, y3, true);
  g.setColour(Colours::lavender); 
  rojue::drawTriangle(g, x1, y1, x2, y2, x3, y3, false);
  */

  drawHeadline(g);
}

void VectorSamplePlayerLfoEditor::resized()
{
  LowFrequencyOscillatorEditor::resized();

  /*
  int x = 0;
  int y = 0;
  int w = getWidth();
  int h = getHeight();

  triggerButton->setVisible(false);  // we don't use it here

  if( layout == RIGHT )
  {
  headlineX = 4;
  headlineY = getHeight()/2 - boldFont16px.getFontHeight()/2;

  x = headlineWidth;
  y = 24;
  w = getWidth()-x;
  h = getHeight()-48;
  waveformDisplay->setBounds(x, y, w, h);
  emptyDisplay->setBounds(   x, y, w, h);

  x = 8;
  y = waveformDisplay->getBottom();
  w = getWidth()-x-72;
  waveformComboBox->setBounds(x, y, w, 16);
  x = waveformComboBox->getRight()-2;
  waveLoadButton->setBounds(x, y, 40, 16);
  waveMinusButton->setBounds(waveLoadButton->getRight()-2, y, 16, 16);
  wavePlusButton->setBounds(waveMinusButton->getRight()-2, y, 16, 16);
  x = waveformComboBox->getX();
  y = waveformComboBox->getY()-14;
  onOffButton->setBounds(x, y, 32, 16);

  x = getWidth()-44;
  y = 4;
  editButton->setBounds(x, y, 40, 16);

  w = editButton->getX();
  cycleLengthSlider->setBounds(4,     y, w/2-4, 16);
  depthSlider->setBounds(      w/2+4, y, w/2-8, 16);

  x = cycleLengthSlider->getX();
  y = cycleLengthSlider->getY()+14;
  tempoSyncButton->setBounds(x, y, 40, 16);
  }
  else
  {
  w = boldFont16px.getTextPixelWidth(getHeadlineString(), boldFont16px.getDefaultKerning());
  headlineX = getWidth()-w-4;
  headlineY = getHeight()/2 - boldFont16px.getFontHeight()/2;

  x = 0;
  y = 24;
  w = getWidth()-headlineWidth;
  h = getHeight()-48;
  waveformDisplay->setBounds(x, y, w, h);
  emptyDisplay->setBounds(   x, y, w, h);

  x = 6;
  y = waveformDisplay->getBottom();
  waveLoadButton->setBounds(x, y, 40, 16);
  waveMinusButton->setBounds(waveLoadButton->getRight()-2, y, 16, 16);
  wavePlusButton->setBounds(waveMinusButton->getRight()-2, y, 16, 16);
  w = getWidth()-x-74; // why not -72 as above
  waveformComboBox->setBounds(wavePlusButton->getRight()-2, y, w, 16);
  x = waveformComboBox->getRight()-32;
  y = waveformComboBox->getY()-14;
  onOffButton->setBounds(x, y, 32, 16);

  x = 4;
  y = 4;
  editButton->setBounds(x, y, 40, 16);

  x = editButton->getRight();
  w = getWidth() - x;
  depthSlider->setBounds(       x+4,     y, w/2-4, 16);
  cycleLengthSlider->setBounds( x+w/2+4, y, w/2-8, 16);

  x = cycleLengthSlider->getRight()-40;
  y = cycleLengthSlider->getY()+14;
  tempoSyncButton->setBounds(x, y, 40, 16);
  }
  */
}

//=================================================================================================
// class VectorSamplePlayerModuleEditor:

VectorSamplePlayerEditor::VectorSamplePlayerEditor(CriticalSection *newPlugInLock, 
  VectorSamplePlayerAudioModule* newVectorSamplePlayerAudioModule) 
  : AudioModuleEditor(newVectorSamplePlayerAudioModule)
{
  //setHeadlineStyle(MAIN_HEADLINE);
  setHeadlineText( juce::String("VectorSamplePlayer") );

  //setHeadlinePosition(AudioModuleEditor::TOP_LEFT);
  //setPresetSectionPosition(AudioModuleEditor::BELOW_HEADLINE);
  //stateWidgetSet->setLayout(StateLoadSaveWidgetSet::LABEL_AND_BUTTONS_ABOVE);

  // assign the pointer to the rosic::VectorSamplePlayer object to be used as aduio engine:
  jassert(newVectorSamplePlayerAudioModule != NULL ); // you must pass a valid module here
  vectorSamplePlayerAudioModule = newVectorSamplePlayerAudioModule;

  //---------------------------------------------------------------------------
  // create and setup the sub-module editors:

  samplePlayerTopLeftEditor = new VectorSamplePlayerSampleEditor(lock, 
    vectorSamplePlayerAudioModule->samplePlayerTopLeftModule);
  samplePlayerTopLeftEditor->setLayout(VectorSamplePlayerSampleEditor::TOP_LEFT);
  samplePlayerTopLeftEditor->setHeadlineText(juce::String("Osc1"));
  addChildEditor( samplePlayerTopLeftEditor );
  samplePlayerTopLeftEditor->addChangeListener(this);
  samplePlayerTopLeftEditor->setDescriptionField(infoField, true);

  samplePlayerTopRightEditor = new VectorSamplePlayerSampleEditor(lock, 
    vectorSamplePlayerAudioModule->samplePlayerTopRightModule);
  samplePlayerTopRightEditor->setLayout(VectorSamplePlayerSampleEditor::TOP_RIGHT);
  samplePlayerTopRightEditor->setHeadlineText(juce::String("Osc2"));
  addChildEditor( samplePlayerTopRightEditor );
  samplePlayerTopRightEditor->addChangeListener(this);
  samplePlayerTopRightEditor->setDescriptionField(infoField, true);

  samplePlayerBottomLeftEditor = new VectorSamplePlayerSampleEditor(lock, 
    vectorSamplePlayerAudioModule->samplePlayerBottomLeftModule);
  samplePlayerBottomLeftEditor->setLayout(VectorSamplePlayerSampleEditor::BOTTOM_LEFT);
  samplePlayerBottomLeftEditor->setHeadlineText(juce::String("Osc3"));
  addChildEditor( samplePlayerBottomLeftEditor );
  samplePlayerBottomLeftEditor->addChangeListener(this);
  samplePlayerBottomLeftEditor->setDescriptionField(infoField, true);

  samplePlayerBottomRightEditor = new VectorSamplePlayerSampleEditor(lock, 
    vectorSamplePlayerAudioModule->samplePlayerBottomRightModule);
  samplePlayerBottomRightEditor->setLayout(VectorSamplePlayerSampleEditor::BOTTOM_RIGHT);
  samplePlayerBottomRightEditor->setHeadlineText(juce::String("Osc4"));
  addChildEditor( samplePlayerBottomRightEditor );
  samplePlayerBottomRightEditor->addChangeListener(this);
  samplePlayerBottomRightEditor->setDescriptionField(infoField, true);

  addAndMakeVisible( vectorMixerPad = new VectorMixerModuleEditor(lock, 
    vectorSamplePlayerAudioModule->vectorMixerModule) );

  xLfoEditor = new VectorSamplePlayerLfoEditor(lock, vectorSamplePlayerAudioModule->xLfoModule);
  xLfoEditor->setHeadlineStyle(Editor::SUB_HEADLINE);
  xLfoEditor->setDescription(juce::String("Low frequency oscillator (LFO) for the x-coordinate"));
  addChildEditor( xLfoEditor );

  yLfoEditor = new VectorSamplePlayerLfoEditor(lock, vectorSamplePlayerAudioModule->yLfoModule);
  yLfoEditor->setDescription(juce::String("Low frequency oscillator (LFO) for the y-coordinate"));
  yLfoEditor->setHeadlineStyle(Editor::SUB_HEADLINE);
  yLfoEditor->setLayout(VectorSamplePlayerLfoEditor::RIGHT);
  addChildEditor( yLfoEditor );

  // set up the widgets:
  updateWidgetsAccordingToState();
}

//-------------------------------------------------------------------------------------------------
// callbacks:

void VectorSamplePlayerEditor::updateWidgetsAccordingToState()
{
  if( vectorSamplePlayerAudioModule == NULL )
    return;
  if( vectorSamplePlayerAudioModule->wrappedVectorSamplePlayer == NULL )
    return;

  // remember if the preset was clean or dirty before making a few calls that may lead to a 
  // dirtification of the preset-state:
  bool presetIsDirty = vectorSamplePlayerAudioModule->isStateDirty();
  const MessageManagerLock mmLock;     
  // the event loop will now be locked so it's safe to make a few calls..

  // update the global widgets and automatable sliders:
  AudioModuleEditor::updateWidgetsAccordingToState();

  // update the sub-editors:
  samplePlayerTopLeftEditor->updateWidgetsAccordingToState();
  samplePlayerTopRightEditor->updateWidgetsAccordingToState();
  samplePlayerBottomLeftEditor->updateWidgetsAccordingToState();
  samplePlayerBottomRightEditor->updateWidgetsAccordingToState();
  vectorMixerPad->updateWidgetsAccordingToState();
  xLfoEditor->updateWidgetsAccordingToState();
  yLfoEditor->updateWidgetsAccordingToState();

  // preserve the clean/dirty state of the preset regardless of any parameter changes that may take
  // place that may take place - note that not sending a change notification from the widgets is 
  // not enough to make that sure because some of the have AutomatableParameters associated with 
  // them which themselves may dirtify the preset:
  if( presetIsDirty )  
    vectorSamplePlayerAudioModule->markStateAsDirty();
  else
    vectorSamplePlayerAudioModule->markStateAsClean();
}

void VectorSamplePlayerEditor::resized()
{
  AudioModuleEditor::resized();
  int x = 0;
  int y = getPresetSectionBottom()+4;
  int w = getWidth()/2;
  int h = getHeight();

  int vi = VectorSamplePlayerSampleEditor::verticalIntrusion;    // 40
  int hi = VectorSamplePlayerSampleEditor::horizontalIntrusion;  // 90
  int padSize = 2*hi;
  int sampleEditorHeight = 100;

  h  = sampleEditorHeight;
  samplePlayerTopLeftEditor->setBounds(    x,   y,   w, h);
  samplePlayerTopRightEditor->setBounds(   x+w, y,   w, h);
  y += padSize - 2*vi;
  samplePlayerBottomLeftEditor->setBounds( x,   y+h, w, h);
  samplePlayerBottomRightEditor->setBounds(x+w, y+h, w, h);

  w = h = padSize;
  x = samplePlayerTopLeftEditor->getRight()  - hi;
  y = samplePlayerTopLeftEditor->getBottom() - vi;
  vectorMixerPad->setBounds(x, y, w, h);

  x = 0;
  y = samplePlayerTopLeftEditor->getBottom();
  w = vectorMixerPad->getX()-x;
  h = samplePlayerBottomLeftEditor->getY()-y;
  xLfoEditor->setBounds(x, y, w, h);
  x += w + vectorMixerPad->getWidth();
  yLfoEditor->setBounds(x, y, w, h);
}