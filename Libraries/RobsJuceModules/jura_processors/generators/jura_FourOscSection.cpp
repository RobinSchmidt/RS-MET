
//-------------------------------------------------------------------------------------------------
// construction/destruction:

FourOscSectionAudioModule::FourOscSectionAudioModule(CriticalSection *newPlugInLock, rosic::FourOscSection *fourOscSectionToWrap)
 : AudioModule(newPlugInLock)
{
  jassert(fourOscSectionToWrap != NULL); // you must pass a valid rosic-object to the constructor
  wrappedFourOscSection = fourOscSectionToWrap;

  setModuleTypeName("FourOscillatorSection");

  osc1Module = new WaveOscModule(lock, &wrappedFourOscSection->osc1);
  osc1Module->setModuleName("Osc1");
  addChildAudioModule(osc1Module);

  osc2Module = new WaveOscModule(lock, &wrappedFourOscSection->osc2);
  osc2Module->setModuleName("Osc2");
  addChildAudioModule(osc2Module);

  osc3Module = new WaveOscModule(lock, &wrappedFourOscSection->osc3);
  osc3Module->setModuleName("Osc3");
  addChildAudioModule(osc3Module);

  osc4Module = new WaveOscModule(lock, &wrappedFourOscSection->osc4);
  osc4Module->setModuleName("Osc4");
  addChildAudioModule(osc4Module);
}

//=================================================================================================

FourOscSectionModuleEditor::FourOscSectionModuleEditor(CriticalSection *newPlugInLock, 
  FourOscSectionAudioModule *newFourOscSectionToEdit) : AudioModuleEditor(newFourOscSectionToEdit) 
{
  jassert(moduleToEdit != NULL); // you must pass a valid object here
  stateWidgetSet->stateLabel->setVisible(false);
  webLink->setVisible(false);
  infoField->setVisible(false);

  setHeadlineText("Oscillators");
  setDescription("This is the oscillator section");

  addChildEditor( osc1Editor = new WaveOscEditor(lock, newFourOscSectionToEdit->osc1Module));
  addChildEditor( osc2Editor = new WaveOscEditor(lock, newFourOscSectionToEdit->osc2Module) );
  addChildEditor( osc3Editor = new WaveOscEditor(lock, newFourOscSectionToEdit->osc3Module) );
  addChildEditor( osc4Editor = new WaveOscEditor(lock, newFourOscSectionToEdit->osc4Module) );

  stateWidgetSet->stateLoadButton->setDescription(   "Load oscillator section setting from file");
  stateWidgetSet->stateSaveButton->setDescription(   "Save oscillator section setting to file");
  stateWidgetSet->statePlusButton->setDescription(   "Skip to next oscillator section setting in current directory");
  stateWidgetSet->stateMinusButton->setDescription(  "Skip to previous oscillator section setting in current directory");
  stateWidgetSet->stateFileNameLabel->setDescription("Name of current preset for the oscillator section (if any)");

  isTopLevelEditor = false;

  //setSize(272, 136);  // widget arrangement is optimized for this size
}

void FourOscSectionModuleEditor::resized()
{
  AudioModuleEditor::resized();

  int y  = getPresetSectionBottom()+4;
  int w  = getWidth();
  int h  = getHeight()-y;
  int w2 = w/2+1;
  int h2 = h/2+1;

  osc1Editor->setBounds( 0,        y, w2, h2);
  osc2Editor->setBounds(w2-2,      y, w2, h2);
  osc3Editor->setBounds( 0,   y+h2-2, w2, h2);
  osc4Editor->setBounds(w2-2, y+h2-2, w2, h2);
}

void FourOscSectionModuleEditor::updateWidgetsAccordingToState()
{
  osc1Editor->updateWidgetsAccordingToState();
  osc2Editor->updateWidgetsAccordingToState();
  osc3Editor->updateWidgetsAccordingToState();
  osc4Editor->updateWidgetsAccordingToState();
  stateWidgetSet->updateStateNameField();
}
