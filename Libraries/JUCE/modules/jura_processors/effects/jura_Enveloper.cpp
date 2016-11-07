
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
  return new jura::BreakpointModulatorEditor(plugInLock, &envGenWrapper);
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

//void Enveloper::reset()
//{
// 
//}

// delete, if it's clear that we won't need it:
////=================================================================================================
//// the GUI editor class for the Ladder:
//
//LadderEditor::LadderEditor(jura::Ladder *newLadderToEdit) : AudioModuleEditor(newLadderToEdit)
//{
//  ScopedLock scopedLock(*plugInLock);
//  // maybe we should avoid this lock here and instead have a function that connects the widgets 
//  // with the parameters where we acquire the lock - but maybe not
//
//  // set the plugIn-headline:
//  setHeadlineText("Ladder");
//
//  // assign the pointer to the edited object:
//  jassert(newLadderToEdit != nullptr ); // you must pass a valid module here
//  ladderToEdit = newLadderToEdit;
//
//  // create the widgets and assign the automatable parameters to them:
//
//  frequencyResponseDisplay = new LadderSpectrumEditor("SpectrumEditor");
//  frequencyResponseDisplay->setFilterToEdit(ladderToEdit);
//  frequencyResponseDisplay->addChangeListener(this); // do we need this?
//  frequencyResponseDisplay->assignParameterFreq( moduleToEdit->getParameterByName("Cutoff"));
//  frequencyResponseDisplay->assignParameterReso( moduleToEdit->getParameterByName("Resonance"));
//  addPlot( frequencyResponseDisplay );
//
//  addWidget( cutoffSlider = new RSlider("CutoffSlider") );
//  cutoffSlider->assignParameter( ladderToEdit->getParameterByName("Cutoff") );
//  cutoffSlider->setSliderName("Cutoff");
//  cutoffSlider->setDescription("Cutoff frequency in Hz");
//  cutoffSlider->setDescriptionField(infoField);
//  cutoffSlider->setStringConversionFunction(&hertzToStringWithUnitTotal5);
//
//  addWidget( resonanceSlider = new RSlider("ResoSlider") );
//  resonanceSlider->assignParameter( ladderToEdit->getParameterByName("Resonance") );
//  resonanceSlider->setSliderName("Resonance");
//  resonanceSlider->setDescription("Amount of feedback");
//  resonanceSlider->setDescriptionField(infoField);
//  resonanceSlider->setStringConversionFunction(&valueToStringTotal5);
//
//  addWidget( spreadSlider = new RSlider("SpreadSlider") );
//  spreadSlider->assignParameter( ladderToEdit->getParameterByName("StereoSpread") );
//  spreadSlider->setSliderName("Spread");
//  spreadSlider->setDescription("Detunes cutoff frequencies of channels");
//  spreadSlider->setDescriptionField(infoField);
//  spreadSlider->setStringConversionFunction(&valueToStringTotal5);
//
//  addWidget( modeComboBox = new RComboBox("ModeComboBox") );
//  modeComboBox->assignParameter( ladderToEdit->getParameterByName("Mode") );
//  modeComboBox->setDescription("Select frequency response type");
//  modeComboBox->setDescriptionField(infoField);
//  modeComboBox->registerComboBoxObserver(this); // to update plot when mode is switched
//   // todo: pass the mode parameter to the frequency response plot, too - when use does 
//   // right-click, the popup opens and the mode can be selected. then we don't need to observe
//   // the combobox here...
//   // maybe, we can allow also different modes for the two filters and different resonances
//
//  // change to midSideButton
//  //addWidget( invertButton = new RButton(juce::String(T("Invert"))) );
//  //invertButton->assignParameter( pitchShifterModuleToEdit->getParameterByName(T("Invert")) );
//  //invertButton->setDescription(juce::String(T("Invert polarity of wet (shifted) signal")));
//  //invertButton->setDescriptionField(infoField);
//  //invertButton->setClickingTogglesState(true);
//
//  // set up the widgets:
//  updateWidgetsAccordingToState();
//
//  setSize(420, 240);  
//  // There's a bug somewhere - if we request a size of 400 x 300, the widgets become invisible.
//  // I guess, 400 x 300 is the size, it already has by default, and when we set it again to this 
//  // size, it will be recognized that nothing changed and some update function is not being called
//  // ...or something. ...need to figure that out
//  // ..i think, resized() is not being called and thus, the widget positions won't be set up.
//  // Maybe, the default size should not be 400 x 300 but something like 0 x 0 or (1 x 1, if that's
//  // not allowed) - something that is assured to be never requested as actual size of an editor.
//}
//
//void LadderEditor::resized()
//{
//  ScopedLock scopedLock(*plugInLock);
//  AudioModuleEditor::resized();
//
//  int x  = 0;
//  int y  = 0;
//  int w  = getWidth();
//  int w2 = w/2;
//  int h  = getHeight();
//
//  y = getPresetSectionBottom();
//
//  frequencyResponseDisplay->setBounds(x, y+4, w, h-y-2*20-8); // 2*20 for 2 widget-rows below it
//  y = frequencyResponseDisplay->getBottom();
//
//  cutoffSlider->setBounds(      4, y+4, w2-4, 16);
//  resonanceSlider->setBounds(w2+4, y+4, w2-8, 16);
//
//  y = cutoffSlider->getBottom();  
//
//  modeComboBox->setBounds(   4, y+4, w2-4, 16);
//  spreadSlider->setBounds(w2+4, y+4, w2-8, 16);
//
//  //y = modeComboBox->getBottom(); 
//
//  // maybe here, we somehow have to also resize our wrapper, if any
//}
//
//void LadderEditor::rComboBoxChanged(RComboBox* comboBoxThatHasChanged)
//{
//  frequencyResponseDisplay->updatePlot();
//}