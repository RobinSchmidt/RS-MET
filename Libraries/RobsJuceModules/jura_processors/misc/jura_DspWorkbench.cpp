

//-------------------------------------------------------------------------------------------------
// construction/destruction:

DspWorkbenchAudioModule::DspWorkbenchAudioModule(CriticalSection *newPlugInLock,
  rosic::DspWorkbench *dspWorkbenchToWrap) : AudioModule(newPlugInLock)
{
  //jassert(dspWorkbenchToWrap != NULL); // you must pass a valid rosic-object to the constructor

  if(dspWorkbenchToWrap != nullptr)
    wrappedDspWorkbench = dspWorkbenchToWrap;
  else
  {
    wrappedDspWorkbench = new  rosic::DspWorkbench;
    wrappedDspWorkbenchIsOwned = true;
  }
  setModuleTypeName("DSPWorkbench");
  initializeAutomatableParameters();
}

DspWorkbenchAudioModule::~DspWorkbenchAudioModule()
{
  if(wrappedDspWorkbenchIsOwned)
    delete wrappedDspWorkbench;
}

AudioModuleEditor* DspWorkbenchAudioModule::createEditor(int type)
{
  return new DspWorkbenchModuleEditor(lock, this);
}

//-------------------------------------------------------------------------------------------------
// state management:

XmlElement* DspWorkbenchAudioModule::getStateAsXml(const juce::String& stateName, bool markAsClean)
{
  // store the inherited controller mappings:
  XmlElement *xmlState = AudioModule::getStateAsXml(stateName, markAsClean);

  if( wrappedDspWorkbench == NULL )
    return xmlState;

  xmlState->setAttribute("Algorithm", getAlgorithmString());

  // add additional attributes and/or child elements, which are specific to DSP Workbench:
  // store the parameter names, ranges, etc.:
  XmlElement* xmlParameters = new XmlElement(juce::String("Parameters") );
  XmlElement* xmlParameterSettings;
  for(int p=0; p<rosic::DspWorkbench::numParameters; p++)
  {
    xmlParameterSettings = new XmlElement(juce::String("Par") + juce::String(p));
    xmlParameterSettings->setAttribute("Name",
      juce::String(wrappedDspWorkbench->getParameterName(p)));
    xmlParameterSettings->setAttribute("Value",
      juce::String(wrappedDspWorkbench->getParameter(p)));
    /*
    xmlParameterSettings->setAttribute(T("Min"),
      juce::String(wrappedDspWorkbench->getParameterMinimum(p)));
    xmlParameterSettings->setAttribute(T("Max"),
      juce::String(wrappedDspWorkbench->getParameterMaximum(p)));
    */
    xmlParameterSettings->setAttribute("Min", juce::String(getParameterByIndex(p)->getMinValue()));
    xmlParameterSettings->setAttribute("Max", juce::String(getParameterByIndex(p)->getMaxValue()));
    xmlParameterSettings->setAttribute("Default",
      juce::String(getParameterByIndex(p)->getDefaultValue()));
    xmlParameterSettings->setAttribute("Scaling", getParameterByIndex(p)->getScalingString());
    xmlParameters->addChildElement(xmlParameterSettings);
  }
  xmlState->addChildElement(xmlParameters);


  /*
  int i;
  for(i=0; i<wrappedDspWorkbench->numModulators; i++)
  {
    childState = new XmlElement(juce::String(T("BreakpointModulatorState")) + juce::String(i));
    breakpointModulatorStateToXml( &(wrappedDspWorkbench->modulators[0]), childState);
    xmlState->addChildElement(childState);
  }
  */

  return xmlState;
}

void DspWorkbenchAudioModule::setStateFromXml(const XmlElement& xmlState,
                                              const juce::String& stateName, bool markAsClean)
{
  // restore the settings of the inherited AudioModule object:
  AudioModule::setStateFromXml(xmlState, stateName, markAsClean);

  if( wrappedDspWorkbench == NULL )
    return;

  // restore the values of the non-automatable parameters (the automatable ones are already taken
  // care of by automatableModuleStateFromXml():
  //wrappedDspWorkbench->setSwitchWetLeftForRight(xmlState.getBoolAttribute(T("ChannelSwitch"), false));

  // set up the algorithm string:
  setAlgorithmString( xmlState.getStringAttribute("Algorithm", "outL=inL;\noutR=inR") );
  /*
  juce::String algorithmjuce::String = xmlState.getStringAttribute(T("Algorithm"), T("outL=inL;\noutR=inR"));
  char* algorithmjuce::StringC = toZeroTerminatedString(algorithmjuce::String);
  wrappedDspWorkbench->setAlgorithmjuce::String(algorithmjuce::StringC);
  if( algorithmjuce::StringC != NULL )
    delete[] algorithmjuce::StringC;
  */

  // set up the parameters:
  XmlElement* xmlParameters = xmlState.getChildByName("Parameters");
  XmlElement* xmlParameterSettings;
  double min, max, value;
  juce::String nameString;
  char*  nameStringC = NULL;
  if( xmlParameters != NULL )
  {
    for(int p=0; p<rosic::DspWorkbench::numParameters; p++)
    {
      xmlParameterSettings = xmlParameters->getChildByName(juce::String("Par") + juce::String(p));
      if( xmlParameterSettings != NULL )
      {
        nameString =
          xmlParameterSettings->getStringAttribute("Name", juce::String("Par")+juce::String(p));
        nameStringC = toZeroTerminatedString(nameString);
        if( nameStringC != NULL ) // should not go wrong actually, however
        {
          wrappedDspWorkbench->setParameterName(p, nameStringC);
          delete nameStringC;
          nameStringC = NULL;
        }

        min   = xmlParameterSettings->getDoubleAttribute("Min",   0.0);
        max   = xmlParameterSettings->getDoubleAttribute("Max",   1.0);
        value = xmlParameterSettings->getDoubleAttribute("Value", 0.5);

        //wrappedDspWorkbench->setParameterRange(p, min, max);
        wrappedDspWorkbench->setParameter(p, value);
        getParameterByIndex(p)->setRange(min, max);
        getParameterByIndex(p)->setValue(value, true, true);
        getParameterByIndex(p)->setDefaultValue(
          xmlParameterSettings->getDoubleAttribute("Default", 0.5*(min+max)));
        getParameterByIndex(p)->setScalingFromString(
          xmlParameterSettings->getStringAttribute("Scaling", "Linear"));
      }
    }
  }
}

bool DspWorkbenchAudioModule::setAlgorithmString(const juce::String& newAlgorithmString)
{
  algorithmWithComments    = newAlgorithmString;
  algorithmWithoutComments = stripOffComments (algorithmWithComments);

  char* algorithmStringC = toZeroTerminatedString(algorithmWithoutComments);
  bool success = wrappedDspWorkbench->setAlgorithmString(algorithmStringC);
  if( algorithmStringC != NULL )
    delete[] algorithmStringC;

  return success;
}

juce::String DspWorkbenchAudioModule::getAlgorithmString()
{
  return algorithmWithComments;
}

juce::String DspWorkbenchAudioModule::stripOffComments(const juce::String &inputText)
{
  return inputText; //preliminray
}

//-------------------------------------------------------------------------------------------------
// automation:

void DspWorkbenchAudioModule::parameterChanged(Parameter* parameterThatHasChanged)
{
  if( wrappedDspWorkbench == NULL )
    return;

  // find out the index in the vector of the parameter that has been changed:
  int parameterIndex = getIndexOfParameter(parameterThatHasChanged);

  // parameterIndex now contains the index in the array of the parameter that has changed now set
  // up the signal processing:
  double value = parameterThatHasChanged->getValue();

  wrappedDspWorkbench->setParameter(parameterIndex, value);

  //int dummy = 0;
}

//-------------------------------------------------------------------------------------------------
// internal functions:

void DspWorkbenchAudioModule::initializeAutomatableParameters()
{
  // create the automatable parameters and add them to the list - note that the order of the adds
  // is important because in parameterChanged(), the index (position in the array) will be used to
  // identify which particlua parameter has changed.

  // this pointer will be used to temporarily store the addresses of the created
  // Parameter-objects:
  AutomatableParameter* parameter;

  rosic::DspWorkbench* b = wrappedDspWorkbench;
  int p;
  for(p=0; p<rosic::DspWorkbench::numParameters; p++)
  {
    /*
    parameter = new Parameter(b->getParameterName(p), b->getParameterMinimum(p),
      b->getParameterMaximum(p), 0.0, b->getParameterDefault(p) ); //, b->getParameterScaling(), b->getParameterController(), true);
      */
    parameter = new AutomatableParameter(lock, b->getParameterName(p), 0.0, 1.0, 0.0, 0.5,
      Parameter::LINEAR, p);
    addObservedParameter(parameter);
  }

  // make a call to parameterChanged for each parameter in order to set up the DSP-core to reflect
  // the values the automatable parameters:
  for(p=0; p < (int) parameters.size(); p++ )
    parameterChanged(parameters[p]);
}

//=================================================================================================

//-------------------------------------------------------------------------------------------------
// construction/destruction:

DspWorkbenchModuleEditor::DspWorkbenchModuleEditor(CriticalSection *newPlugInLock, DspWorkbenchAudioModule* newDspWorkbenchAudioModule)
  : AudioModuleEditor(newDspWorkbenchAudioModule)
{
  // set the plugIn-headline:
  setHeadlineText( juce::String("DSP Workbench") );

  // assign the pointer to the rosic::DspWorkbench object to be used as aduio engine:
  jassert(dspWorkbenchAudioModule != NULL ); // you must pass a valid module here
  dspWorkbenchAudioModule = newDspWorkbenchAudioModule;



  //-----------------------------------------------------------------------------------------------
  // the widgets for the user parameters:

  for(int p=0; p<rosic::DspWorkbench::numParameters; p++)
  {
    addAndMakeVisible( parameterLabels[p] = new RTextEntryField(
      juce::String("Par")+juce::String(p)+juce::String(":")) );
    parameterLabels[p]->setDescription(juce::String("User parameter ") + juce::String("p"));

    addAndMakeVisible( parameterExpButtons[p] = new RButton(juce::String("Exp")) );
    parameterExpButtons[p]->setDescription(juce::String("Use exponential scaling for this parameter."));
    parameterExpButtons[p]->addRButtonListener(this);

    addAndMakeVisible( parameterSliders[p] = new RSlider("parameterSlider") );
    parameterSliders[p]->addListener(this);
    parameterSliders[p]->setRange(0.0, 1.0, 0.0, 0.5);
    parameterSliders[p]->assignParameter(dspWorkbenchAudioModule->getParameterByIndex(p));
    parameterSliders[p]->setSliderName(juce::String());
    parameterSliders[p]->setDescription(juce::String("Adjusts the value of the parameter"));
    parameterSliders[p]->setStringConversionFunction(&valueToString3);

    addAndMakeVisible( parameterMinFields[p] = new RTextEntryField( juce::String("0.0")) );
    parameterMinFields[p]->registerTextEntryFieldObserver(this);
    parameterMinFields[p]->setDescription("Enter minimum value for the parameter");
    parameterMinFields[p]->setDescriptionField(infoField);

    addAndMakeVisible( parameterMaxFields[p] = new RTextEntryField( juce::String("1.0")) );
    parameterMaxFields[p]->registerTextEntryFieldObserver(this);
    parameterMaxFields[p]->setDescription("Enter maximum value for the parameter");
    parameterMaxFields[p]->setDescriptionField(infoField);

    addAndMakeVisible( parameterDefaultFields[p] = new RTextEntryField( juce::String("0.5")) );
    parameterDefaultFields[p]->registerTextEntryFieldObserver(this);
    parameterDefaultFields[p]->setDescription("Enter default value for the parameter");
    parameterDefaultFields[p]->setDescriptionField(infoField);

    addAndMakeVisible( parameterNameFields[p] = new RTextEntryField( juce::String("Par01")) );
    parameterNameFields[p]->registerTextEntryFieldObserver(this);
    parameterNameFields[p]->setDescription("Enter name value for the parameter");
    parameterNameFields[p]->setDescriptionField(infoField);
  }

  addAndMakeVisible( codeEditor = new RTextEditor(juce::String("CodeEditor")) );
  codeEditor->addListener(this);
  codeEditor->setMultiLine(true, false);
  codeEditor->setReturnKeyStartsNewLine(true);
  //codeEditor->setFont(Font(Font::getDefaultMonospacedFontName(), 14, Font::plain) );
  //codeEditor->setDescription(T("This is the editor for the DSP code."));
  //codeEditor->setDescriptionField(infoField);

  // create the BreakpointModulatorEditor:
  //breakpointModulatorEditor = new BreakpointModulatorEditor( &dspWorkbenchAudioModule->wrappedDspWorkbench->modulators[0] );
  breakpointModulatorEditor = new BreakpointModulatorEditor(lock, NULL);
  breakpointModulatorEditor->addChangeListener(this);
  addAndMakeVisible(breakpointModulatorEditor);

  /*
  sampleModulatorEditor = new SampleModulatorEditor(juce::String(T("SampleModulatorEditor")));
  //sampleModulatorEditor->enablePresetFileManagement(false);
  sampleModulatorEditor->addChangeListener(this);
  addChildComponent(sampleModulatorEditor);
  */

  //---------------------------------------------------------------------------
  // create and setup the filter-section editor:

  //filterEditor = new MultiModeFilterModuleEditor(juce::String(T("FilterEditor")));
  filterEditor = new MultiModeFilterModuleEditor(lock, NULL);
  //filterEditor->addChangeListener(this);
  //filterEditor->setFilterToEdit( &(dspWorkbenchAudioModule->wrappedDspWorkbench->filters[0]) );
  filterEditor->setDescriptionField(infoField, true);
  /*
  filterEditor->freqSlider->assignParameter(
  dspWorkbenchAudioModule->getParameterByName("FilterFrequency") );
  filterEditor->freqKeyTrackSlider->assignParameter(
  dspWorkbenchAudioModule->getParameterByName("FilterFrequencyByKey") );
  filterEditor->freqVelTrackSlider->assignParameter(
  dspWorkbenchAudioModule->getParameterByName("FilterFrequencyByVel") );
  filterEditor->resoSlider->assignParameter(
  dspWorkbenchAudioModule->getParameterByName("FilterResonance") );
  filterEditor->qSlider->assignParameter(
  dspWorkbenchAudioModule->getParameterByName("FilterQ") );
  filterEditor->preAllpassSlider->assignParameter(
  dspWorkbenchAudioModule->getParameterByName("FilterPreAllpass") );
  filterEditor->makeUpSlider->assignParameter(
  dspWorkbenchAudioModule->getParameterByName("FilterMakeUp") );
  filterEditor->gainSlider->assignParameter(
  dspWorkbenchAudioModule->getParameterByName("FilterGain") );
  filterEditor->driveSlider->assignParameter(
  dspWorkbenchAudioModule->getParameterByName("FilterDrive") );
  filterEditor->morphSlider->assignParameter(
  dspWorkbenchAudioModule->getParameterByName("FilterMorph") );
  filterEditor->frequencyResponseDisplay->assignParameterFreq(
  dspWorkbenchAudioModule->getParameterByName("FilterFrequency") );
  filterEditor->frequencyResponseDisplay->assignParameterReso(
  dspWorkbenchAudioModule->getParameterByName("FilterResonance") );
  filterEditor->frequencyResponseDisplay->assignParameterQ(
  dspWorkbenchAudioModule->getParameterByName("FilterQ") );
  filterEditor->frequencyResponseDisplay->assignParameterGain(
  dspWorkbenchAudioModule->getParameterByName("FilterGain") );
  filterEditor->frequencyResponseDisplay->assignParameterMorph(
  dspWorkbenchAudioModule->getParameterByName("FilterMorph") );
  */
  addAndMakeVisible( filterEditor );

  /*
  addAndMakeVisible( sourceComboBox
  = new RComboBox(juce::String(T("ModulationSourceComboBox"))) );
  sourceComboBox->addItem(juce::String(T("Off")),               1);
  sourceComboBox->addItem(juce::String(T("Breakpoints")),       2);
  sourceComboBox->addItem(juce::String(T("Audio-File")),        3);
  //sourceComboBox->addItem(juce::String(T("Envelope-Follower")), 4);
  //sourceComboBox->addItem(juce::String(T("Pitch-Detector")),    5);
  sourceComboBox->setSelectedId(2, true);
  sourceComboBox->addListener(this);
  */

  // set up the widgets:
  updateWidgetsAccordingToState();

  setSize(600, 600);
}

/*
DspWorkbenchModuleEditor::~DspWorkbenchModuleEditor()
{
deleteAllChildren();
}
*/

//-------------------------------------------------------------------------------------------------
// callbacks:

void DspWorkbenchModuleEditor::rButtonClicked(RButton *buttonThatWasClicked)
{
  if( dspWorkbenchAudioModule == NULL )
    return;
  if( dspWorkbenchAudioModule->wrappedDspWorkbench == NULL )
    return;

  for(int p=0; p<rosic::DspWorkbench::numParameters; p++)
  {
    if( buttonThatWasClicked == parameterExpButtons[p] )
    {
      if( parameterExpButtons[p]->getToggleState() == true )
      {
        dspWorkbenchAudioModule->getParameterByIndex(p)->setScaling(Parameter::EXPONENTIAL);
        parameterSliders[p]->setScaling(Parameter::EXPONENTIAL);
      }
      else
      {
        dspWorkbenchAudioModule->getParameterByIndex(p)->setScaling(Parameter::LINEAR);
        parameterSliders[p]->setScaling(Parameter::LINEAR);
      }
      return;
    }
  }

  /*
  //else
  {
  AudioModuleEditor::rButtonClicked(buttonThatWasClicked);
  return;
  }
  */
  dspWorkbenchAudioModule->markStateAsDirty();
}


void DspWorkbenchModuleEditor::comboBoxChanged(ComboBox *comboBoxThatHasChanged)
{
  if( dspWorkbenchAudioModule == NULL )
    return;
  if( dspWorkbenchAudioModule->wrappedDspWorkbench == NULL )
    return;

  dspWorkbenchAudioModule->markStateAsDirty();
}

void DspWorkbenchModuleEditor::rTextFieldChanged(RTextField *textFieldThatHasChanged)
{
  if( dspWorkbenchAudioModule == NULL )
    return;
  if( dspWorkbenchAudioModule->wrappedDspWorkbench == NULL )
    return;
  rosic::DspWorkbench* dspWorkbench = dspWorkbenchAudioModule->wrappedDspWorkbench;

  double min, max, value, defaultValue;
  for(int p=0; p<rosic::DspWorkbench::numParameters; p++)
  {
    if( textFieldThatHasChanged == parameterMinFields[p]
      ||  textFieldThatHasChanged == parameterMaxFields[p] )
    {
      // set up the parameter to the new range:
      min = parameterMinFields[p]->getText().getDoubleValue();
      max = parameterMaxFields[p]->getText().getDoubleValue();
      dspWorkbenchAudioModule->getParameterByIndex(p)->setRange(min, max);

      // re-retrieve the values because changing the range may also change the default- and
      // current value, also min or max may have been invalid and rejected therefore:
      min          = dspWorkbenchAudioModule->getParameterByIndex(p)->getMinValue();
      max          = dspWorkbenchAudioModule->getParameterByIndex(p)->getMaxValue();
      value        = dspWorkbenchAudioModule->getParameterByIndex(p)->getValue();
      defaultValue = dspWorkbenchAudioModule->getParameterByIndex(p)->getDefaultValue();

      // update the fields and the slider according to the range:
      parameterMinFields[p]->setText(    juce::String(min)         );
      parameterMaxFields[p]->setText(    juce::String(max)         );
      parameterDefaultFields[p]->setText(juce::String(defaultValue));
      parameterSliders[p]->setRange(min, max, 0.0, defaultValue,     false);
      parameterSliders[p]->setValue(            value,         false);

      return;
      /*
      dspWorkbench->setParameterRange(p, min, max);
      parameterMinFields[p]->setText(    juce::String(dspWorkbench->getParameterMinimum(p)), false);
      parameterMaxFields[p]->setText(    juce::String(dspWorkbench->getParameterMaximum(p)), false);
      parameterDefaultFields[p]->setText(juce::String(dspWorkbench->getParameterDefault(p)), false);
      parameterSliders[p]->setValue(dspWorkbench->getParameter(p), false, false);
      */
    }
    else if( textFieldThatHasChanged == parameterDefaultFields[p] )
    {
      defaultValue = parameterDefaultFields[p]->getText().getDoubleValue();
      dspWorkbenchAudioModule->getParameterByIndex(p)->setDefaultValue(defaultValue);
      defaultValue = dspWorkbenchAudioModule->getParameterByIndex(p)->getDefaultValue();
      parameterSliders[p]->setDefaultValue(defaultValue);
      parameterDefaultFields[p]->setText(juce::String(defaultValue));
    }
    else if( textFieldThatHasChanged == parameterNameFields[p] )
    {
      juce::String newName = parameterNameFields[p]->getText();
      char* newNameC = toZeroTerminatedString(newName);
      dspWorkbench->setParameterName(p, newNameC);
      if( newNameC != NULL )
        delete[] newNameC;
    }
  }

  dspWorkbenchAudioModule->markStateAsDirty();
}

void DspWorkbenchModuleEditor::rSliderValueChanged(RSlider* sliderThatHasChanged)
{
  if( dspWorkbenchAudioModule == NULL )
    return;
  if( dspWorkbenchAudioModule->wrappedDspWorkbench == NULL )
    return;

  for(int p=0; p<rosic::DspWorkbench::numParameters; p++)
  {
    if( sliderThatHasChanged == parameterSliders[p] )
    {
      dspWorkbenchAudioModule->wrappedDspWorkbench->setParameter(
        p, parameterSliders[p]->getValue());
      return;
    }
  }

  dspWorkbenchAudioModule->markStateAsDirty();
}


void DspWorkbenchModuleEditor::rTextEditorTextChanged(RTextEditor &editor)
{

}

void DspWorkbenchModuleEditor::rTextEditorReturnKeyPressed(RTextEditor &editor)
{
  rTextEditorFocusLost(editor);
}

void DspWorkbenchModuleEditor::rTextEditorEscapeKeyPressed(RTextEditor &editor)
{
  rTextEditorFocusLost(editor);
}

void DspWorkbenchModuleEditor::rTextEditorFocusLost(RTextEditor &editor)
{
  if( dspWorkbenchAudioModule == NULL )
    return;
  if( dspWorkbenchAudioModule->wrappedDspWorkbench == NULL )
    return;

  bool success = dspWorkbenchAudioModule->setAlgorithmString(juce::String(editor.getText()));

  //char* algorithmjuce::StringC = toZeroTerminatedString(editor.getText());
  //bool success = dspWorkbenchAudioModule->wrappedDspWorkbench->setAlgorithmjuce::String(algorithmjuce::StringC);
  //if( algorithmjuce::StringC != NULL )
  //  delete[] algorithmjuce::StringC;

  if( success )
    editor.setColour(TextEditor::backgroundColourId, Colours::white);
  else
    editor.setColour(TextEditor::backgroundColourId, Colours::red.brighter(2.0));

  //setPresetDirty();
}


void DspWorkbenchModuleEditor::updateWidgetsAccordingToState()
{
  if( dspWorkbenchAudioModule == NULL )
    return;
  if( dspWorkbenchAudioModule->wrappedDspWorkbench == NULL )
    return;
  rosic::DspWorkbench* dspWorkbench = dspWorkbenchAudioModule->wrappedDspWorkbench;


  // update the global widgets and automatable sliders:
  AudioModuleEditor::updateWidgetsAccordingToState();

  // update the non-automatable widgets:
  codeEditor->setText(dspWorkbenchAudioModule->getAlgorithmString(), false);

  double min, max, value, defaultValue;
  for(int p=0; p<rosic::DspWorkbench::numParameters; p++)
  {
    // parameterExpButtons...
    min          = dspWorkbenchAudioModule->getParameterByIndex(p)->getMinValue();
    max          = dspWorkbenchAudioModule->getParameterByIndex(p)->getMaxValue();
    value        = dspWorkbenchAudioModule->getParameterByIndex(p)->getValue();
    defaultValue = dspWorkbenchAudioModule->getParameterByIndex(p)->getDefaultValue();
    parameterNameFields[p]->setText(   juce::String(dspWorkbench->getParameterName(p)));
    parameterMinFields[p]->setText(    juce::String(min));
    parameterMaxFields[p]->setText(    juce::String(max));
    parameterDefaultFields[p]->setText(juce::String(defaultValue));

    parameterSliders[p]->setRange(min, max, 0.0, defaultValue, false);
    parameterSliders[p]->setValue(            value,         false);

    bool expScaling
      = dspWorkbenchAudioModule->getParameterByIndex(p)->getScaling() == Parameter::EXPONENTIAL;
    parameterExpButtons[p]->setToggleState(expScaling, false);
    if( expScaling == true )
      parameterSliders[p]->setScaling(Parameter::EXPONENTIAL);
    else
      parameterSliders[p]->setScaling(Parameter::LINEAR);


    /*
    parameterMinFields[p]->setText(    juce::String(dspWorkbench->getParameterMinimum(p)), false);
    parameterMaxFields[p]->setText(    juce::String(dspWorkbench->getParameterMaximum(p)), false);
    parameterDefaultFields[p]->setText(juce::String(dspWorkbench->getParameterDefault(p)), false);

    parameterSliders[p]->setValue(            dspWorkbench->getParameter(p),         false, false);
    */
  }
}

void DspWorkbenchModuleEditor::showOrHideTargetSpecificWidgets()
{
  if( dspWorkbenchAudioModule == NULL )
    return;
  if( dspWorkbenchAudioModule->wrappedDspWorkbench == NULL )
    return;


}

void DspWorkbenchModuleEditor::paint(Graphics &g)
{
  fillRectWithBilinearGradient(g, topLeftRectangle, editorColourScheme.topLeft,  editorColourScheme.topRight,
    editorColourScheme.bottomLeft, editorColourScheme.bottomRight);
  fillRectWithBilinearGradient(g, topRightRectangle, editorColourScheme.topLeft, editorColourScheme.topRight,
    editorColourScheme.bottomLeft, editorColourScheme.bottomRight);
  fillRectWithBilinearGradient(g, bottomLeftRectangle, editorColourScheme.topLeft, editorColourScheme.topRight,
    editorColourScheme.bottomLeft, editorColourScheme.bottomRight);
  fillRectWithBilinearGradient(g, bottomRightRectangle, editorColourScheme.topLeft, editorColourScheme.topRight,
    editorColourScheme.bottomLeft, editorColourScheme.bottomRight);

  g.setColour(editorColourScheme.outline);
  g.drawRect(topLeftRectangle);
  g.drawRect(topRightRectangle);
  g.drawRect(bottomLeftRectangle);
  g.drawRect(bottomRightRectangle);
}

void DspWorkbenchModuleEditor::resized()
{
  int x = 0;
  int y = 0;
  int w = getWidth();
  int h = getHeight();

  //webLink->setBounds(w-112, 0, 112-4, 20);
  //infoLabel->setBounds(0, getHeight()-20, 40, 20);
  //infoField->setBounds(infoLabel->getRight(), getHeight()-20, getWidth()-infoLabel->getRight(),20);

  y = getHeadlineBottom();
  h = (infoField->getY()-y)/2;
  w = getWidth()/2;

  topLeftRectangle.setBounds(   x+4, y,   w-4, h);
  topRightRectangle.setBounds(    w, y,   w-4, h);
  bottomLeftRectangle.setBounds(x+4, y+h, w-4, h);
  bottomRightRectangle.setBounds( w, y+h, w-4, h);


  breakpointModulatorEditor->setBounds(bottomLeftRectangle);

  codeEditor->setBounds(topRightRectangle);

  //---------------------------------------------------------------------------

  x = topLeftRectangle.getX();
  y = topLeftRectangle.getY();
  w = topLeftRectangle.getWidth();
  h = topLeftRectangle.getHeight();
  stateWidgetSet->stateLabel->setBounds(x+4, y+4, 52, 20);
  stateWidgetSet->statePlusButton->setBounds(x+w-4-20, y+4, 20, 20);
  stateWidgetSet->stateMinusButton->setBounds(stateWidgetSet->statePlusButton->getX()-20, y+4, 20, 20);
  stateWidgetSet->stateLoadButton->setBounds(stateWidgetSet->stateMinusButton->getX()-40-4, y+4, 40, 20);
  stateWidgetSet->stateSaveButton->setBounds(stateWidgetSet->stateLoadButton->getX()-40-4, y+4, 40, 20);
  x = stateWidgetSet->stateLabel->getRight();
  w = stateWidgetSet->stateSaveButton->getX()-x;
  stateWidgetSet->stateFileNameLabel->setBounds(x+4, y+4, stateWidgetSet->stateSaveButton->getX()-x-8, 20);

  // set up the bounds of the user parameter sliders:
  w = topLeftRectangle.getWidth()/2;
  x = topLeftRectangle.getX()+w;
  y = stateWidgetSet->stateFileNameLabel->getBottom()+4;
  int w2 = (w-8)/4;
  for(int p=0; p<rosic::DspWorkbench::numParameters; p++)
  {
    parameterSliders[p]->setBounds(                                     x+4,      y+16, w-8,  16);
    parameterLabels[p]->setBounds(       parameterSliders[p]->getX(),                y, 56,   16);
    parameterExpButtons[p]->setBounds(   parameterSliders[p]->getRight()-32,      y+32, 32,   16);
    parameterNameFields[p]->setBounds(   parameterLabels[p]->getRight(), y,
      parameterSliders[p]->getRight()-parameterLabels[p]->getRight(), 16);
    parameterMinFields[p]->setBounds(    parameterSliders[p]->getX(),             y+32, w2-4, 16);
    parameterDefaultFields[p]->setBounds(parameterMinFields[p]->getRight()+4,     y+32, w2-4, 16);
    parameterMaxFields[p]->setBounds(    parameterDefaultFields[p]->getRight()+4, y+32, w2-4, 16);
    y += 52;
  }

}
