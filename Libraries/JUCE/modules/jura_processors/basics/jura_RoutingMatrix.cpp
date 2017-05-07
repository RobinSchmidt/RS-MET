
//-------------------------------------------------------------------------------------------------
// construction/destruction:

RoutingMatrixAudioModule::RoutingMatrixAudioModule(CriticalSection *newPlugInLock, 
  rosic::RoutingMatrix *newRoutingMatrixToWrap)
: AudioModule(newPlugInLock)
{
  jassert( newRoutingMatrixToWrap != NULL ); // you must pass a valid rosic-object to the constructor
  wrappedRoutingMatrix = newRoutingMatrixToWrap;
  moduleName = juce::String("RoutingMatrix");

  // create and initialize the automatable parameters:
  initializeAutomatableParameters();
}

//-------------------------------------------------------------------------------------------------
// automation:

void RoutingMatrixAudioModule::parameterChanged(Parameter* parameterThatHasChanged)
{
  if( wrappedRoutingMatrix == NULL )
    return;

  double value  = parameterThatHasChanged->getValue();
  int    index  = getIndexOfParameter(parameterThatHasChanged);
  int    input  = index / wrappedRoutingMatrix->getNumOutputs();
  int    output = index % wrappedRoutingMatrix->getNumOutputs();

  wrappedRoutingMatrix->setMatrixEntry(input, output, value);
}

//-------------------------------------------------------------------------------------------------
// internal functions:

void RoutingMatrixAudioModule::initializeAutomatableParameters()
{
  if( wrappedRoutingMatrix == NULL )
    return;

  AutomatableParameter* p;
  for(int i=0; i<wrappedRoutingMatrix->getNumInputs(); i++)
  {
    for(int o=0; o<wrappedRoutingMatrix->getNumOutputs(); o++)
    {
      juce::String name = juce::String("M_") + juce::String(i+1) + 
        juce::String("_") + juce::String(o+1);
      p = new AutomatableParameter(lock, name, -1.0, 1.0, 0.01, 0.0, Parameter::LINEAR);
      addObservedParameter(p);
    }
  }

  for(int i=0; i < (int) parameters.size(); i++ )
    parameterChanged(parameters[i]);
}

//=================================================================================================

// construction/destruction:

RoutingMatrixModuleEditor::RoutingMatrixModuleEditor(CriticalSection *newPlugInLock, 
  RoutingMatrixAudioModule* newRoutingMatrixAudioModule) 
  : AudioModuleEditor(newRoutingMatrixAudioModule)
{
  jassert( newRoutingMatrixAudioModule != NULL ); // you must pass a valid module here

  numInputs  = newRoutingMatrixAudioModule->wrappedRoutingMatrix->getNumInputs();
  numOutputs = newRoutingMatrixAudioModule->wrappedRoutingMatrix->getNumOutputs();

  for(int i=0; i<numInputs; i++)
  {
    RTextField *label = new RTextField(juce::String(i+1));
    rowLabels.add(label);
    addWidget(label);
    label->setJustification(juce::Justification::centred);
    label->setDescription("Input " + juce::String(i+1));
    label->setDescriptionField(infoField);
  }
  for(int i=0; i<numOutputs; i++)
  {
    RTextField *label = new RTextField(juce::String(i+1));
    columnLabels.add(label);
    addWidget(label);
    label->setJustification(juce::Justification::centred);
    label->setDescription("Output " + juce::String(i+1));
    label->setDescriptionField(infoField);
  }
  for(int i=0; i<numInputs; i++)
  {
    for(int o=0; o<numOutputs; o++)
    {
      RDraggableNumber *entryField = new RDraggableNumber(juce::String::empty);
      matrixFields.add(entryField);
      addWidget(entryField);
      entryField->assignParameter( newRoutingMatrixAudioModule->getParameterByName(
        juce::String("M_") + juce::String(i+1) + juce::String("_") + juce::String(o+1) ) );
      entryField->setSliderName(juce::String::empty);
      entryField->setDescription(  juce::String("Amount by which input ") + juce::String(i+1) 
        + juce::String(" goes to output ") + juce::String(o+1) );
      entryField->setDescriptionField(infoField);
      entryField->setStringConversionFunction(valueToString2);
    }
  }

  setLinkPosition(AudioModuleEditor::INVISIBLE);
  stateWidgetSet->setLayout(StateLoadSaveWidgetSet::LABEL_AND_BUTTONS_ABOVE);
  stateWidgetSet->stateLabel->setText(juce::String("Matrix"));

  //setHeadlinePosition(AudioModuleEditor::NO_HEADLINE);

  // set up the widgets:
  updateWidgetsAccordingToState();
}

//-------------------------------------------------------------------------------------------------
// callbacks:

void RoutingMatrixModuleEditor::resized()
{
  AudioModuleEditor::resized();

  int x = 0;
  int y = 0;
  int w = getWidth();
  int h = getHeight();

  stateWidgetSet->setBounds(4, 4, w-8, 32);

  // \todo make the arrangement more general (has been done ad hoc to work for quadrifex)

  x         = 4;
  y         = getPresetSectionBottom();
  w         = 36; // test
  h         = 16;

  for(int i=0; i<numOutputs; i++)
    columnLabels[i]->setBounds(x+i*(w-2)+12, y, w, h);

  y = columnLabels[0]->getBottom()-2;

  for(int i=0; i<numInputs; i++)
    rowLabels[i]->setBounds(x, y+i*(h-2), 14, h);

  int xLeft = rowLabels[0]->getRight()-2;
  y = rowLabels[0]->getY();
  for(int i=0; i<numInputs; i++)
  {
    x = xLeft;
    for(int o=0; o<numOutputs; o++)
    {
      matrixFields[i*numOutputs+o]->setBounds(x, y, w, h);
      x += w-2;
    }
    y += h-2;
  }
}
