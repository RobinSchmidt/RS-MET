
// construction/destruction:

FuncShaperAudioModule::FuncShaperAudioModule(CriticalSection *newPlugInLock,
  rosic::FuncShaper *funcShaperToWrap, MetaParameterManager* metaManagerToUse,
  ModulationManager* modManagerToUse) 
  : ModulatableAudioModule(newPlugInLock, metaManagerToUse, modManagerToUse)
{
  jassert(funcShaperToWrap != NULL); // you must pass a valid rosic-object to the constructor
  wrappedFuncShaper = funcShaperToWrap;
  init();
}

FuncShaperAudioModule::FuncShaperAudioModule(CriticalSection *newPlugInLock, 
  MetaParameterManager* metaManagerToUse, ModulationManager* modManagerToUse)
  : ModulatableAudioModule(newPlugInLock, metaManagerToUse, modManagerToUse)
{
  wrappedFuncShaper = new rosic::FuncShaper;
  wrappedFuncShaperIsOwned = true;
  init();
}

void FuncShaperAudioModule::init()
{
  moduleName  = juce::String("FuncShaper");
  setActiveDirectory(getApplicationDirectory() + juce::String("/Presets/FuncShaper") );
  createParameters();

  // use initial value for "a" that is different from the default value:
  setFormulaParameterMaxValue("a", 4.0);
  getParameterByName("a")->setValue(2.0, true, true);
}

FuncShaperAudioModule::~FuncShaperAudioModule()
{
  if(wrappedFuncShaperIsOwned)
    delete wrappedFuncShaper;
}

AudioModuleEditor* FuncShaperAudioModule::createEditor()
{
  return new jura::FuncShaperModuleEditor(lock, this); // get rid of passing the lock
}

// state management:

XmlElement* FuncShaperAudioModule::getStateAsXml(const juce::String& stateName, bool markAsClean)
{
  XmlElement *xmlState = AudioModule::getStateAsXml(stateName, markAsClean);
  xmlState->setAttribute("FunctionString", juce::String(wrappedFuncShaper->getFunctionString()));
  return xmlState;
}

void FuncShaperAudioModule::setStateFromXml(const XmlElement& xmlState,
                                            const juce::String& stateName, bool markAsClean)
{
  // set up min/max values for a,b,c,d - we can't use inherited behavior because min and max must be set simultaneously:
  setFormulaParameterRange("a", xmlState.getDoubleAttribute("aMin", -1.0), xmlState.getDoubleAttribute("aMax", +1.0));
  setFormulaParameterRange("b", xmlState.getDoubleAttribute("bMin", -1.0), xmlState.getDoubleAttribute("bMax", +1.0));
  setFormulaParameterRange("c", xmlState.getDoubleAttribute("cMin", -1.0), xmlState.getDoubleAttribute("cMax", +1.0));
  setFormulaParameterRange("d", xmlState.getDoubleAttribute("dMin", -1.0), xmlState.getDoubleAttribute("dMax", +1.0));

  // use basclass implementation to restore numeric parameters
  AudioModule::setStateFromXml(xmlState, stateName, markAsClean);

  // restore the function-string:
  juce::String functionString = xmlState.getStringAttribute("FunctionString");
  char* functionStringC = toZeroTerminatedString(functionString);
  //bool stringIsValid = wrappedFuncShaper->setFunctionString(functionStringC, false);
  bool stringIsValid = wrappedFuncShaper->setFunctionString(functionStringC, true);
  if(functionStringC)
    delete functionStringC;

  if( markAsClean == true )
    markStateAsClean();

  sendChangeMessage();
}

//-------------------------------------------------------------------------------------------------
// automation:

void FuncShaperAudioModule::parameterChanged(Parameter* p)
{
  juce::String name = p->getName();
  if(name == "a" || name == "b" || name == "c" || name == "d")
    sendChangeMessage();
  markStateAsDirty();
}

void FuncShaperAudioModule::setFormulaParameterMinValue(const juce::String& name, double newMin)
{
  getParameterByName(name + "Min")->setValue(newMin, false, false);
  Parameter *p = getParameterByName(name);
  p->setMinValue(newMin);
  p->setDefaultValue(0.5 * (p->getMinValue() + p->getMaxValue()) );
    // maybe use a member useMidValueAsDefaultValue in Parameter -> may take mapping into account
    // or maybe use a useProportionAsDefaultValue(double proportion) and call with 0.5
}

void FuncShaperAudioModule::setFormulaParameterMaxValue(const juce::String& name, double newMax)
{
  getParameterByName(name + "Max")->setValue(newMax, false, false);
  Parameter *p = getParameterByName(name);
  p->setMaxValue(newMax);
  p->setDefaultValue(0.5 * (p->getMinValue() + p->getMaxValue()) );
}

void FuncShaperAudioModule::setFormulaParameterRange(const juce::String& name, double newMin, 
  double newMax)
{
  getParameterByName(name + "Min")->setValue(newMin, false, false);
  getParameterByName(name + "Max")->setValue(newMax, false, false);

  Parameter *p = getParameterByName(name);
  p->setRange(newMin, newMax);
  p->setDefaultValue(0.5 * (p->getMinValue() + p->getMaxValue()) );
}

//-------------------------------------------------------------------------------------------------
// internal functions:

void FuncShaperAudioModule::createParameters()
{
  typedef ModulatableParameter Param;
  Param* p;

  typedef rosic::FuncShaper FS;
  FS* fs = wrappedFuncShaper;

  // create automatable parameters:

  p = new Param("InputFilterUsed", 0.0, 1.0, 0.0, Parameter::BOOLEAN);
  //p->setDefaultModParameters(0, 1, 0, 1, ModulationConnection::ABSOLUTE); // does not compile - why?
  p->setValueChangeCallback<FS>(fs, &FS::useInputFilter);
  addObservedParameter(p);

  p = new Param("InputHighpass", 20.0, 20000.0, 20.0, Parameter::EXPONENTIAL);
  //p->setDefaultModParameters(20, 20000, -1, +1, ModulationConnection::EXPONENTIAL);
  p->setValueChangeCallback<FS>(fs, &FS::setInHighpassCutoff);
  addObservedParameter(p);

  p = new Param("InputLowpass", 20.0, 20000.0, 20000.0, Parameter::EXPONENTIAL);
  //p->setDefaultModParameters(20, 20000, -1, +1, ModulationConnection::EXPONENTIAL);
  p->setValueChangeCallback<FS>(fs, &FS::setInLowpassCutoff);
  addObservedParameter(p);

  p = new Param("Drive", -48.0, 48.0, 0.0, Parameter::LINEAR);
  p->setValueChangeCallback<FS>(fs, &FS::setDrive);
  addObservedParameter(p);

  p = new Param("DC", -1.0, 1.0, 0.0, Parameter::LINEAR);
  p->setValueChangeCallback<FS>(fs, &FS::setDcOffset);
  addObservedParameter(p);

  p = new Param("OutputFilterUsed", 0.0, 1.0, 0.0, Parameter::BOOLEAN);
  //p->setDefaultModParameters(0, 1, 0, 1, ModulationConnection::ABSOLUTE); // does not compile - why?
  p->setValueChangeCallback<FS>(fs, &FS::useOutputFilter);
  addObservedParameter(p);

  p = new Param("OutputHighpass", 20.0, 20000.0, 20.0, Parameter::EXPONENTIAL);
  //p->setDefaultModParameters(20, 20000, -1, +1, ModulationConnection::EXPONENTIAL);
  p->setValueChangeCallback<FS>(fs, &FS::setOutHighpassCutoff);
  addObservedParameter(p);

  p = new Param("OutputLowpass", 20.0, 20000.0, 20000.0, Parameter::EXPONENTIAL);
  //p->setDefaultModParameters(20, 20000, -1, +1, ModulationConnection::EXPONENTIAL);
  p->setValueChangeCallback<FS>(fs, &FS::setOutLowpassCutoff);
  addObservedParameter(p);

  p = new Param("OutLevel", -24.0, 24.0, 0.0, Parameter::LINEAR);
  p->setValueChangeCallback<FS>(fs, &FS::setOutVol);
  addObservedParameter(p);

  p = new Param("DryWet", 0.0, 100.0, 100.0, Parameter::LINEAR);
  p->setValueChangeCallback<FS>(fs, &FS::setDryWet);
  addObservedParameter(p);


  // create non-automatable parameters:

  typedef jura::FuncShaperAudioModule FSM;

  Parameter* q;

  q = new Parameter("Oversampling", 1.0, 16.0, 4.0, Parameter::LINEAR, 1.0);
  q->setValueChangeCallback<FS>(fs, &FS::setOversampling);
  addObservedParameter(q);

  q = new Parameter("a", -1.0, +1.0, 0.0, Parameter::LINEAR);
  q->setValueChangeCallback<FSM>(this, &FSM::setA);
  addObservedParameter(q);

  q = new Parameter("b", -1.0, +1.0, 0.0, Parameter::LINEAR);
  q->setValueChangeCallback<FSM>(this, &FSM::setB);
  addObservedParameter(q);

  q = new Parameter("c", -1.0, +1.0, 0.0, Parameter::LINEAR);
  q->setValueChangeCallback<FSM>(this, &FSM::setC);
  addObservedParameter(q);

  q = new Parameter("d", -1.0, +1.0, 0.0, Parameter::LINEAR);
  q->setValueChangeCallback<FSM>(this, &FSM::setD);
  addObservedParameter(q);


  q = new Parameter("aMin", -INF, +INF, -1.0, Parameter::IDENTITY);
  addObservedParameter(q);
  q->setValueChangeCallback<FSM>(this, &FSM::setMinA);

  q = new Parameter("bMin", -INF, +INF, -1.0, Parameter::IDENTITY);
  addObservedParameter(q);
  q->setValueChangeCallback<FSM>(this, &FSM::setMinB);

  q = new Parameter("cMin", -INF, +INF, -1.0, Parameter::IDENTITY);
  addObservedParameter(q);
  q->setValueChangeCallback<FSM>(this, &FSM::setMinC);

  q = new Parameter("dMin", -INF, +INF, -1.0, Parameter::IDENTITY);
  addObservedParameter(q);
  q->setValueChangeCallback<FSM>(this, &FSM::setMinD);

  q = new Parameter("aMax", -INF, +INF, +1.0, Parameter::IDENTITY);
  addObservedParameter(q);
  q->setValueChangeCallback<FSM>(this, &FSM::setMaxA);

  q = new Parameter("bMax", -INF, +INF, +1.0, Parameter::IDENTITY);
  addObservedParameter(q);
  q->setValueChangeCallback<FSM>(this, &FSM::setMaxB);

  q = new Parameter("cMax", -INF, +INF, +1.0, Parameter::IDENTITY);
  addObservedParameter(q);
  q->setValueChangeCallback<FSM>(this, &FSM::setMaxC);

  q = new Parameter("dMax", -INF, +INF, +1.0, Parameter::IDENTITY);
  addObservedParameter(q);
  q->setValueChangeCallback<FSM>(this, &FSM::setMaxD);
}

//=================================================================================================

// construction/destruction:

FuncShaperModuleEditor::FuncShaperModuleEditor(CriticalSection *newPlugInLock,
  FuncShaperAudioModule* newFuncShaperAudioModule)
  : AudioModuleEditor(newFuncShaperAudioModule)
{
  // set the plugIn-headline:
  setHeadlineText( juce::String("FuncShaper") );

  // assign the pointer to the rosic::FuncShaper object to be used as aduio engine:
  jassert(newFuncShaperAudioModule != NULL ); // you must pass a valid module here
  funcShaperAudioModule = newFuncShaperAudioModule;

  addWidget( formulaLabel = new RTextField(juce::String("Formula:")) );
  formulaLabel->setDescription("Expression for waveshaping transfer function");
  formulaLabel->setDescriptionField(infoField);

  addWidget( formulaField = new RTextEntryField( juce::String("tanh(a*x);")) );
  formulaField->registerTextEntryFieldObserver(this);
  formulaField->setDescription(formulaLabel->getDescription());
  formulaField->setDescriptionField(infoField);

  addWidget( aMinField = new RTextEntryField( juce::String("0.0")) );
  aMinField->assignParameter(funcShaperAudioModule->getParameterByName("aMin"));
  aMinField->setDescription("Minimum value for a-parameter");
  aMinField->setDescriptionField(infoField);

  addWidget( aMaxField = new RTextEntryField( juce::String("1.0")) );
  aMaxField->assignParameter(funcShaperAudioModule->getParameterByName("aMax"));
  aMaxField->setDescription("Maximum value for a-parameter");
  aMaxField->setDescriptionField(infoField);

  addWidget( aSlider = new RSlider("aSlider") );
  aSlider->assignParameter(funcShaperAudioModule->getParameterByName("a"));
  aSlider->setSliderName(juce::String("a = "));
  aSlider->setDescription(juce::String("Value of a-parameter in the formula"));
  aSlider->setDescriptionField(infoField);
  aSlider->setStringConversionFunction(&valueToString3);

  addWidget( bMinField = new RTextEntryField( juce::String("0.0")) );
  bMinField->assignParameter(funcShaperAudioModule->getParameterByName("bMin"));
  bMinField->setDescription("Minimum value for b-parameter");
  bMinField->setDescriptionField(infoField);

  addWidget( bMaxField = new RTextEntryField( juce::String("1.0")) );
  bMaxField->assignParameter(funcShaperAudioModule->getParameterByName("bMax"));
  bMaxField->setDescription("Maximum value for b-parameter");
  bMaxField->setDescriptionField(infoField);

  addWidget( bSlider = new RSlider("bSlider") );
  bSlider->assignParameter(funcShaperAudioModule->getParameterByName("b"));
  bSlider->setSliderName(juce::String("b = "));
  bSlider->setDescription(juce::String("Value of b-parameter in the formula"));
  bSlider->setDescriptionField(infoField);
  bSlider->setStringConversionFunction(&valueToString3);

  addWidget( cMinField = new RTextEntryField( juce::String("0.0")) );
  cMinField->assignParameter(funcShaperAudioModule->getParameterByName("cMin"));
  cMinField->setDescription("Minimum value for c-parameter");
  cMinField->setDescriptionField(infoField);

  addWidget( cMaxField = new RTextEntryField( juce::String("1.0")) );
  cMaxField->assignParameter(funcShaperAudioModule->getParameterByName("cMax"));
  cMaxField->setDescription("Maximum value for c-parameter");
  cMaxField->setDescriptionField(infoField);

  addWidget( cSlider = new RSlider("cSlider") );
  cSlider->assignParameter(funcShaperAudioModule->getParameterByName("c"));
  cSlider->setSliderName(juce::String("c = "));
  cSlider->setDescription(juce::String("Value of c-parameter in the formula"));
  cSlider->setDescriptionField(infoField);
  cSlider->setStringConversionFunction(&valueToString3);

  addWidget( dMinField = new RTextEntryField( juce::String("0.0")) );
  dMinField->assignParameter(funcShaperAudioModule->getParameterByName("dMin"));
  dMinField->setDescription("Minimum value for d-parameter");
  dMinField->setDescriptionField(infoField);

  addWidget( dMaxField = new RTextEntryField( juce::String("1.0")) );
  dMaxField->assignParameter(funcShaperAudioModule->getParameterByName("dMax"));
  dMaxField->setDescription("Maximum value for d-parameter");
  dMaxField->setDescriptionField(infoField);

  addWidget( dSlider = new RSlider("dSlider") );
  dSlider->assignParameter(funcShaperAudioModule->getParameterByName("d"));
  dSlider->setSliderName(juce::String("d = "));
  dSlider->setDescription(juce::String("Value of d-parameter in the formula"));
  dSlider->setDescriptionField(infoField);
  dSlider->setStringConversionFunction(&valueToString3);

  addWidget( inputLabel = new RTextField(juce::String("Input Signal")) );
  inputLabel->setDescription("Waveshaper input signal manipulations");
  inputLabel->setDescriptionField(infoField);

  addWidget( preFilterButton = new RButton(juce::String("Filter")) );
  preFilterButton->assignParameter(funcShaperAudioModule->getParameterByName("InputFilterUsed"));
  preFilterButton->setDescription(juce::String("Switch input filter on/off"));
  preFilterButton->setDescriptionField(infoField);
  preFilterButton->setClickingTogglesState(true);

  addWidget( inHighpassSlider = new RSlider ("InHighpassSlider") );
  inHighpassSlider->setRange(20.0, 20000.0, 0.001, 20.0);
  inHighpassSlider->setScaling(Parameter::EXPONENTIAL);
  inHighpassSlider->assignParameter(funcShaperAudioModule->getParameterByName("InputHighpass") );
  inHighpassSlider->setSliderName(juce::String("HPF"));
  inHighpassSlider->setDescription(juce::String("Highpass cutoff for input signal"));
  inHighpassSlider->setDescriptionField(infoField);
  inHighpassSlider->setStringConversionFunction(&hertzToStringWithUnitTotal5);

  addWidget( inLowpassSlider = new RSlider("InLowpassSlider") );
  inLowpassSlider->setRange(20.0, 20000.0, 0.001, 20000.0);
  inLowpassSlider->setScaling(Parameter::EXPONENTIAL);
  inLowpassSlider->assignParameter(
    funcShaperAudioModule->getParameterByName("InputLowpass") );
  inLowpassSlider->setSliderName(juce::String("LPF"));
  inLowpassSlider->setDescription(juce::String("Lowpass cutoff for input signal"));
  inLowpassSlider->setDescriptionField(infoField);
  inLowpassSlider->setStringConversionFunction(&hertzToStringWithUnitTotal5);

  //addWidget( driveSlider = new RSlider("DriveSlider") );
  addWidget( driveSlider = new ModulatableSlider );
  driveSlider->setRange(-48.0, 48.0, 0.01, 0.0);
  driveSlider->assignParameter(funcShaperAudioModule->getParameterByName("Drive") );
  driveSlider->setSliderName(juce::String("Drive"));
  driveSlider->setDescription(juce::String("Gain at the waveshaper's input"));
  driveSlider->setDescriptionField(infoField);
  driveSlider->setStringConversionFunction(&decibelsToStringWithUnit2);

  //addWidget( dcSlider = new RSlider("DCSlider") );
  addWidget( dcSlider = new ModulatableSlider );
  dcSlider->setRange(-1.0, 1.0, 0.001, 0.0);
  dcSlider->assignParameter(funcShaperAudioModule->getParameterByName("DC") );
  dcSlider->setSliderName(juce::String("DC"));
  dcSlider->setDescription(juce::String("DC offset at the waveshaper's input (post \"Drive\")"));
  dcSlider->setDescriptionField(infoField);
  dcSlider->setStringConversionFunction(&valueToString3);

  addWidget( oversamplingSlider = new RSlider("oversamplingSlider") );
  oversamplingSlider->setRange(1.0, 16.0, 1.0, 4.0);
  oversamplingSlider->assignParameter(funcShaperAudioModule->getParameterByName("Oversampling") );
  oversamplingSlider->setSliderName(juce::String("Oversampling"));
  oversamplingSlider->setDescription(juce::String("Oversampling factor for internal processing"));
  oversamplingSlider->setDescriptionField(infoField);
  oversamplingSlider->setStringConversionFunction(&valueToString0);

  addWidget( outputLabel = new RTextField(juce::String("Output Signal")) );
  outputLabel->setDescription("Waveshaper output signal manipulations");
  outputLabel->setDescriptionField(infoField);

  addWidget( postFilterButton = new RButton(juce::String("Filter")) );
  postFilterButton->assignParameter(funcShaperAudioModule->getParameterByName("OutputFilterUsed"));
  postFilterButton->setDescription(juce::String("Switch output filter on/off"));
  postFilterButton->setDescriptionField(infoField);
  postFilterButton->setClickingTogglesState(true);

  addWidget( outHighpassSlider = new RSlider("OutHighpassSlider") );
  outHighpassSlider->setRange(20.0, 20000.0, 0.001, 20.0);
  outHighpassSlider->setScaling(Parameter::EXPONENTIAL);
  outHighpassSlider->assignParameter(
    funcShaperAudioModule->getParameterByName("OutputHighpass") );
  outHighpassSlider->setSliderName(juce::String("HPF"));
  outHighpassSlider->setDescription(juce::String("Highpass cutoff for output signal"));
  outHighpassSlider->setDescriptionField(infoField);
  outHighpassSlider->setStringConversionFunction(&hertzToStringWithUnitTotal5);

  addWidget( outLowpassSlider = new RSlider("OutLowpassSlider") );
  outLowpassSlider->setRange(20.0, 20000.0, 0.001, 20000.0);
  outLowpassSlider->setScaling(Parameter::EXPONENTIAL);
  outLowpassSlider->assignParameter(
    funcShaperAudioModule->getParameterByName("OutputLowpass") );
  outLowpassSlider->setSliderName(juce::String("LPF"));
  outLowpassSlider->setDescription(juce::String("Lowpass cutoff for output signal"));
  outLowpassSlider->setDescriptionField(infoField);
  outLowpassSlider->setStringConversionFunction(&hertzToStringWithUnitTotal5);

  addWidget( outVolumeSlider = new RSlider("OutVolumeSlider") );
  outVolumeSlider->setRange(-24.0, 24.0, 0.01, 0.0);
  outVolumeSlider->assignParameter(funcShaperAudioModule->getParameterByName("OutLevel") );
  outVolumeSlider->setSliderName(juce::String("Volume"));
  outVolumeSlider->setDescription(juce::String("Gain at the waveshaper's output"));
  outVolumeSlider->setDescriptionField(infoField);
  outVolumeSlider->setStringConversionFunction(&decibelsToStringWithUnit2);

  addWidget( dryWetSlider = new RSlider("DryWetSlider") );
  dryWetSlider->setRange(0.0, 100.0, 0.01, 100.0);
  dryWetSlider->assignParameter(funcShaperAudioModule->getParameterByName("DryWet") );
  dryWetSlider->setSliderName(juce::String("Dry/Wet"));
  dryWetSlider->setDescription(juce::String("Mix between clean and distorted signal"));
  dryWetSlider->setDescriptionField(infoField);
  dryWetSlider->setStringConversionFunction(&percentToStringWithUnit1);

  // initialize the plot:
  shaperPlot = new CurveFamilyPlotOld();
  shaperPlot->setCurrentRange(-1.5, +1.5, -1.5, +1.5);
  shaperPlot->setAxisLabelX(juce::String("in"));
  shaperPlot->setAxisLabelY(juce::String("out"));
  shaperPlot->setHorizontalCoarseGrid(1.0, true);
  shaperPlot->setHorizontalFineGrid(0.1, true);
  shaperPlot->setVerticalCoarseGrid(1.0, true);
  shaperPlot->setVerticalFineGrid(0.1, true);
  addPlot(shaperPlot);

  // observe funcShaperAudioModule in order to redraw plot when formula changes:
  funcShaperAudioModule->addChangeListener(this);

  // generate the x-values for the plot:
  int    numValues = funcShaperAudioModule->wrappedFuncShaper->distortionCurve.getTableSize();
  double minX      = funcShaperAudioModule->wrappedFuncShaper->distortionCurve.getLowerLimitX();
  double maxX      = funcShaperAudioModule->wrappedFuncShaper->distortionCurve.getUpperLimitX();
  double slope     = (maxX-minX) / (double) (numValues-1);
  xValues = new double[numValues];
  for(int i=0; i<numValues; i++)
    xValues[i] = minX + slope * (double) i;

  // pass the adresses with x- and y- values over to the plot:
  yValues = new double*[1];
  *yValues = funcShaperAudioModule->wrappedFuncShaper->distortionCurve.getTableAdress();
  shaperPlot->setFunctionFamilyValues(numValues, 1, xValues, yValues);

  // set up the widgets:
  updateWidgetsAccordingToState();

  setSize(480, 300);
}

FuncShaperModuleEditor::~FuncShaperModuleEditor()
{
  funcShaperAudioModule->removeChangeListener(this);
  if( xValues )
    delete[] xValues;
  if( yValues )
    delete[] yValues;
}

//-------------------------------------------------------------------------------------------------
// callbacks:

void FuncShaperModuleEditor::changeListenerCallback(ChangeBroadcaster *source)
{
  if(source == funcShaperAudioModule)
    shaperPlot->updatePlotImage(); 
  else
    AudioModuleEditor::changeListenerCallback(source);
}

void FuncShaperModuleEditor::textChanged(RTextEntryField *rTextEntryFieldThatHasChanged)
{
  if( funcShaperAudioModule == NULL )
    return;
  if( funcShaperAudioModule->wrappedFuncShaper == NULL )
    return;

  if( rTextEntryFieldThatHasChanged == formulaField )
  {
    juce::String functionString  = formulaField->getText();
    char*        functionStringC = toZeroTerminatedString(functionString);
    bool stringIsValid = funcShaperAudioModule->wrappedFuncShaper
      ->setFunctionString(functionStringC, false); // why false, shouldn't we recalculate the table?
    if(functionStringC)
      delete functionStringC;

    // indicate, whether the string was valid or not by the background colour of the editor:
    if( stringIsValid )
    {
      funcShaperAudioModule->wrappedFuncShaper->calculateTable();
      formulaField->markTextAsInvalid(false);
    }
    else
      formulaField->markTextAsInvalid(true);

    shaperPlot->updatePlotImage();
  }

  funcShaperAudioModule->markStateAsDirty();
}

void FuncShaperModuleEditor::updateWidgetsAccordingToState()
{
  AudioModuleEditor::updateWidgetsAccordingToState();
  formulaField->setText(juce::String(funcShaperAudioModule->wrappedFuncShaper
    ->getFunctionString()));
  shaperPlot->updatePlotImage();
}

void FuncShaperModuleEditor::paint(Graphics &g)
{
  AudioModuleEditor::paint(g);

  fillRectWithBilinearGradient(g, formulaRectangle, editorColourScheme.topLeft, 
    editorColourScheme.topRight, editorColourScheme.bottomLeft, editorColourScheme.bottomRight);
  fillRectWithBilinearGradient(g, inputRectangle, editorColourScheme.topLeft,
    editorColourScheme.topRight, editorColourScheme.bottomLeft, editorColourScheme.bottomRight);
  fillRectWithBilinearGradient(g, outputRectangle, editorColourScheme.topLeft, 
    editorColourScheme.topRight, editorColourScheme.bottomLeft, editorColourScheme.bottomRight);

  g.setColour(editorColourScheme.outline);
  g.drawRect(formulaRectangle);
  g.drawRect(inputRectangle);
  g.drawRect(outputRectangle);
}

void FuncShaperModuleEditor::resized()
{
  AudioModuleEditor::resized();
  int x = 0;
  int y = 0;
  int w = getWidth();
  int h = getHeight();

  y = getPresetSectionBottom();
  formulaRectangle.setBounds(4, y+4, getWidth()-8, 64);

  // the top rectangle:

  x = formulaRectangle.getX();
  y = formulaRectangle.getY();
  w = formulaRectangle.getWidth();

  formulaLabel->setBounds(x+4, y+4, 60, 20);
  formulaField->setBounds(formulaLabel->getRight()+4, y+4, w-formulaLabel->getRight()-4, 20);

  y += 24;
  w  = formulaRectangle.getWidth()/4;

  int w2 = w/2-8;
  aSlider->setBounds(x+4, y+4, w-8, 16);
  aMinField->setBounds(aSlider->getX(),        aSlider->getBottom()-2, w2, 16);
  aMaxField->setBounds(aSlider->getRight()-w2, aSlider->getBottom()-2, w2, 16);
  x += w;
  bSlider->setBounds(x+4, y+4, w-8, 16);
  bMinField->setBounds(bSlider->getX(),        bSlider->getBottom()-2, w2, 16);
  bMaxField->setBounds(bSlider->getRight()-w2, bSlider->getBottom()-2, w2, 16);
  x += w;
  cSlider->setBounds(x+4, y+4, w-8, 16);
  cMinField->setBounds(cSlider->getX(),        cSlider->getBottom()-2, w2, 16);
  cMaxField->setBounds(cSlider->getRight()-w2, cSlider->getBottom()-2, w2, 16);
  x += w;
  dSlider->setBounds(x+4, y+4, w-8, 16);
  dMinField->setBounds(dSlider->getX(),        dSlider->getBottom()-2, w2, 16);
  dMaxField->setBounds(dSlider->getRight()-w2, dSlider->getBottom()-2, w2, 16);

  // the plot:

  x = getWidth()/4;
  y = aMinField->getBottom();
  w = getWidth()/2;
  h = w;

  shaperPlot->setBounds(x, y+4, w, h);

  // the left rectangle:

  inputRectangle.setBounds(4, shaperPlot->getY(), shaperPlot->getX()-4, shaperPlot->getHeight() );

  x = inputRectangle.getX();
  y = inputRectangle.getY();
  w = inputRectangle.getWidth();

  inputLabel->setBounds(x+4, y+4, w, 16);
  y += 24;
  preFilterButton->setBounds(x+4, y+4, 48, 16);
  y += 14;
  inHighpassSlider->setBounds(x+4, y+4, w-8, 16);
  y += 14;
  inLowpassSlider->setBounds(x+4, y+4, w-8, 16);
  y += 32;
  driveSlider->setBounds(x+4, y+4, w-8, 16);
  y += 24;
  dcSlider->setBounds(x+4, y+4, w-8, 16);
  y += 32;
  oversamplingSlider->setBounds(x+4, y+4, w-8, 16);

  // the right rectangle:

  outputRectangle.setBounds(shaperPlot->getRight(), shaperPlot->getY() ,
    shaperPlot->getX()-4, shaperPlot->getHeight() );

  x = outputRectangle.getX();
  y = outputRectangle.getY();
  w = outputRectangle.getWidth();

  outputLabel->setBounds(x+4, y+4, w, 16);
  y += 24;
  postFilterButton->setBounds(x+4, y+4, 48, 16);
  y += 14;
  outHighpassSlider->setBounds(x+4, y+4, w-8, 16);
  y += 14;
  outLowpassSlider->setBounds(x+4, y+4, w-8, 16);
  y += 32;
  outVolumeSlider->setBounds(x+4, y+4, w-8, 16);
  y += 32;
  dryWetSlider->setBounds(x+4, y+4, w-8, 16);
}


