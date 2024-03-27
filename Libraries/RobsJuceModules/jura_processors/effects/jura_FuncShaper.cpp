
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
  setModuleTypeName("FuncShaper");
  createParameters();

  // use initial value for "a" that is different from the default value:
  setFormulaParameterMax("a", 4.0);
  getParameterByName("a")->setValue(2.0, true, true);
}

FuncShaperAudioModule::~FuncShaperAudioModule()
{
  if(wrappedFuncShaperIsOwned)
    delete wrappedFuncShaper;
}

AudioModuleEditor* FuncShaperAudioModule::createEditor(int type)
{
  return new jura::FuncShaperModuleEditor(lock, this); // get rid of passing the lock
}

// state management:

XmlElement* FuncShaperAudioModule::getStateAsXml(const juce::String& stateName, bool markAsClean)
{
  XmlElement *xmlState = AudioModule::getStateAsXml(stateName, markAsClean);
  //juce::String str = xmlState->toString();  // For debug
  xmlState->setAttribute("FunctionString", juce::String(wrappedFuncShaper->getFunctionString()));
  return xmlState;
}

void FuncShaperAudioModule::setStateFromXml(const XmlElement& xml,
                                            const juce::String& stateName, bool markAsClean)
{
  autoRecalcTable = false;

  // restore the function-string:
  juce::String functionString = xml.getStringAttribute("FunctionString");
  char* functionStringC = toZeroTerminatedString(functionString);
  //bool stringIsValid = wrappedFuncShaper->setFunctionString(functionStringC, false);
  wrappedFuncShaper->setFunctionString(functionStringC, false);
  if(functionStringC)
    delete functionStringC;

  // recall formula a,b,c,d parameters (and their ranges):
  recallFormulaParameterFromXml(xml, "a");
  recallFormulaParameterFromXml(xml, "b");
  recallFormulaParameterFromXml(xml, "c");
  recallFormulaParameterFromXml(xml, "d");

  wrappedFuncShaper->calculateTable();
  autoRecalcTable = true;

  // use basclass implementation to restore other numeric parameters:
  AudioModule::setStateFromXml(xml, stateName, markAsClean);

  if( markAsClean == true )
    markStateAsClean();

  sendChangeMessage();
}

void FuncShaperAudioModule::recallFormulaParameterFromXml(const XmlElement& xml, 
  const juce::String& name)
{
  // Retrieve formula parameters (e.g. "a", "aMin", "aMax", "b", "bMin", ...):
  double val = xml.getDoubleAttribute(name,          0.0);
  double min = xml.getDoubleAttribute(name + "Min", -1.0);
  double max = xml.getDoubleAttribute(name + "Max", +1.0);

  // Sanitize and set up:
  max = RAPT::rsMax(min, max);          // Enforce max >= min
  val = RAPT::rsClip(val, min, max);    // Enforce min <= val <= max
  setFormulaParameterAndRange(name, val, min, max);
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

void FuncShaperAudioModule::setFormulaParameterMin(const juce::String& name, double newMin)
{
  getParameterByName(name + "Min")->setValue(newMin, false, false);
  Parameter *p = getParameterByName(name);
  p->setMinValue(newMin);
  p->setDefaultValue(0.5 * (p->getMinValue() + p->getMaxValue()) );
    // maybe use a member useMidValueAsDefaultValue in Parameter -> may take mapping into account
    // or maybe use a useProportionAsDefaultValue(double proportion) and call with 0.5
}

void FuncShaperAudioModule::setFormulaParameterMax(const juce::String& name, double newMax)
{
  getParameterByName(name + "Max")->setValue(newMax, false, false);
  Parameter *p = getParameterByName(name);
  p->setMaxValue(newMax);
  p->setDefaultValue(0.5 * (p->getMinValue() + p->getMaxValue()) );
}

void FuncShaperAudioModule::setFormulaParameterAndRange(const juce::String& name, double newValue,
  double newMin, double newMax)
{
  Parameter *p = getParameterByName(name);
  p->setRangeAndValue(newMin, newMax, newValue, true, true);
  p->setDefaultValue(0.5 * (p->getMinValue() + p->getMaxValue()) );
  getParameterByName(name + "Min")->setValue(newMin, false, false);
  getParameterByName(name + "Max")->setValue(newMax, false, false);
}

//-------------------------------------------------------------------------------------------------
// internal functions:

void FuncShaperAudioModule::createParameters()
{
  //typedef ModulatableParameter Param;
  typedef MetaControlledParameter Param;
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
  // assign the pointer to the rosic::FuncShaper object to be used as aduio engine:
  jassert(newFuncShaperAudioModule != NULL ); // you must pass a valid module here
  funcShaperAudioModule = newFuncShaperAudioModule;
  createWidgets();
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

void FuncShaperModuleEditor::createWidgets()
{
  typedef rsAutomatableSlider Sld;
  //typedef rsAutomatableComboBox Box;
  typedef rsAutomatableButton Btn;
  typedef RTextField Lbl;
  Sld* s;
  //Box* c;
  Btn* b;
  Lbl* l;


  addWidget( formulaLabel = l = new Lbl("Formula:") );
  l->setDescription("Expression for waveshaping transfer function");
  l->setDescriptionField(infoField);

  addWidget( formulaField = new RTextEntryField("tanh(a*x);") );
  formulaField->registerTextEntryFieldObserver(this);
  formulaField->setDescription(formulaLabel->getDescription());
  formulaField->setDescriptionField(infoField);

  addWidget( aMinField = new RTextEntryField("0.0") );
  aMinField->assignParameter(funcShaperAudioModule->getParameterByName("aMin"));
  aMinField->setDescription("Minimum value for a-parameter");
  aMinField->setDescriptionField(infoField);

  addWidget( aMaxField = new RTextEntryField("1.0") );
  aMaxField->assignParameter(funcShaperAudioModule->getParameterByName("aMax"));
  aMaxField->setDescription("Maximum value for a-parameter");
  aMaxField->setDescriptionField(infoField);

  addWidget( aSlider = s = new Sld );
  s->assignParameter(funcShaperAudioModule->getParameterByName("a"));
  s->setSliderName("a = ");
  s->setDescription("Value of a-parameter in the formula");
  s->setDescriptionField(infoField);
  s->setStringConversionFunction(&valueToString3);

  addWidget( bMinField = new RTextEntryField("0.0") );
  bMinField->assignParameter(funcShaperAudioModule->getParameterByName("bMin"));
  bMinField->setDescription("Minimum value for b-parameter");
  bMinField->setDescriptionField(infoField);

  addWidget( bMaxField = new RTextEntryField("1.0") );
  bMaxField->assignParameter(funcShaperAudioModule->getParameterByName("bMax"));
  bMaxField->setDescription("Maximum value for b-parameter");
  bMaxField->setDescriptionField(infoField);

  addWidget( bSlider = s = new Sld );
  s->assignParameter(funcShaperAudioModule->getParameterByName("b"));
  s->setSliderName("b = ");
  s->setDescription("Value of b-parameter in the formula");
  s->setDescriptionField(infoField);
  s->setStringConversionFunction(&valueToString3);

  addWidget( cMinField = new RTextEntryField("0.0") );
  cMinField->assignParameter(funcShaperAudioModule->getParameterByName("cMin"));
  cMinField->setDescription("Minimum value for c-parameter");
  cMinField->setDescriptionField(infoField);

  addWidget( cMaxField = new RTextEntryField( juce::String("1.0")) );
  cMaxField->assignParameter(funcShaperAudioModule->getParameterByName("cMax"));
  cMaxField->setDescription("Maximum value for c-parameter");
  cMaxField->setDescriptionField(infoField);

  addWidget( cSlider = s = new Sld );
  s->assignParameter(funcShaperAudioModule->getParameterByName("c"));
  s->setSliderName("c = ");
  s->setDescription("Value of c-parameter in the formula");
  s->setDescriptionField(infoField);
  s->setStringConversionFunction(&valueToString3);

  addWidget( dMinField = new RTextEntryField("0.0") );
  dMinField->assignParameter(funcShaperAudioModule->getParameterByName("dMin"));
  dMinField->setDescription("Minimum value for d-parameter");
  dMinField->setDescriptionField(infoField);

  addWidget( dMaxField = new RTextEntryField("1.0") );
  dMaxField->assignParameter(funcShaperAudioModule->getParameterByName("dMax"));
  dMaxField->setDescription("Maximum value for d-parameter");
  dMaxField->setDescriptionField(infoField);

  addWidget( dSlider = s = new Sld );
  s->assignParameter(funcShaperAudioModule->getParameterByName("d"));
  s->setSliderName("d = ");
  s->setDescription("Value of d-parameter in the formula");
  s->setDescriptionField(infoField);
  s->setStringConversionFunction(&valueToString3);

  addWidget( inputLabel = new RTextField("Input Signal") );
  inputLabel->setDescription("Waveshaper input signal manipulations");
  inputLabel->setDescriptionField(infoField);

  addWidget( preFilterButton = b = new Btn("Filter") );
  //b->setPainter(&buttonPainter);  // test
  b->assignParameter(funcShaperAudioModule->getParameterByName("InputFilterUsed"));
  b->setDescription("Switch input filter on/off");
  b->setDescriptionField(infoField);
  b->setClickingTogglesState(true);

  addWidget( inHighpassSlider = s = new Sld );
  s->setRange(20.0, 20000.0, 0.001, 20.0);
  s->setScaling(Parameter::EXPONENTIAL);
  s->assignParameter(funcShaperAudioModule->getParameterByName("InputHighpass") );
  s->setSliderName("HPF");
  s->setDescription("Highpass cutoff for input signal");
  s->setDescriptionField(infoField);
  s->setStringConversionFunction(&hertzToStringWithUnitTotal5);

  addWidget( inLowpassSlider = s = new Sld );
  s->setRange(20.0, 20000.0, 0.001, 20000.0);
  s->setScaling(Parameter::EXPONENTIAL);
  s->assignParameter(funcShaperAudioModule->getParameterByName("InputLowpass") );
  s->setSliderName("LPF");
  s->setDescription("Lowpass cutoff for input signal");
  s->setDescriptionField(infoField);
  s->setStringConversionFunction(&hertzToStringWithUnitTotal5);

  addWidget( driveSlider = s = new Sld );
  s->setRange(-48.0, 48.0, 0.01, 0.0);
  s->assignParameter(funcShaperAudioModule->getParameterByName("Drive") );
  s->setSliderName("Drive");
  s->setDescription("Gain at the waveshaper's input");
  s->setDescriptionField(infoField);
  s->setStringConversionFunction(&decibelsToStringWithUnit2);

  addWidget( dcSlider = s = new Sld );
  s->setRange(-1.0, 1.0, 0.001, 0.0);
  s->assignParameter(funcShaperAudioModule->getParameterByName("DC") );
  s->setSliderName("DC");
  s->setDescription("DC offset at the waveshaper's input (post \"Drive\")");
  s->setDescriptionField(infoField);
  s->setStringConversionFunction(&valueToString3);

  addWidget( oversamplingSlider = s = new Sld );
  s->setRange(1.0, 16.0, 1.0, 4.0);
  s->assignParameter(funcShaperAudioModule->getParameterByName("Oversampling") );
  s->setSliderName("Oversampling");
  s->setDescription("Oversampling factor for internal processing");
  s->setDescriptionField(infoField);
  s->setStringConversionFunction(&valueToString0);

  addWidget( outputLabel = new Lbl("Output Signal") );
  outputLabel->setDescription("Waveshaper output signal manipulations");
  outputLabel->setDescriptionField(infoField);

  addWidget( postFilterButton = b = new Btn("Filter") );
  b->assignParameter(funcShaperAudioModule->getParameterByName("OutputFilterUsed"));
  b->setDescription("Switch output filter on/off");
  b->setDescriptionField(infoField);
  b->setClickingTogglesState(true);

  addWidget( outHighpassSlider = s = new Sld );
  s->setRange(20.0, 20000.0, 0.001, 20.0);
  s->setScaling(Parameter::EXPONENTIAL);
  s->assignParameter( funcShaperAudioModule->getParameterByName("OutputHighpass") );
  s->setSliderName("HPF");
  s->setDescription("Highpass cutoff for output signal");
  s->setDescriptionField(infoField);
  s->setStringConversionFunction(&hertzToStringWithUnitTotal5);

  addWidget( outLowpassSlider = s = new Sld );
  s->setRange(20.0, 20000.0, 0.001, 20000.0);
  s->setScaling(Parameter::EXPONENTIAL);
  s->assignParameter( funcShaperAudioModule->getParameterByName("OutputLowpass") );
  s->setSliderName("LPF");
  s->setDescription("Lowpass cutoff for output signal");
  s->setDescriptionField(infoField);
  s->setStringConversionFunction(&hertzToStringWithUnitTotal5);

  addWidget( outVolumeSlider = s = new Sld );
  s->setRange(-24.0, 24.0, 0.01, 0.0);
  s->assignParameter(funcShaperAudioModule->getParameterByName("OutLevel") );
  s->setSliderName("Volume");
  s->setDescription("Gain at the waveshaper's output");
  s->setDescriptionField(infoField);
  s->setStringConversionFunction(&decibelsToStringWithUnit2);

  addWidget( dryWetSlider = s = new Sld );
  s->setRange(0.0, 100.0, 0.01, 100.0);
  s->assignParameter(funcShaperAudioModule->getParameterByName("DryWet") );
  s->setSliderName("Dry/Wet");
  s->setDescription("Mix between clean and distorted signal");
  s->setDescriptionField(infoField);
  s->setStringConversionFunction(&percentToStringWithUnit1);

  // initialize the plot:
  shaperPlot = new rsDataPlot();
  shaperPlot->setCurrentRange(-1.5, +1.5, -1.5, +1.5);
  shaperPlot->setAxisLabelX("in");
  shaperPlot->setAxisLabelY("out");
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
  formulaRectangle.setBounds(0, y+4, getWidth(), 64);

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

  inputRectangle.setBounds(0, shaperPlot->getY(), shaperPlot->getX()+2, shaperPlot->getHeight() );

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

  outputRectangle.setBounds(shaperPlot->getRight()-2, shaperPlot->getY(),
    shaperPlot->getX(), shaperPlot->getHeight() );

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
