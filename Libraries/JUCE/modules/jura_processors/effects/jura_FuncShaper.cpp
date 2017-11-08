
// construction/destruction:

FuncShaperAudioModule::FuncShaperAudioModule(CriticalSection *newPlugInLock,
  rosic::FuncShaper *funcShaperToWrap, MetaParameterManager* metaManagerToUse,
  ModulationManager* modManagerToUse) 
  : ModulatableAudioModule(newPlugInLock, metaManagerToUse, modManagerToUse)
{
  jassert(funcShaperToWrap != NULL); // you must pass a valid rosic-object to the constructor
  wrappedFuncShaper = funcShaperToWrap;


  moduleName  = juce::String("FuncShaper");
  setActiveDirectory(getApplicationDirectory() + juce::String("/Presets/FuncShaper") );
  createParameters();

  // use initial value for "a" that is different from the default value:
  setFormulaParameterMaxValue("aMax", 4.0);
  getParameterByName("a")->setValue(2.0, true, true);
}

FuncShaperAudioModule::FuncShaperAudioModule(CriticalSection *newPlugInLock, 
  MetaParameterManager* metaManagerToUse, ModulationManager* modManagerToUse)
  : ModulatableAudioModule(newPlugInLock, metaManagerToUse, modManagerToUse)
{
  wrappedFuncShaper = new rosic::FuncShaper;
  wrappedFuncShaperIsOwned = true;

  // move to init function, shaer with other cosntructor:
  moduleName  = juce::String("FuncShaper");
  setActiveDirectory(getApplicationDirectory() + juce::String("/Presets/FuncShaper") );
  createParameters();

  // use initial value for "a" that is different from the default value:
  setFormulaParameterMaxValue("aMax", 4.0);
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
  setFormulaParameterRange("a", xmlState.getDoubleAttribute("aMin", 0.0), xmlState.getDoubleAttribute("aMax", 1.0));
  setFormulaParameterRange("b", xmlState.getDoubleAttribute("bMin", 0.0), xmlState.getDoubleAttribute("bMax", 1.0));
  setFormulaParameterRange("c", xmlState.getDoubleAttribute("cMin", 0.0), xmlState.getDoubleAttribute("cMax", 1.0));
  setFormulaParameterRange("d", xmlState.getDoubleAttribute("dMin", 0.0), xmlState.getDoubleAttribute("dMax", 1.0));

  // restore the function-string:
  juce::String functionString = xmlState.getStringAttribute("FunctionString");
  char* functionStringC = toZeroTerminatedString(functionString);
  //bool stringIsValid = wrappedFuncShaper->setFunctionString(functionStringC, false);
  if(functionStringC)
    delete functionStringC;

  // use basclass implementation to restore numeric parameters
  AudioModule::setStateFromXml(xmlState, stateName, markAsClean);

  if( markAsClean == true )
    markStateAsClean();
}

//-------------------------------------------------------------------------------------------------
// automation:

void FuncShaperAudioModule::parameterChanged(Parameter* parameterThatHasChanged)
{
  // function will be obsolete after completing updating the parameter handling

  if( wrappedFuncShaper == NULL )
    return;

  int    index = getIndexOfParameter(parameterThatHasChanged);
  double value = parameterThatHasChanged->getValue();
  switch( index )
  {
  //case   0: wrappedFuncShaper->useInputFilter(value >= 0.5);  break;
  //case   1: wrappedFuncShaper->setInHighpassCutoff(value);    break;
  //case   2: wrappedFuncShaper->setInLowpassCutoff(value);     break;
  //case   3: wrappedFuncShaper->setDrive(value);               break;
  //case   4: wrappedFuncShaper->setDcOffset(value);            break;
  //case   5: wrappedFuncShaper->useOutputFilter(value >= 0.5); break;
  //case   6: wrappedFuncShaper->setOutHighpassCutoff(value);   break;
  //case   7: wrappedFuncShaper->setOutLowpassCutoff(value);    break;
  //case   8: wrappedFuncShaper->setOutVol(value);              break;
  //case   9: wrappedFuncShaper->setDryWet(value);              break;

  //case  10: wrappedFuncShaper->setOversampling((int)value);   break;
  //case  11: wrappedFuncShaper->setA(value, true);             break;
  //case  12: wrappedFuncShaper->setB(value, true);             break;
  //case  13: wrappedFuncShaper->setC(value, true);             break;
  //case  14: wrappedFuncShaper->setD(value, true);             break;
  //default:
  //  {
  //    if( index <= 18 )        // handle changes of min-values
  //      setFormulaParameterMinValue(parameterThatHasChanged->getName(), value);
  //    else if( index <= 22 )   // handle changes of max-values
  //      setFormulaParameterMaxValue(parameterThatHasChanged->getName(), value);
  //    else
  //      DEBUG_BREAK; // unknown parameter
  //  }
  } // end of switch( parameterIndex )


  markStateAsDirty();
}


void FuncShaperAudioModule::setFormulaParameterMinValue(const juce::String& augmentedName, double newMinValue)
{
  /*
  getParameterByName(augmentedName)->setValue(newMinValue, false, false);
  Parameter *p = getParameterByName(augmentedName.substring(0, 1));
  p->setMinValue(newMinValue);
  p->setDefaultValue(0.5 * (p->getMinValue() + p->getMaxValue()) );
    // maybe use a member useMidValueAsDefaultValue in Parameter -> may take mapping into account
    // or maybe use a useProportionAsDefaultValue(double proportion) and call with 0.5
    */
}

void FuncShaperAudioModule::setFormulaParameterMaxValue(const juce::String& augmentedName, double newMaxValue)
{
  /*
  getParameterByName(augmentedName)->setValue(newMaxValue, false, false);
  Parameter *p = getParameterByName(augmentedName.substring(0, 1));
  p->setMaxValue(newMaxValue);
  p->setDefaultValue(0.5 * (p->getMinValue() + p->getMaxValue()) );
  */
}

void FuncShaperAudioModule::setFormulaParameterRange(const juce::String& augmentedName, double newMinValue, double newMaxValue)
{
  /*
  getParameterByName(augmentedName.substring(0, 1) + "Min")->setValue(newMinValue, false, false);
  getParameterByName(augmentedName.substring(0, 1) + "Max")->setValue(newMaxValue, false, false);
  Parameter *p = getParameterByName(augmentedName.substring(0, 1));
  p->setRange(newMinValue, newMaxValue);
  p->setDefaultValue(0.5 * (p->getMinValue() + p->getMaxValue()) );
  */
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

  // old:
  //addObservedParameter(new AutomatableParameter(lock, "InputFilterUsed",   0.0,     1.0, 1.0,     0.0, Parameter::BOOLEAN    ));
  //addObservedParameter(new AutomatableParameter(lock, "InputHighpass",    20.0, 20000.0, 0.0,    20.0, Parameter::EXPONENTIAL));
  //addObservedParameter(new AutomatableParameter(lock, "InputLowpass",     20.0, 20000.0, 0.0, 20000.0, Parameter::EXPONENTIAL));
  //addObservedParameter(new AutomatableParameter(lock, "Drive",           -48.0,    48.0, 0.0,     0.0, Parameter::LINEAR)     );
  //addObservedParameter(new AutomatableParameter(lock, "DC",               -1.0,     1.0, 0.0,     0.0, Parameter::LINEAR)     );
  //addObservedParameter(new AutomatableParameter(lock, "OutputFilterUsed",  0.0,     1.0, 1.0,     0.0, Parameter::BOOLEAN    ));
  //addObservedParameter(new AutomatableParameter(lock, "OutputHighpass",   20.0, 20000.0, 0.0,    20.0, Parameter::EXPONENTIAL));
  //addObservedParameter(new AutomatableParameter(lock, "OutputLowpass",    20.0, 20000.0, 0.0, 20000.0, Parameter::EXPONENTIAL));
  //addObservedParameter(new AutomatableParameter(lock, "OutLevel",        -24.0,    24.0, 0.0,     0.0, Parameter::LINEAR)     );
  //addObservedParameter(new AutomatableParameter(lock, "DryWet",            0.0,   100.0, 0.0,   100.0, Parameter::LINEAR)     );


  // create non-automatable parameters:

  typedef jura::FuncShaperAudioModule FSM;

  Parameter* q;

  q = new Parameter("Oversampling", 1.0, 16.0, 4.0, Parameter::LINEAR, 1.0);
  q->setValueChangeCallback<FS>(fs, &FS::setOversampling);
  addObservedParameter(q);

  q = new Parameter("a", 0.0, 1.0, 0.5, Parameter::LINEAR, 0.01);
  q->setValueChangeCallback<FSM>(this, &FSM::setA);
  addObservedParameter(q);

  q = new Parameter("b", 0.0, 1.0, 0.5, Parameter::LINEAR, 0.01);
  q->setValueChangeCallback<FSM>(this, &FSM::setB);
  addObservedParameter(q);

  q = new Parameter("c", 0.0, 1.0, 0.5, Parameter::LINEAR, 0.01);
  q->setValueChangeCallback<FSM>(this, &FSM::setC);
  addObservedParameter(q);

  q = new Parameter("d", 0.0, 1.0, 0.5, Parameter::LINEAR, 0.01);
  q->setValueChangeCallback<FSM>(this, &FSM::setD);
  addObservedParameter(q);


  q = new Parameter("aMin", -INF, +INF, 0.0, Parameter::LINEAR);
  addObservedParameter(q);
  q->setValueChangeCallback<FSM>(this, &FSM::setMinA);

  q = new Parameter("bMin", -INF, +INF, 0.0, Parameter::LINEAR);
  addObservedParameter(q);
  q->setValueChangeCallback<FSM>(this, &FSM::setMinB);

  q = new Parameter("cMin", -INF, +INF, 0.0, Parameter::LINEAR);
  addObservedParameter(q);
  q->setValueChangeCallback<FSM>(this, &FSM::setMinC);

  q = new Parameter("dMin", -INF, +INF, 0.0, Parameter::LINEAR);
  addObservedParameter(q);
  q->setValueChangeCallback<FSM>(this, &FSM::setMinD);

  q = new Parameter("aMax", -INF, +INF, 0.0, Parameter::LINEAR);
  addObservedParameter(q);
  q->setValueChangeCallback<FSM>(this, &FSM::setMaxA);

  q = new Parameter("bMax", -INF, +INF, 0.0, Parameter::LINEAR);
  addObservedParameter(q);
  q->setValueChangeCallback<FSM>(this, &FSM::setMaxB);

  q = new Parameter("cMax", -INF, +INF, 0.0, Parameter::LINEAR);
  addObservedParameter(q);
  q->setValueChangeCallback<FSM>(this, &FSM::setMaxC);

  q = new Parameter("dMax", -INF, +INF, 0.0, Parameter::LINEAR);
  addObservedParameter(q);
  q->setValueChangeCallback<FSM>(this, &FSM::setMaxD);


  // old:
  //addObservedParameter(new Parameter(lock, "Oversampling", 1.0, 16.0, 1.0, 4.0, Parameter::LINEAR));
  //addObservedParameter(new Parameter(lock, "a",     0.0, 1.0, 0.01, 0.5, Parameter::LINEAR));
  //addObservedParameter(new Parameter(lock, "b",     0.0, 1.0, 0.01, 0.5, Parameter::LINEAR));
  //addObservedParameter(new Parameter(lock, "c",     0.0, 1.0, 0.01, 0.5, Parameter::LINEAR));
  //addObservedParameter(new Parameter(lock, "d",     0.0, 1.0, 0.01, 0.5, Parameter::LINEAR));
  //addObservedParameter(new Parameter(lock, "aMin", -INF, INF, 0.0, 0.0, Parameter::LINEAR));
  //addObservedParameter(new Parameter(lock, "bMin", -INF, INF, 0.0, 0.0, Parameter::LINEAR));
  //addObservedParameter(new Parameter(lock, "cMin", -INF, INF, 0.0, 0.0, Parameter::LINEAR));
  //addObservedParameter(new Parameter(lock, "dMin", -INF, INF, 0.0, 0.0, Parameter::LINEAR));
  //addObservedParameter(new Parameter(lock, "aMax", -INF, INF, 0.0, 1.0, Parameter::LINEAR));
  //addObservedParameter(new Parameter(lock, "bMax", -INF, INF, 0.0, 1.0, Parameter::LINEAR));
  //addObservedParameter(new Parameter(lock, "cMax", -INF, INF, 0.0, 1.0, Parameter::LINEAR));
  //addObservedParameter(new Parameter(lock, "dMax", -INF, INF, 0.0, 1.0, Parameter::LINEAR));

  //// make a call to parameterChanged for each parameter in order to set up the DSP-core to reflect
  //// the values the automatable parameters:
  //for(int i=0; i < (int) parameters.size(); i++ )
  //  parameterChanged(parameters[i]);
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

  // we observe the a,b,c,d parameters in order to redraw the plot when one of them changes:
  funcShaperAudioModule->getParameterByName("a")->registerParameterObserver(this);
  funcShaperAudioModule->getParameterByName("b")->registerParameterObserver(this);
  funcShaperAudioModule->getParameterByName("c")->registerParameterObserver(this);
  funcShaperAudioModule->getParameterByName("d")->registerParameterObserver(this);
  setLocalAutomationSwitch(true);
  setIsGuiElement(true);

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
  funcShaperAudioModule->getParameterByName("a")->deRegisterParameterObserver(this);
  funcShaperAudioModule->getParameterByName("b")->deRegisterParameterObserver(this);
  funcShaperAudioModule->getParameterByName("c")->deRegisterParameterObserver(this);
  funcShaperAudioModule->getParameterByName("d")->deRegisterParameterObserver(this);

  if( xValues )
    delete[] xValues;
  if( yValues )
    delete[] yValues;
  //deleteAllChildren();
}

//-------------------------------------------------------------------------------------------------
// callbacks:

void FuncShaperModuleEditor::parameterIsGoingToBeDeleted(Parameter* parameterThatWillBeDeleted)
{

}

void FuncShaperModuleEditor::parameterChanged(Parameter* parameterThatHasChanged)
{
  shaperPlot->updatePlotImage();  // it was one of the a,b,c,d parameters (these are the ones we observe)
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
    bool stringIsValid = funcShaperAudioModule->wrappedFuncShaper->setFunctionString(functionStringC, false);
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
  formulaField->setText(juce::String(funcShaperAudioModule->wrappedFuncShaper->getFunctionString()));
  shaperPlot->updatePlotImage();
}

void FuncShaperModuleEditor::paint(Graphics &g)
{
  AudioModuleEditor::paint(g);

  fillRectWithBilinearGradient(g, formulaRectangle, editorColourScheme.topLeft, editorColourScheme.topRight,
    editorColourScheme.bottomLeft, editorColourScheme.bottomRight);
  fillRectWithBilinearGradient(g, inputRectangle, editorColourScheme.topLeft,
    editorColourScheme.topRight, editorColourScheme.bottomLeft, editorColourScheme.bottomRight);
  fillRectWithBilinearGradient(g, outputRectangle, editorColourScheme.topLeft, editorColourScheme.topRight,
    editorColourScheme.bottomLeft, editorColourScheme.bottomRight);

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


