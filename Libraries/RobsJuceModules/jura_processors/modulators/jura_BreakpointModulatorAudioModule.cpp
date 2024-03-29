//-------------------------------------------------------------------------------------------------
// construction/destruction:

BreakpointModulatorAudioModule::BreakpointModulatorAudioModule(CriticalSection *newPlugInLock,
  rosic::BreakpointModulator *newBreakpointModulatorToWrap)
  : /*AudioModule(newPlugInLock)*/ AudioModuleWithMidiIn(newPlugInLock)
{
  //jassert( newBreakpointModulatorToWrap != NULL ); // you must pass a valid object to the constructor
  if(newBreakpointModulatorToWrap == nullptr)
  {
    newBreakpointModulatorToWrap = new rosic::BreakpointModulator;
    wrappedBreakpointModulatorIsOwned = true;
  }
  wrappedBreakpointModulator = newBreakpointModulatorToWrap;
  setModuleTypeName("BreakpointModulator");
  createParameters();
}

BreakpointModulatorAudioModule::~BreakpointModulatorAudioModule()
{
  if(wrappedBreakpointModulatorIsOwned)
    delete wrappedBreakpointModulator;
}

//-------------------------------------------------------------------------------------------------
// automation:

void BreakpointModulatorAudioModule::parameterChanged(Parameter* parameterThatHasChanged)
{
  /*
  // old parameter handling (obsolete, now using callbacks)
  if( wrappedBreakpointModulator == NULL )
    return;

  double value = parameterThatHasChanged->getValue();
  switch( getIndexOfParameter(parameterThatHasChanged) )
  {
  case  0: wrappedBreakpointModulator->setTimeScale(     value); break;
  case  1: wrappedBreakpointModulator->setTimeScaleByKey(value); break;
  case  2: wrappedBreakpointModulator->setTimeScaleByVel(value); break;
  case  3: wrappedBreakpointModulator->setDepth(         value); break;
  case  4: wrappedBreakpointModulator->setDepthByKey(    value); break;
  case  5: wrappedBreakpointModulator->setDepthByVel(    value); break;
  } // end of switch( parameterIndex )
  */

  //markStateAsDirty(); // leads to always marked as dirty even after loading preset
}

//-------------------------------------------------------------------------------------------------
// state saving and recall:

int stringToModBreakpointShapeIndex(const juce::String &shapeString)
{
  // use a std::map for this:
  typedef rosic::ModBreakpoint BP;
  if(      shapeString == "Stairstep" )  return BP::STAIRSTEP;
  else if( shapeString == "Linear" )     return BP::LINEAR;
  else if( shapeString == "Smooth" )     return BP::SMOOTH;
  else if( shapeString == "Analog" )     return BP::ANALOG;
  else if( shapeString == "AntiAnalog" ) return BP::GROWING;
  else if( shapeString == "Sigmoid" )    return BP::SIGMOID;
  else if( shapeString == "Spikey" )     return BP::SPIKEY;
  else if( shapeString == "Sine 1" )     return BP::SINE_1;
  else if( shapeString == "Sine 2" )     return BP::SINE_2;
  else                                   return BP::LINEAR;
}

const juce::String modBreakpointShapeIndexToString(int shapeIndex)
{
  typedef rosic::ModBreakpoint BP;
  if(      shapeIndex == BP::STAIRSTEP ) return "Stairstep";
  else if( shapeIndex == BP::LINEAR )    return "Linear";
  else if( shapeIndex == BP::SMOOTH )    return "Smooth";
  else if( shapeIndex == BP::ANALOG )    return "Analog";
  else if( shapeIndex == BP::GROWING )   return "AntiAnalog";
  else if( shapeIndex == BP::SIGMOID )   return "Sigmoid";
  else if( shapeIndex == BP::SPIKEY )    return "Spikey";
  else if( shapeIndex == BP::SINE_1 )    return "Sine 1";
  else if( shapeIndex == BP::SINE_2 )    return "Sine 2";
  else                                   return juce::String();
}

XmlElement* BreakpointModulatorAudioModule::getStateAsXml(const juce::String& stateName,
  bool markAsClean)
{
  // store the inherited controller mappings:
  XmlElement *xmlState = AudioModule::getStateAsXml(stateName, markAsClean);

  // add some attributes (todo: have parameter objects for this, so it can be done by baseclass):
  xmlState->setAttribute(juce::String("ScaleFactor"),      wrappedBreakpointModulator->getScaleFactor() );
  xmlState->setAttribute(juce::String("Offset"),           wrappedBreakpointModulator->getOffset() );

  // these 6 should be superfluous now:
  xmlState->setAttribute(juce::String("TimeScale"),        wrappedBreakpointModulator->getTimeScale() );
  xmlState->setAttribute(juce::String("TimeScaleByKey"),   wrappedBreakpointModulator->getTimeScaleByKey() );
  xmlState->setAttribute(juce::String("TimeScaleByVel"),   wrappedBreakpointModulator->getTimeScaleByVel() );
  xmlState->setAttribute(juce::String("Depth"),            wrappedBreakpointModulator->getDepth() );
  xmlState->setAttribute(juce::String("DepthByKey"),       wrappedBreakpointModulator->getDepthByKey() );
  xmlState->setAttribute(juce::String("DepthByVel"),       wrappedBreakpointModulator->getDepthByVel() );


  xmlState->setAttribute(juce::String("SyncMode"),         wrappedBreakpointModulator->isInSyncMode());
  xmlState->setAttribute(juce::String("LoopStartIndex"),   wrappedBreakpointModulator->getLoopStartIndex() );
  xmlState->setAttribute(juce::String("LoopEndIndex"),     wrappedBreakpointModulator->getLoopEndIndex() );

  if( wrappedBreakpointModulator->getLoopMode() == rosic::BreakpointModulator::FORWARD_LOOP )
    xmlState->setAttribute(juce::String("LoopMode"), juce::String("Forward") );
  else
    xmlState->setAttribute(juce::String("LoopMode"), juce::String("Off") );

  // create an XmlElement for each breakpoint and add it as child-XmlElement:
  for(int p = 0; p <= wrappedBreakpointModulator->lastBreakpointIndex(); p++)
  {
    XmlElement* breakpointState = new XmlElement(juce::String("Breakpoint"));

    breakpointState->setAttribute(juce::String("Time"),        wrappedBreakpointModulator->getBreakpointTime(p));
    breakpointState->setAttribute(juce::String("Level"),       wrappedBreakpointModulator->getBreakpointLevel(p));
    breakpointState->setAttribute(juce::String("Shape"),       modBreakpointShapeIndexToString(wrappedBreakpointModulator->getBreakpointShape(p) ) );
    breakpointState->setAttribute(juce::String("ShapeAmount"), wrappedBreakpointModulator->getBreakpointShapeAmount(p));
    xmlState->addChildElement(breakpointState);
  }

  return xmlState;
}

void BreakpointModulatorAudioModule::setStateFromXml(const XmlElement& xmlState,
  const juce::String& stateName, bool markAsClean)
{
  // restore the inherited controller mappings:
  AudioModule::setStateFromXml(xmlState, stateName, markAsClean);

  rosic::BreakpointModulator *modulator = wrappedBreakpointModulator; // use a shorter name here...

  // restore the settings:
  modulator->setScaleFactor(xmlState.getDoubleAttribute(   "ScaleFactor",    1.0));
  modulator->setOffset(xmlState.getDoubleAttribute(        "Offset",         0.0));
  modulator->setTimeScale(xmlState.getDoubleAttribute(     "TimeScale",      1.0));
  modulator->setTimeScaleByKey(xmlState.getDoubleAttribute("TimeScaleByKey", 0.0));
  modulator->setTimeScaleByVel(xmlState.getDoubleAttribute("TimeScaleByVel", 0.0));
  modulator->setDepth(xmlState.getDoubleAttribute(         "Depth",          1.0));
  modulator->setDepthByKey(xmlState.getDoubleAttribute(    "DepthByKey",     0.0));
  modulator->setDepthByVel(xmlState.getDoubleAttribute(    "DepthByVel",     0.0));

  juce::String stringValue; // for temporary storage

  // check if there are at least two child elements corresponding to two breakpoints (this is the
  // minimum allowed number):
  if( xmlState.getNumChildElements() < 2 )
    jassertfalse;

  // initialize our preview-BreakpointModulator (this will result in two breakpoints at times 0.0
  // and 1.0 with levels 0.0 both):
  modulator->initialize();

  double time, level, shapeAmount;
  int    shape = 0;
  // for temporary storage of retrieved values

  // retrieve all child elements with tag 'Breakpoint' and put them into an array:
  juce::Array<XmlElement*> breakpointElements;
  findChildElementsWithTagName(breakpointElements, xmlState, "Breakpoint");
  if( breakpointElements.size() < 2 )
    jassertfalse;

  // restore the data of the first breakpoint:
  XmlElement* breakpointState;
  breakpointState = breakpointElements[0];
  if( breakpointState->hasAttribute("Time")      &&
    breakpointState->hasAttribute(  "Level")     &&
    breakpointState->hasAttribute(  "Shape")     &&
    breakpointState->hasAttribute(  "ShapeAmount") )
  {
    level       = breakpointState->getDoubleAttribute("Level", 0.0);
    shapeAmount = breakpointState->getDoubleAttribute("ShapeAmount", 1.0);
    stringValue = breakpointState->getStringAttribute("Shape", juce::String());
    shape       = stringToModBreakpointShapeIndex(stringValue);
    modulator->modifyBreakpoint(0, 0.0, level, shape, shapeAmount);
  }
  else
    jassertfalse;

  // restore the data of the last breakpoint:
  breakpointState = breakpointElements[breakpointElements.size()-1];
  if( breakpointState->hasAttribute("Time")        &&
    breakpointState->hasAttribute(  "Level")       &&
    breakpointState->hasAttribute(  "Shape")       &&
    breakpointState->hasAttribute(  "ShapeAmount") )
  {
    time        = breakpointState->getDoubleAttribute("Time",        1.0);
    level       = breakpointState->getDoubleAttribute("Level",       0.0);
    shapeAmount = breakpointState->getDoubleAttribute("ShapeAmount", 1.0);
    stringValue = breakpointState->getStringAttribute("Shape", juce::String());
    shape       = stringToModBreakpointShapeIndex(stringValue);
    modulator->modifyBreakpoint(1, time, level, shape, shapeAmount);
  }
  else
    jassertfalse;

  // insert all the intermediate breakpoints with their data:
  int breakpointIndex = 1;
  for( breakpointIndex=1; breakpointIndex<=breakpointElements.size()-2; breakpointIndex++ )
  {
    breakpointState = breakpointElements[breakpointIndex];
    if( breakpointState->hasAttribute("Time")        &&
      breakpointState->hasAttribute(  "Level")       &&
      breakpointState->hasAttribute(  "Shape")       &&
      breakpointState->hasAttribute(  "ShapeAmount")     )
    {
      time        = breakpointState->getDoubleAttribute("Time", 1.0);
      level       = breakpointState->getDoubleAttribute("Level", 0.0);
      shapeAmount = breakpointState->getDoubleAttribute("ShapeAmount", 1.0);
      stringValue = breakpointState->getStringAttribute("Shape", juce::String());
      shape       = stringToModBreakpointShapeIndex(stringValue);
      modulator->insertBreakpoint(time, level, shape, shapeAmount);
    }
    else
      jassertfalse;
  }

  // restore the loop-settings:
  if( xmlState.hasAttribute("LoopMode") )
  {
    stringValue = xmlState.getStringAttribute("LoopMode" );
    if( stringValue == "Off" )
      modulator->setLoopMode(false);
    else if( stringValue == "Forward" )
      modulator->setLoopMode(true);
  }
  else
    jassertfalse;

  int index = 0;
  if( xmlState.hasAttribute("LoopStartIndex") )
  {
    index = xmlState.getIntAttribute("LoopStartIndex" );
    modulator->setLoopStartIndex(index);
  }
  else
    jassertfalse;

  if( xmlState.hasAttribute("LoopEndIndex") )
  {
    index = xmlState.getIntAttribute("LoopEndIndex" );
    modulator->setLoopEndIndex(index);
  }
  else
    jassertfalse;

  // restore the sync-setting:
  if( xmlState.hasAttribute("SyncMode") )
  {
    modulator->setSyncMode(xmlState.getBoolAttribute("SyncMode", true));
  }
  else
    jassertfalse;

  //if(markAsClean)
  //  markStateAsClean(); // does not work

  sendChangeMessage();
}

void BreakpointModulatorAudioModule::setStateToDefaults()
{
  if( wrappedBreakpointModulator != NULL )
    wrappedBreakpointModulator->initialize();
  sendChangeMessage();
}

AudioModuleEditor* BreakpointModulatorAudioModule::createEditor(int type)
{ 
  return new BreakpointModulatorEditor(lock, this);
  //return nullptr; 
}

//-------------------------------------------------------------------------------------------------
// internal functions:

void BreakpointModulatorAudioModule::createParameters()
{
  std::vector<double> defaultValues;

  typedef ModulatableParameter Param;
  Param* p;

  typedef rosic::BreakpointModulator BM;
  BM* bm = wrappedBreakpointModulator;

  p = new Param("TimeScale", 0.0625, 16.0, 1.0, Parameter::EXPONENTIAL);
  defaultValues.push_back(1.0/16.0);
  defaultValues.push_back(1.0/12.0);
  defaultValues.push_back(1.0/8.0);
  defaultValues.push_back(1.0/6.0);
  defaultValues.push_back(1.0/4.0);
  defaultValues.push_back(1.0/3.0);
  defaultValues.push_back(1.0/2.0);
  defaultValues.push_back(3.0/4.0);
  defaultValues.push_back(1.0);
  defaultValues.push_back(3.0/2.0);
  defaultValues.push_back(2.0);
  defaultValues.push_back(3.0);
  defaultValues.push_back(4.0);
  defaultValues.push_back(6.0);
  defaultValues.push_back(8.0);
  defaultValues.push_back(12.0);
  defaultValues.push_back(16.0);
  p->setDefaultValues(defaultValues);
  p->setValueChangeCallback<BM>(bm, &BM::setTimeScale);
  addObservedParameter(p);

  p = new Param("TimeScaleByKey", -150.0, 150.0, 0.0, Parameter::LINEAR, 1.0);
  defaultValues.clear();
  defaultValues.push_back(-100.0);
  defaultValues.push_back(0.0);
  defaultValues.push_back(100.0);
  p->setDefaultValues(defaultValues);
  p->setValueChangeCallback<BM>(bm, &BM::setTimeScaleByKey);
  addObservedParameter(p);

  p = new Param("TimeScaleByVel", -150.0, 150.0, 0.0, Parameter::LINEAR, 1.0);
  defaultValues.clear();
  defaultValues.push_back(-100.0);
  defaultValues.push_back(0.0);
  defaultValues.push_back(100.0);
  p->setDefaultValues(defaultValues);
  p->setValueChangeCallback<BM>(bm, &BM::setTimeScaleByVel);
  addObservedParameter(p);

  p = new Param("Depth", 0.0, 4.0, 1.0, Parameter::LINEAR);
  defaultValues.clear();
  defaultValues.push_back(1.0/16.0);
  defaultValues.push_back(1.0/12.0);
  defaultValues.push_back(1.0/8.0);
  defaultValues.push_back(1.0/6.0);
  defaultValues.push_back(1.0/4.0);
  defaultValues.push_back(1.0/3.0);
  defaultValues.push_back(1.0/2.0);
  defaultValues.push_back(3.0/4.0);
  defaultValues.push_back(1.0);
  defaultValues.push_back(3.0/2.0);
  defaultValues.push_back(2.0);
  defaultValues.push_back(3.0);
  defaultValues.push_back(4.0);
  defaultValues.push_back(6.0);
  defaultValues.push_back(8.0);
  defaultValues.push_back(12.0);
  defaultValues.push_back(16.0);
  p->setDefaultValues(defaultValues);
  p->setValueChangeCallback<BM>(bm, &BM::setDepth);
  addObservedParameter(p);

  p = new Param("DepthByKey", -150.0, 150.0, 0.0, Parameter::LINEAR, 1.0);
  defaultValues.clear();
  defaultValues.push_back(-100.0);
  defaultValues.push_back(0.0);
  defaultValues.push_back(100.0);
  p->setDefaultValues(defaultValues);
  p->setValueChangeCallback<BM>(bm, &BM::setDepthByKey);
  addObservedParameter(p);

  p = new Param("DepthByVel", -150.0, 150.0, 0.0, Parameter::LINEAR, 1.0);
  defaultValues.clear();
  defaultValues.push_back(-100.0);
  defaultValues.push_back(0.0);
  defaultValues.push_back(100.0);
  p->setDefaultValues(defaultValues);
  p->setValueChangeCallback<BM>(bm, &BM::setDepthByVel);
  addObservedParameter(p);

  // new parameters added 2017 for save/recall grid/snap settings:
  addObservedParameter(new ParameterGridInterval("GridX"));
  addObservedParameter(new ParameterGridInterval("GridY"));
  addObservedParameter(new Parameter("SnapX", 0.0, 1.0, 0.0, Parameter::BOOLEAN, 1.0));
  addObservedParameter(new Parameter("SnapY", 0.0, 1.0, 0.0, Parameter::BOOLEAN, 1.0));

  /*
  // make a call to setValue for each parameter in order to set up all the slave voices (still
  // necessary?):
  for(int i=0; i < (int) parameters.size(); i++ )
    parameterChanged(parameters[i]);
  */
}
