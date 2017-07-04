
RayBouncerAudioModule::RayBouncerAudioModule(CriticalSection *lockToUse) : AudioModule(lockToUse)
{
  ScopedLock scopedLock(*lock);
  moduleName = "RayBouncer";
  setActiveDirectory(getApplicationDirectory() + "Presets/RayBouncer");

  createParameters();
}

void RayBouncerAudioModule::createParameters()
{
  ScopedLock scopedLock(*lock);

  std::vector<double> defaultValues;
  typedef MetaControlledParameter Param;
  Param* p;

  p = new Param("Frequency", 20.0, 20000.0, 1000.0, Parameter::EXPONENTIAL);
  defaultValues.clear();
  defaultValues.push_back(125.0);
  defaultValues.push_back(250.0);
  defaultValues.push_back(500.0);
  defaultValues.push_back(1000.0);
  defaultValues.push_back(2000.0);
  defaultValues.push_back(4000.0);
  defaultValues.push_back(8000.0);
  defaultValues.push_back(16000.0);
  p->setDefaultValues(defaultValues);
  addObservedParameter(p);
  p->setValueChangeCallback<RayBouncerAudioModule>(this, &RayBouncerAudioModule::setFrequency);

  typedef RayBouncerAudioModule RBAM;

  p = new Param("Size", 0.25, 1.0, 1.0, Parameter::EXPONENTIAL);
  addObservedParameter(p);
  p->setValueChangeCallback<RBAM>(this, &RBAM::setEllipseSize);

  p = new Param("AspectRatio", 0.5, 2.0, 1.0, Parameter::EXPONENTIAL);
  addObservedParameter(p);
  p->setValueChangeCallback<RBAM>(this, &RBAM::setEllipseAspectRatio);

  p = new Param("Angle", 0.0, 360.0, 0.0, Parameter::LINEAR);
  addObservedParameter(p);
  p->setValueChangeCallback<RBAM>(this, &RBAM::setEllipseAngleDegrees);

  p = new Param("CenterX", -0.2, 0.2, 0.0, Parameter::LINEAR);
  addObservedParameter(p);
  p->setValueChangeCallback<RBAM>(this, &RBAM::setEllipseCenterX);

  p = new Param("CenterY", -0.2, 0.2, 0.0, Parameter::LINEAR);
  addObservedParameter(p);
  p->setValueChangeCallback<RBAM>(this, &RBAM::setEllipseCenterY);

  p = new Param("StartX", -0.2, 0.2, 0.0, Parameter::LINEAR);
  addObservedParameter(p);
  p->setValueChangeCallback<RBAM>(this, &RBAM::setStartX);

  p = new Param("StartY", -0.2, 0.2, 0.0, Parameter::LINEAR);
  addObservedParameter(p);
  p->setValueChangeCallback<RBAM>(this, &RBAM::setStartY);

  p = new Param("LaunchAngle", 0.0, 360.0, 45.0, Parameter::LINEAR);
  addObservedParameter(p);
  p->setValueChangeCallback<RBAM>(this, &RBAM::setLaunchAngle);

  p = new Param("Bending", -1.0, 1.0, 0.0, Parameter::LINEAR);
  addObservedParameter(p);
  p->setValueChangeCallback<RBAM>(this, &RBAM::setBending);


  p = new Param("BendX2Y", -0.01, 0.01, 0.0, Parameter::LINEAR);
  addObservedParameter(p);
  p->setValueChangeCallback<RBAM>(this, &RBAM::setBendX2Y);

  p = new Param("BendY2X", -0.01, 0.01, 0.0, Parameter::LINEAR);
  addObservedParameter(p);
  p->setValueChangeCallback<RBAM>(this, &RBAM::setBendY2X);


  p = new Param("BendXX2X", -2.0, 2.0, 0.0, Parameter::LINEAR);
  addObservedParameter(p);
  p->setValueChangeCallback<RBAM>(this, &RBAM::setBendXX2X);

  p = new Param("BendXX2Y", -2.0, 2.0, 0.0, Parameter::LINEAR);
  addObservedParameter(p);
  p->setValueChangeCallback<RBAM>(this, &RBAM::setBendXX2Y);


  p = new Param("BendXY2X", -2.0, 2.0, 0.0, Parameter::LINEAR);
  addObservedParameter(p);
  p->setValueChangeCallback<RBAM>(this, &RBAM::setBendXY2X);

  p = new Param("BendXY2Y", -2.0, 2.0, 0.0, Parameter::LINEAR);
  addObservedParameter(p);
  p->setValueChangeCallback<RBAM>(this, &RBAM::setBendXY2Y);


  p = new Param("BendYY2X", -2.0, 2.0, 0.0, Parameter::LINEAR);
  addObservedParameter(p);
  p->setValueChangeCallback<RBAM>(this, &RBAM::setBendYY2X);

  p = new Param("BendYY2Y", -2.0, 2.0, 0.0, Parameter::LINEAR);
  addObservedParameter(p);
  p->setValueChangeCallback<RBAM>(this, &RBAM::setBendYY2Y);

  // maybe eventually we only want to retain xy2x, xy2y and let them be equal for symmetry
  // but maybe x2y, and y2x make also sense



  p = new Param("AutoReset", 0.0, 1.0, 0.0, Parameter::BOOLEAN);
  addObservedParameter(p);
  p->setValueChangeCallback<RBAM>(this, &RBAM::setAutoReset);

  // maybe we need to call the valueChangeCallback in setValueChangeCallback to make everything
  // initially consistent...

  // add initial condition parameters - the output is very sensitive to the conditions
  // swicthing presets back and forth always sounds different because of different initial 
  // conditions
}

// Editor creation:

//AudioModuleEditor* RayBouncerAudioModule::createEditor()
//{
//  return new jura::RayBouncerAudioModuleEditor(this);
//}

// audio processing:

void RayBouncerAudioModule::processBlock(double **inOutBuffer, int numChannels, int numSamples)
{
  jassert(numChannels == 2);
  double* x = inOutBuffer[0];
  double* y = inOutBuffer[1];
  for(int n = 0; n < numSamples; n++)
    rayBouncer.getSampleFrame(x[n], y[n]);
}

void RayBouncerAudioModule::setSampleRate(double newSampleRate)
{
  sampleRate = newSampleRate;
  rayBouncer.setFrequencyAndSampleRate(frequency, sampleRate);
}

void RayBouncerAudioModule::reset()
{
  rayBouncer.reset();
}

void RayBouncerAudioModule::setStateFromXml(const XmlElement& xmlState, 
  const juce::String& stateName, bool markAsClean)
{
  AudioModule::setStateFromXml(xmlState, stateName, markAsClean);
  rayBouncer.reset();
}

// parameter setters (callback targets for the Parameter objects):

void RayBouncerAudioModule::setFrequency(double newFrequency)
{
  frequency = newFrequency;
  rayBouncer.setFrequencyAndSampleRate(frequency, sampleRate);
}
void RayBouncerAudioModule::setEllipseSize(double newSize)
{
  rayBouncer.setEllipseSize(newSize);
  autoResetIfDesired();
}
void RayBouncerAudioModule::setEllipseAspectRatio(double newRatio)
{
  rayBouncer.setEllipseAspectRatio(newRatio);
  autoResetIfDesired();
}
void RayBouncerAudioModule::setEllipseAngleDegrees(double newAngle)
{
  rayBouncer.setEllipseAngleDegrees(newAngle);
  autoResetIfDesired();
}
void RayBouncerAudioModule::setEllipseCenterX(double newX)
{
  rayBouncer.setEllipseCenterX(newX);
  autoResetIfDesired();
}
void RayBouncerAudioModule::setEllipseCenterY(double newY)
{
  rayBouncer.setEllipseCenterY(newY);
  autoResetIfDesired();
}
void RayBouncerAudioModule::setStartX(double newX)
{
  rayBouncer.setStartX(newX);
  autoResetIfDesired();
}
void RayBouncerAudioModule::setStartY(double newY)
{
  rayBouncer.setStartY(newY);
  autoResetIfDesired();
}
void RayBouncerAudioModule::setLaunchAngle(double newAngle)
{
  rayBouncer.setLaunchAngleDegrees(newAngle);
  autoResetIfDesired();
}
void RayBouncerAudioModule::setBending(double newValue)
{
  rayBouncer.rayBouncer.setBend(newValue);
  autoResetIfDesired();
}

void RayBouncerAudioModule::setBendX2Y(double newValue)
{
  rayBouncer.rayBouncer.setBendX2Y(newValue);
  autoResetIfDesired();
}
void RayBouncerAudioModule::setBendY2X(double newValue)
{
  rayBouncer.rayBouncer.setBendY2X(newValue);
  autoResetIfDesired();
}

void RayBouncerAudioModule::setBendXX2X(double newValue)
{
  rayBouncer.rayBouncer.setBendXX2X(newValue);
  autoResetIfDesired();
}
void RayBouncerAudioModule::setBendXX2Y(double newValue)
{
  rayBouncer.rayBouncer.setBendXX2Y(newValue);
  autoResetIfDesired();
}

void RayBouncerAudioModule::setBendXY2X(double newValue)
{
  rayBouncer.rayBouncer.setBendXY2X(newValue);
  autoResetIfDesired();
}
void RayBouncerAudioModule::setBendXY2Y(double newValue)
{
  rayBouncer.rayBouncer.setBendXY2Y(newValue);
  autoResetIfDesired();
}

void RayBouncerAudioModule::setBendYY2X(double newValue)
{
  rayBouncer.rayBouncer.setBendYY2X(newValue);
  autoResetIfDesired();
}
void RayBouncerAudioModule::setBendYY2Y(double newValue)
{
  rayBouncer.rayBouncer.setBendYY2Y(newValue);
  autoResetIfDesired();
}

//void setBendXX2Y(double newValue);
//void setBendXY2X(double newValue);
//void setBendXY2Y(double newValue);
//void setBendYY2X(double newValue);
//void setBendYY2Y(double newValue);


void RayBouncerAudioModule::setAutoReset(bool shouldReset)
{
  autoReset = shouldReset;
}

// misc:

void RayBouncerAudioModule::autoResetIfDesired()
{
  if(autoReset)
    rayBouncer.reset();
}

//=================================================================================================
// the GUI editor class for the RayBouncer:

//RayBouncerAudioModuleEditor::RayBouncerAudioModuleEditor(
//  jura::RayBouncerAudioModule *newRayBouncerToEdit) 
//  : AudioModuleEditor(newRayBouncerToEdit)
//{
//  ScopedLock scopedLock(*lock);
//
//  setHeadlineText("RayBouncer");
//
//  // assign the pointer to the edited object:
//  jassert(newRayBouncerToEdit != nullptr ); // you must pass a valid module here
//  rayBouncerToEdit = newRayBouncerToEdit;
//
//  // create the widgets and assign the automatable parameters to them:
//  addWidget( cutoffSlider = new AutomatableSlider() );
//  cutoffSlider->assignParameter( ladderToEdit->getParameterByName("Cutoff") );
//  cutoffSlider->setSliderName("Cutoff");
//  cutoffSlider->setDescription("Cutoff frequency in Hz");
//  cutoffSlider->setDescriptionField(infoField);
//  cutoffSlider->setStringConversionFunction(&hertzToStringWithUnitTotal5);
//}

//void LadderEditor::resized()
//{
//  ScopedLock scopedLock(*lock);
//  AudioModuleEditor::resized();
//
//  int x  = 0;
//  int y  = 0;
//  int w  = getWidth();
//  int h  = getHeight();
//  y = getPresetSectionBottom();
//
//  cutoffSlider->setBounds(      4, y+4, w2-4, 16);
//}
