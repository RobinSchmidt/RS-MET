

TrackMeterAudioModule::TrackMeterAudioModule(CriticalSection *newPlugInLock, rosic::TrackMeter *trackMeterToWrap)
: AudioModule(newPlugInLock)
{
  jassert(trackMeterToWrap != NULL); // you must pass a valid rosic-object to the constructor
  wrappedTrackMeter = trackMeterToWrap;
  moduleName = juce::String("TrackMeter");
  initializeAutomatableParameters();
}

void TrackMeterAudioModule::parameterChanged(Parameter* parameterThatHasChanged)
{
  if( wrappedTrackMeter == NULL )
    return;

  double value = parameterThatHasChanged->getValue();
  switch( getIndexOfParameter(parameterThatHasChanged) )
  {
  case   0: wrappedTrackMeter->setAttackTimeInMilliseconds(value);   break;
  case   1: wrappedTrackMeter->setReleaseTimeInMilliseconds(value);  break;
  } // end of switch( parameterIndex )
  markStateAsDirty();
}

void TrackMeterAudioModule::initializeAutomatableParameters()
{
  // create the automatable parameters and add them to the list - note that the order of the adds
  // is important because in parameterChanged(), the index (position in the array) will be used to
  // identify which particlua parameter has changed.

  AutomatableParameter* p;
    
  // calculate the rise- and fall- time constants for VU-Meter to be used as default values:
  double time  = 0.3;  // 0.3 sec rise and fall time for VU-meters
  double level = 0.99; // after 0.3 sec, the level should have dropped to 1% (for fall) or 
                       // reached 99% for rise  
  double tau   = -1000.0*time / log(1.0-level);

  // #000:
  p = new AutomatableParameter(lock, "RiseTime", 1.0, 300.0, 0.0, tau, Parameter::EXPONENTIAL);
  addObservedParameter(p);

  // #001:
  p = new AutomatableParameter(lock, "FallTime", 1.0, 3000.0, 0.0, tau, Parameter::EXPONENTIAL);
  addObservedParameter(p);

  // make a call to parameterChanged for each parameter in order to set up the DSP-core to reflect 
  // the values the automatable parameters:
  for(int i=0; i < (int) parameters.size(); i++ )
    parameterChanged(parameters[i]);
}

//=================================================================================================

// construction/destruction:

TrackMeterModuleEditor::TrackMeterModuleEditor(CriticalSection *newPlugInLock, 
  TrackMeterAudioModule* newTrackMeterAudioModule) : AudioModuleEditor(newTrackMeterAudioModule)
{
  // set the plugIn-headline:
  setHeadlineText( juce::String(("TrackMeter")) );

  // assign the pointer to the rosic::TrackMeter object to be used as aduio engine:
  jassert(newTrackMeterAudioModule != NULL ); // you must pass a valid module here
  trackMeterModuleToEdit = newTrackMeterAudioModule;

  rangeMin = -60.0;
  rangeMax = +6.0;

  // create the widgets and assign the automatable parameters to them:
  addWidget( riseSlider = new RSlider (("RiseSlider")) );
  riseSlider->assignParameter( trackMeterModuleToEdit->getParameterByName(("RiseTime")) );
  riseSlider->setSliderName(juce::String(("Rise")));
  riseSlider->setDescription(juce::String(("Rise/Attack time-constant")));
  riseSlider->setDescriptionField(infoField);
  riseSlider->setStringConversionFunction(&millisecondsToStringWithUnit2);

  addWidget( fallSlider = new RSlider (("FallSlider")) );
  fallSlider->assignParameter( trackMeterModuleToEdit->getParameterByName(("FallTime")) );
  fallSlider->setSliderName(juce::String(("Fall")));
  fallSlider->setDescription(juce::String(("Fall/Release time-constant")));
  fallSlider->setDescriptionField(infoField);
  fallSlider->setStringConversionFunction(&millisecondsToStringWithUnit2);

  addWidget( vuButton = new RButton(juce::String(("VU"))) );
  vuButton->assignParameter( trackMeterModuleToEdit->getParameterByName(("VU")) );
  vuButton->setDescription(juce::String(("Set ballistics to VU mode")));
  vuButton->setDescriptionField(infoField);
  vuButton->setClickingTogglesState(false);
  vuButton->addRButtonListener(this);

  addWidget( ppmButton = new RButton(juce::String(("PPM"))) );
  ppmButton->assignParameter( trackMeterModuleToEdit->getParameterByName(("PPM")) );
  ppmButton->setDescription(juce::String(("Set ballistics to PPM mode")));
  ppmButton->setDescriptionField(infoField);
  ppmButton->setClickingTogglesState(false);
  ppmButton->addRButtonListener(this);

  addWidget( leftLevelLabel = new RTextField(juce::String(("L"))) );
  leftLevelLabel->setDescription(juce::String(("Level of left channel")));
  leftLevelLabel->setDescriptionField(infoField);

  addWidget( leftLevelMeter = new MeteringDisplay(juce::String(("L"))) );
  leftLevelMeter->setDescription(leftLevelLabel->getDescription());
  leftLevelMeter->setDescriptionField(infoField);
  leftLevelMeter->setMeterStyle(MeteringDisplay::levelMeterStyle);
  leftLevelMeter->setRange(rangeMin, rangeMax);
  leftLevelMeter->setReferenceValue(0.0);
  leftLevelMeter->setCurrentValue(0.0);

  addWidget( rightLevelLabel = new RTextField(juce::String(("R"))) );
  rightLevelLabel->setDescription(juce::String(("Level of right channel")));
  rightLevelLabel->setDescriptionField(infoField);

  addWidget( rightLevelMeter = new MeteringDisplay(juce::String(("R"))) );
  rightLevelMeter->setDescription(rightLevelLabel->getDescription());
  rightLevelMeter->setDescriptionField(infoField);
  rightLevelMeter->setMeterStyle(MeteringDisplay::levelMeterStyle);
  rightLevelMeter->setRange(rangeMin, rangeMax);
  rightLevelMeter->setReferenceValue(0.0);
  rightLevelMeter->setCurrentValue(0.0);

  addWidget( midLevelLabel = new RTextField(juce::String(("M"))) );
  midLevelLabel->setDescription(juce::String(("Level of mid channel")));
  midLevelLabel->setDescriptionField(infoField);

  midLevelMeter = new MeteringDisplay(juce::String(("M")) );
  midLevelMeter->setDescription(midLevelLabel->getDescription());
  midLevelMeter->setDescriptionField(infoField);
  midLevelMeter->setMeterStyle(MeteringDisplay::levelMeterStyle);
  midLevelMeter->setRange(rangeMin, rangeMax);
  midLevelMeter->setReferenceValue(0.0);
  midLevelMeter->setCurrentValue(0.0);
  addAndMakeVisible(midLevelMeter);

  sideLevelLabel = new RTextField(juce::String(("S")) );
  sideLevelLabel->setDescription(juce::String(("Level of side channel")));
  sideLevelLabel->setDescriptionField(infoField);
  addAndMakeVisible(sideLevelLabel);

  sideLevelMeter = new MeteringDisplay(juce::String(("S")) );
  sideLevelMeter->setDescription(sideLevelLabel->getDescription());
  sideLevelMeter->setDescriptionField(infoField);
  sideLevelMeter->setMeterStyle(MeteringDisplay::levelMeterStyle);
  sideLevelMeter->setRange(rangeMin, rangeMax);
  sideLevelMeter->setReferenceValue(0.0);
  sideLevelMeter->setCurrentValue(0.0);
  addAndMakeVisible(sideLevelMeter);

  correlationLabel = new RTextField(juce::String(("C")) );
  correlationLabel->setDescription(juce::String(("L/R cross-correlation")));
  correlationLabel->setDescriptionField(infoField);
  addAndMakeVisible(correlationLabel);

  correlationMeter = new MeteringDisplay(juce::String(("C")) );
  correlationMeter->setDescription(correlationLabel->getDescription());
  correlationMeter->setDescriptionField(infoField);
  correlationMeter->setMeterStyle(MeteringDisplay::triangularPointerStyle);
  correlationMeter->setRange(-1.0, 1.0);
  correlationMeter->setReferenceValue(0.0);
  correlationMeter->setCurrentValue(0.0);
  addAndMakeVisible(correlationMeter);

  // set up the widgets:
  updateWidgetsAccordingToState();

  startTimer(20);
}

//-------------------------------------------------------------------------------------------------
// callbacks:

void TrackMeterModuleEditor::rButtonClicked(RButton* buttonThatWasClicked)
{
  double time, level, tauRise, tauFall;

  if( buttonThatWasClicked == vuButton )
  {
    // calculate the rise- and fall- time constants:
    time  = 0.3;    // 0.3 sec rise and fall time for VU-meters
    level = 0.99;   // after 0.3 sec, the level should have dropped to 1% (for fall) or reached 
                    // 99% for rise
    tauRise = tauFall = -time / log(1.0-level);
  }
  else if( buttonThatWasClicked == ppmButton )
  {
    // calculate the rise-time:
    time    = 0.01;    
    level   = 0.8;    
    tauRise = -time / log(1.0-level);

    // calculate the fall-time:
    time    = 2.8;    
    level   = 1.0-dB2amp(-24.0);    
    tauFall = -time / log(1.0-level);
  }

  // convert time constants from seconds to milliseconds:
  tauRise *= 1000.0;
  tauFall *= 1000.0;

  // set the slider values (which also causes an update in the audio engine):
  riseSlider->setValue(tauRise, true, false);
  fallSlider->setValue(tauFall, true, false);
}

void TrackMeterModuleEditor::paint(Graphics &g)
{
  AudioModuleEditor::paint(g);
  drawMeterScales(g);
}

void TrackMeterModuleEditor::resized()
{
  presetSectionPosition = INVISIBLE;
  linkPosition          = RIGHT_TO_HEADLINE;
  AudioModuleEditor::resized();
  int x = 0;
  int y = 0;
  int w = getWidth();
  int h = getHeight();


  y = infoField->getY()-36;

  vuButton->setBounds(x+4,  y,    32, 16);
  ppmButton->setBounds(x+4, y+20, 32, 16);

  x = vuButton->getRight();
  w = getWidth()-x;
  riseSlider->setBounds(x+4, y,    w-8, 16);
  fallSlider->setBounds(x+4, y+20, w-8, 16);

  y = getHeadlineBottom();
  h = riseSlider->getY()-y;

  leftLevelLabel->setBounds(40, y+4, 20, 16);
  leftLevelMeter->setBounds(40, leftLevelLabel->getBottom(), 16, h-32);
  rightLevelLabel->setBounds(64, y+4, 16, 16);
  rightLevelMeter->setBounds(64, rightLevelLabel->getBottom(), 16, h-32);

  midLevelLabel->setBounds(96, y+4, 16, 16);
  midLevelMeter->setBounds(96, midLevelLabel->getBottom(), 16, h-32);
  sideLevelLabel->setBounds(120, y+4, 16, 16);
  sideLevelMeter->setBounds(120, sideLevelLabel->getBottom(), 16, h-32);

  correlationLabel->setBounds(152, y+4, 16, 16);
  correlationMeter->setBounds(152, correlationLabel->getBottom(), 16, h-32);

  y = correlationMeter->getBottom()+16;
}

void TrackMeterModuleEditor::drawMeterScales(Graphics &g)
{
  juce::String numberString;
  g.setColour(Colours::black);
  int i;

  const BitmapFontRoundedBoldA10D0* font = &BitmapFontRoundedBoldA10D0::instance;

  // draw the scale for the level-meters:
  float  x1       = (float) leftLevelMeter->getX()-8;
  float  x2       = (float) sideLevelMeter->getRight();
  double y        = (float) leftLevelMeter->getY();
  double stepSize = 6.0;
  double yStep    = stepSize * leftLevelMeter->getHeight() / (rangeMax - rangeMin);
  int    numSteps = roundToInt( (rangeMax-rangeMin) / stepSize ) + 1;
  for(i=1; i<=numSteps; i++)
  {
    g.drawLine(x1, (float) y, x2, (float) y, 2.f); 
    numberString = valueToStringWithSign0(rangeMax-(i-1)*stepSize);
    //g.drawText(numberjuce::String, (int) (x1-24), (int) (y-8), 24, 16, Justification::centredRight, false);
    drawBitmapFontText(g, (int) x1 - 8, (int) y, numberString, font, 
      Colours::black, -1, Justification::centredRight);
    y += yStep;
  }

  // draw the scale for the correlation-meter:
  x1       = (float) correlationMeter->getX();
  x2       = (float) correlationMeter->getRight()+8;
  y        = (float) correlationMeter->getY();
  stepSize = 0.2;
  yStep    = stepSize * correlationMeter->getHeight() / 2.0;
  numSteps = roundToInt( 2.0 / stepSize ) + 1;
  for(i=1; i<=numSteps; i++)
  {
    g.drawLine(x1, (float) y, x2, (float) y, 2.f);
    numberString = valueToStringWithSign1(1.0-(i-1)*stepSize);
    //g.drawText(numberjuce::String, (int) x2, (int) (y-8), 30, 16, Justification::centredLeft, false);
    drawBitmapFontText(g, (int) x1 + 28, (int) y, numberString, font, 
      Colours::black, -1, Justification::centredLeft);
    y += yStep;
  }
}

void TrackMeterModuleEditor::timerCallback()
{
  if( trackMeterModuleToEdit == NULL )
    return;
  if( trackMeterModuleToEdit->wrappedTrackMeter == NULL )
    return;

  SignalMeasures measures = trackMeterModuleToEdit->wrappedTrackMeter->getCurrentMeasurement();

  leftLevelMeter->setCurrentValue(measures.leftLevel);
  rightLevelMeter->setCurrentValue(measures.rightLevel);
  midLevelMeter->setCurrentValue(measures.midLevel);
  sideLevelMeter->setCurrentValue(measures.sideLevel);
  correlationMeter->setCurrentValue(measures.crossCorrelation);
}
