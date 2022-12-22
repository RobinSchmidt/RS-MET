

TrackMeterAudioModule::TrackMeterAudioModule(CriticalSection *newPlugInLock, 
  rosic::TrackMeter *trackMeterToWrap)
: AudioModule(newPlugInLock)
{
  //jassert(trackMeterToWrap != NULL); // you must pass a valid rosic-object to the constructor

  if( trackMeterToWrap != nullptr )
    wrappedTrackMeter = trackMeterToWrap;
  else
  {
    wrappedTrackMeter = new rosic::TrackMeter;
    wrappedTrackMeterIsOwned = true;
  }
  setModuleTypeName("TrackMeter");
  initializeAutomatableParameters();
}

TrackMeterAudioModule::~TrackMeterAudioModule()
{
  if(wrappedTrackMeterIsOwned)
    delete wrappedTrackMeter;
}

AudioModuleEditor* TrackMeterAudioModule::createEditor(int type)
{
  return new TrackMeterModuleEditor(lock, this);
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
  //setHeadlineText( juce::String(("TrackMeter")) );

  // assign the pointer to the rosic::TrackMeter object to be used as aduio engine:
  jassert(newTrackMeterAudioModule != NULL ); // you must pass a valid module here
  trackMeterModuleToEdit = newTrackMeterAudioModule;

  rangeMin = -60.f;
  rangeMax = +6.f;

  // create the widgets and assign the automatable parameters to them:
  addWidget( riseSlider = new RSlider("RiseSlider") );
  riseSlider->assignParameter( trackMeterModuleToEdit->getParameterByName("RiseTime") );
  riseSlider->setSliderName(juce::String("Rise"));
  riseSlider->setDescription(juce::String("Rise/Attack time-constant"));
  riseSlider->setDescriptionField(infoField);
  riseSlider->setStringConversionFunction(&millisecondsToStringWithUnit2);

  addWidget( fallSlider = new RSlider ("FallSlider") );
  fallSlider->assignParameter( trackMeterModuleToEdit->getParameterByName("FallTime") );
  fallSlider->setSliderName(juce::String("Fall"));
  fallSlider->setDescription(juce::String("Fall/Release time-constant"));
  fallSlider->setDescriptionField(infoField);
  fallSlider->setStringConversionFunction(&millisecondsToStringWithUnit2);

  addWidget( vuButton = new RButton("VU") );
  //vuButton->assignParameter( trackMeterModuleToEdit->getParameterByName(("VU")) );
  vuButton->setDescription(juce::String("Set ballistics to VU mode"));
  vuButton->setDescriptionField(infoField);
  vuButton->setClickingTogglesState(false);
  vuButton->addRButtonListener(this);

  addWidget( ppmButton = new RButton("PPM") );
  //ppmButton->assignParameter( trackMeterModuleToEdit->getParameterByName(("PPM")) );
  ppmButton->setDescription(juce::String("Set ballistics to PPM mode"));
  ppmButton->setDescriptionField(infoField);
  ppmButton->setClickingTogglesState(false);
  ppmButton->addRButtonListener(this);

  addWidget( leftLevelLabel = new RTextField(juce::String(("L"))) );
  leftLevelLabel->setDescription(juce::String(("Level of left channel")));
  leftLevelLabel->setDescriptionField(infoField);

  addWidget( leftLevelMeter = new MeteringDisplay() );
  leftLevelMeter->setDescription(leftLevelLabel->getDescription());
  leftLevelMeter->setDescriptionField(infoField);
  leftLevelMeter->setMeterStyle(MeteringDisplay::levelMeterStyle);
  leftLevelMeter->setRange(rangeMin, rangeMax);
  leftLevelMeter->setReferenceValue(0.0);
  leftLevelMeter->setCurrentValue(0.0);

  addWidget( rightLevelLabel = new RTextField(juce::String(("R"))) );
  rightLevelLabel->setDescription(juce::String(("Level of right channel")));
  rightLevelLabel->setDescriptionField(infoField);

  addWidget( rightLevelMeter = new MeteringDisplay() );
  rightLevelMeter->setDescription(rightLevelLabel->getDescription());
  rightLevelMeter->setDescriptionField(infoField);
  rightLevelMeter->setMeterStyle(MeteringDisplay::levelMeterStyle);
  rightLevelMeter->setRange(rangeMin, rangeMax);
  rightLevelMeter->setReferenceValue(0.0);
  rightLevelMeter->setCurrentValue(0.0);

  addWidget( midLevelLabel = new RTextField(juce::String(("M"))) );
  midLevelLabel->setDescription(juce::String(("Level of mid channel")));
  midLevelLabel->setDescriptionField(infoField);

  midLevelMeter = new MeteringDisplay();
  midLevelMeter->setDescription(midLevelLabel->getDescription());
  midLevelMeter->setDescriptionField(infoField);
  midLevelMeter->setMeterStyle(MeteringDisplay::levelMeterStyle);
  midLevelMeter->setRange(rangeMin, rangeMax);
  midLevelMeter->setReferenceValue(0.0);
  midLevelMeter->setCurrentValue(0.0);
  addWidget(midLevelMeter);

  sideLevelLabel = new RTextField(juce::String(("S")) );
  sideLevelLabel->setDescription(juce::String(("Level of side channel")));
  sideLevelLabel->setDescriptionField(infoField);
  addWidget(sideLevelLabel);

  sideLevelMeter = new MeteringDisplay();
  sideLevelMeter->setDescription(sideLevelLabel->getDescription());
  sideLevelMeter->setDescriptionField(infoField);
  sideLevelMeter->setMeterStyle(MeteringDisplay::levelMeterStyle);
  sideLevelMeter->setRange(rangeMin, rangeMax);
  sideLevelMeter->setReferenceValue(0.0);
  sideLevelMeter->setCurrentValue(0.0);
  addWidget(sideLevelMeter);

  correlationLabel = new RTextField(juce::String(("C")) );
  correlationLabel->setDescription(juce::String(("L/R cross-correlation")));
  correlationLabel->setDescriptionField(infoField);
  addWidget(correlationLabel);

  correlationMeter = new MeteringDisplay();
  correlationMeter->setDescription(correlationLabel->getDescription());
  correlationMeter->setDescriptionField(infoField);
  correlationMeter->setMeterStyle(MeteringDisplay::triangularPointerStyle);
  correlationMeter->setRange(-1.0, 1.0);
  correlationMeter->setReferenceValue(0.0);
  correlationMeter->setCurrentValue(0.0);
  addWidget(correlationMeter);

  // set up the widgets:
  updateWidgetsAccordingToState();

  //startTimer(20);  // 20 ms -> 50 FPS
  startTimerHz(60);  // 60 FPS

  setSize(180, 400); // changing this seems to have has no effect so it probably doesn't work
}

//-------------------------------------------------------------------------------------------------
// callbacks:

void TrackMeterModuleEditor::rButtonClicked(RButton* buttonThatWasClicked)
{
  double time, level, tauRise = 0, tauFall = 0;

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
    level   = 1.0-RAPT::rsDbToAmp(-24.0);
    tauFall = -time / log(1.0-level);
  }

  // convert time constants from seconds to milliseconds:
  tauRise *= 1000.0;
  tauFall *= 1000.0;

  // set the slider values (which also causes an update in the audio engine):
  riseSlider->setValue(tauRise, true);
  fallSlider->setValue(tauFall, true);

  // Rise = 1 ms, Fall = 100 ms seems to be quite nice for a fast responding meter
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

  const int t = 16;       // thickness of the meters
  const int m = 4;        // margin

  if(isVertical())
  {
    int x = 0;
    int y = infoField->getY()-2*(t+m);

    vuButton->setBounds(  x+m, y,     t*2,   t);
    ppmButton->setBounds( x+m, y+t+m, t*2,   t);

    x = vuButton->getRight();
    int w = getWidth()-x;
    riseSlider->setBounds(x+m, y,     w-m*2, t);
    fallSlider->setBounds(x+m, y+t+m, w-m*2, t);

    x = 40;
    y = getHeadlineBottom();
    int h = riseSlider->getY()-y;

    // Set up labels:
    leftLevelLabel->setBounds(  x, y+m,                           t+m, t    );
    leftLevelMeter->setBounds(  x, leftLevelLabel->getBottom(),   t,   h-t*2);
    x += t+m;
    rightLevelLabel->setBounds( x, y+m,                           t,   t    );
    rightLevelMeter->setBounds( x, rightLevelLabel->getBottom(),  t,   h-t*2);
    x += t*2;
    midLevelLabel->setBounds(   x, y+m,                           t,   t    );
    midLevelMeter->setBounds(   x, midLevelLabel->getBottom(),    t,   h-t*2); // maybe use y+t
    x += t+m;
    sideLevelLabel->setBounds(  x, y+m,                           t,   t    );
    sideLevelMeter->setBounds(  x, sideLevelLabel->getBottom(),   t,   h-t*2);
    x += t*2;
    correlationLabel->setBounds(x, y+m,                           t,   t    );
    correlationMeter->setBounds(x, correlationLabel->getBottom(), t,   h-t*2);
  }
  else
  {
    int x = m;    // start coordinate of the meters to the left
    int y = 2*t;  // not sure about that value. maybe use a hardcoded number - we'll see
    int w = getWidth() - x - 48; // 48: width of right section for the ballistics controls

    leftLevelLabel->setBounds(  x,   y, t, t);
    leftLevelMeter->setBounds(  x+t, y, w, t);
    y += t+m;
    rightLevelLabel->setBounds( x,   y, t, t);
    rightLevelMeter->setBounds( x+t, y, w, t);
    y += t*2;
    midLevelLabel->setBounds(   x,   y, t, t);
    midLevelMeter->setBounds(   x+t, y, w, t);
    y += t+m;
    sideLevelLabel->setBounds(  x,   y, t, t);
    sideLevelMeter->setBounds(  x+t, y, w, t);
    y += t*2;
    correlationLabel->setBounds(x,   y, t, t);
    correlationMeter->setBounds(x+t, y, w, t);
  }

  // ToDo:
  // -Maybe make the thickness variable depending on the available width (in vertical mode) or 
  //  height (in horizontal mode)
}

void TrackMeterModuleEditor::drawMeterScales(Graphics &g)
{
  double stepSize = 6.0;
  int    numSteps = roundToInt( (rangeMax-rangeMin) / stepSize ) + 1;

  const BitmapFontRoundedBoldA10D0* font = &BitmapFontRoundedBoldA10D0::instance;
  juce::String numberString;
  juce::Colour textColor = getTextColour();
  g.setColour(textColor);  // is this needed?
  // Color for the tick marks. Maybe darken for bright-on-dark schemes and brighten for 
  // dark-on-bright schemes maybe have a function in widgetColourScheme tonedDown that 
  // darkens or brightens and can be used.

  if(isVertical())
  {
    // Draw the scale for the level-meters:
    float x1  = (float)leftLevelMeter->getX()-8;
    float x2  = (float)sideLevelMeter->getRight();
    double y  = (float)leftLevelMeter->getY();
    double dy = stepSize * leftLevelMeter->getHeight() / (rangeMax - rangeMin);

    for(int i = 1; i <= numSteps; i++)
    {
      g.drawLine(x1, (float)y, x2, (float)y, 2.f);
      numberString = valueToStringWithSign0(rangeMax-(i-1)*stepSize);
      drawBitmapFontText(g, (int)x1 - 8, (int)y, numberString, font,
        textColor, -1, Justification::centredRight);
      y += dy;
    }

    // Draw the scale for the correlation-meter:
    x1       = (float)correlationMeter->getX();
    x2       = (float)correlationMeter->getRight()+8;
    y        = (float)correlationMeter->getY();
    stepSize = 0.2;
    dy       = stepSize * correlationMeter->getHeight() / 2.0;
    numSteps = roundToInt(2.0 / stepSize) + 1;
    for(int i = 1; i <= numSteps; i++)
    {
      g.drawLine(x1, (float)y, x2, (float)y, 2.f);
      numberString = valueToStringWithSign1(1.0-(i-1)*stepSize);
      //g.drawText(numberjuce::String, (int) x2, (int) (y-8), 30, 16, Justification::centredLeft, false);
      drawBitmapFontText(g, (int)x1 + 28, (int)y, numberString, font,
        textColor, -1, Justification::centredLeft);
      y += dy;
    }
  }
  else
  {
    //RAPT::rsError("not yet implemented");
    //g.fillAll(Colours::darkgrey);
    //drawBitmapFontText(g, 10, 10, "Horizontal mode not yet implemented", font,
    //  textColor, -1, Justification::centredLeft);

    // Draw the scale for the level-meters:
    float y1  = (float)leftLevelMeter->getY()-8;
    float y2  = (float)sideLevelMeter->getBottom();
    double x  = (float)leftLevelMeter->getX();
    double dx = stepSize * leftLevelMeter->getWidth() / (rangeMax - rangeMin);


    for(int i = 1; i <= numSteps; i++)
    {
      g.drawLine((float)x, y1, (float)x, y2, 2.f);
      //numberString = valueToStringWithSign0(rangeMax-(i-1)*stepSize);
      numberString = valueToStringWithSign0(rangeMin + (i-1)*stepSize);
      drawBitmapFontText(g, (int)x, (int)y1 - 12, numberString, font,
        textColor, -1, Justification::centredRight);
      x += dx;

    }

  }
}

void TrackMeterModuleEditor::timerCallback()
{
  if( trackMeterModuleToEdit == nullptr )
    return;
  if( trackMeterModuleToEdit->wrappedTrackMeter == nullptr )
    return;
  // Maybe wrap into a function isNull or arePointersAssigned()

  SignalMeasures measures = trackMeterModuleToEdit->wrappedTrackMeter->getCurrentMeasurement();

  leftLevelMeter->setCurrentValue(measures.leftLevel);
  rightLevelMeter->setCurrentValue(measures.rightLevel);
  midLevelMeter->setCurrentValue(measures.midLevel);
  sideLevelMeter->setCurrentValue(measures.sideLevel);
  correlationMeter->setCurrentValue(measures.crossCorrelation);
}

/*

Ideas:
-Let it optionally show the meters horizontally. Maybe it should automatically switch to horizontal
 mode when width > height. Rational: when recoding videos on a 1920x1080 resolution, we can put 
 Straightliner with 2 MultiAnalyzers by its (right) side - but an additional TackMeter fits 
 nowhere. A horizontal one, however, could fit below.
-Maybe don't put the ballistics adjustmenst onto the main GUI. Hide them behind a context menu.
-Maybe get rid of the headline "Slot X - TrackMeter". It's a bit pointless and silly

Bugs:
-when playing notes with the sampler, it hickups on note-off
-VU and PPM buttons don't flash on click
-PPM button too far down

*/