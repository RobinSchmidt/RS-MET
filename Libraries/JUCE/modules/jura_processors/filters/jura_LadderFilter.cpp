Ladder::Ladder(CriticalSection *lockToUse, MetaParameterManager* metaManagerToUse,
  ModulationManager* modManagerToUse) 
  : ModulatableAudioModule(lockToUse, metaManagerToUse, modManagerToUse)
{
  ScopedLock scopedLock(*lock);
  setModuleTypeName("Ladder");
  createStaticParameters();
}

void Ladder::createStaticParameters()
{
  ScopedLock scopedLock(*lock);

  std::vector<double> defaultValues;
  typedef ModulatableParameter Param;
  Param* p;

  p = new Param("Cutoff", 20.0, 20000.0, 1000.0, Parameter::EXPONENTIAL);
  p->setDefaultModParameters(20, 20000, -2, 3, true);
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
    // maybe we should have a subclass ModulatableAudioModule of AudioModule and a function
    // addModulatableParameter - this function could register the ModulationTarget for the
    // passed ModulatableParameter
  p->setValueChangeCallback<Ladder>(this, &Ladder::setCutoff);

  p = new Param("Resonance", 0.0, 1.0, 0.2, Parameter::LINEAR, 0.01);
  addObservedParameter(p);
  p->setValueChangeCallback<Ladder>(this, &Ladder::setResonance);

  p = new Param("StereoSpread", -24.0, +24.0, 0.0, Parameter::LINEAR_BIPOLAR);
  //p->setMapper(new rsParameterMapperSinh(-24, +24, 2.0)); // for test
  addObservedParameter(p);
  p->setValueChangeCallback<Ladder>(this, &Ladder::setStereoSpread);

  p = new Param("Mode", 0.0, 14.0, 4.0, Parameter::STRING);
  p->addStringValue("Flat");
  p->addStringValue("Lowpass 6 dB/oct");
  p->addStringValue("Lowpass 12 dB/oct");
  p->addStringValue("Lowpass 18 dB/oct");
  p->addStringValue("Lowpass 24 dB/oct");
  p->addStringValue("Highpass 6 dB/oct");
  p->addStringValue("Highpass 12 dB/oct");
  p->addStringValue("Highpass 18 dB/oct");
  p->addStringValue("Highpass 24 dB/oct");
  p->addStringValue("Bandpass 6/6 dB/oct");
  p->addStringValue("Bandpass 6/12 dB/oct");
  p->addStringValue("Bandpass 6/18 dB/oct");
  p->addStringValue("Bandpass 12/6 dB/oct");
  p->addStringValue("Bandpass 12/12 dB/oct");
  p->addStringValue("Bandpass 18/6 dB/oct");
  addObservedParameter(p);
  p->setValueChangeCallback<Ladder>(this, &Ladder::setMode);
}

//-------------------------------------------------------------------------------------------------
// Editor creation:

AudioModuleEditor* Ladder::createEditor()
{
  return new jura::LadderEditor(this);
}

//-------------------------------------------------------------------------------------------------
// the audio processing callbacks:

void Ladder::processBlock(double **inOutBuffer, int numChannels, int numSamples)
{
  jassert(numChannels == 2);
  for(int n = 0; n < numSamples; n++)
    wrappedLadder.getSampleFrameStereo(&inOutBuffer[0][n], &inOutBuffer[1][n]);
}

void Ladder::processStereoFrame(double *left, double *right)
{
  wrappedLadder.getSampleFrameStereo(left, right);
}

//-------------------------------------------------------------------------------------------------
// parameter setters (callback targets for the Parameter objects):

void Ladder::setSampleRate(double newSampleRate){ wrappedLadder.setSampleRate(newSampleRate); }
void Ladder::reset(){ wrappedLadder.reset(); }

// get rid of these:
void Ladder::setCutoff(double newCutoff)        { cutoff = newCutoff; wrappedLadder.setCutoff(cutoff); }
void Ladder::setResonance(double newResonance)  { wrappedLadder.setResonance(newResonance); }
void Ladder::setMode(int newMode)               { wrappedLadder.setMode(newMode); }


inline double pitchOffsetToFreqFactor(double pitchOffset){ // move to RAPT library, eventually
  return exp(0.057762265046662109118102676788181 * pitchOffset);
}
void Ladder::setStereoSpread(double newSpread){
  freqFactorL = pitchOffsetToFreqFactor(+0.5*newSpread);
  freqFactorR = pitchOffsetToFreqFactor(-0.5*newSpread);  // == 1 / freqFactorL -> optimize later
  setCutoff(cutoff);
  // By spreading the cutoff frequencies of left and right channel, we not only get a stereo 
  // effect, but also a formant-filter like sound. Maybe the effectiveness can be improved, if we
  // have another parameter that makes this spreading dependent on the cutoff frequency. This way, 
  // we could control distance of the formants as function of the nominal center frequency between 
  // them. I think that should allow for more realistic formant sweeps. The spread parameter is 
  // also an interesting target for an LFO.
}
void Ladder::setMidSideMode(bool shouldBeInMidSideMode){
  midSideMode = shouldBeInMidSideMode;
}

// inquiry:

double Ladder::getDecibelsAt(double frequency)  // rename to getDecibelsAt
{
  double tmp = wrappedLadder.getMagnitudeResponseAt(frequency);
  tmp = amp2dBWithCheck(tmp, 0.000001);
  tmp = clip(tmp, -120.0, +120.0);
  return tmp;
}

//=================================================================================================
// rsLadderPlotEditor:

rsLadderPlotEditor::rsLadderPlotEditor(jura::Ladder* ladder) : ladderToEdit(ladder)
{
  freqRespPlot = new rsFunctionPlot;
  //freqRespPlot->addMouseListener(this, true);
  //freqRespPlot->setupForDecibelsAgainstLogFrequency(15.625, 32000.0, -60.0, 60.0, 12);
  freqRespPlot->setupForDecibelsAgainstLogFrequency(20.0, 20000.0, -60.0, 60.0, 12); 
    // frequency range must match cutoff parameter range, otherwise the dot-handle and resonance 
    // freq are out of sync
  freqRespPlot->addFunction([this](double f)->double { return ladderToEdit->getDecibelsAt(f); } );    // maybe try to use a member-function pointer without lambda
  addPlot(freqRespPlot);

  cutoffParam = ladderToEdit->getParameterByName("Cutoff");
  cutoffParam->registerParameterObserver(this);

  resoParam = ladderToEdit->getParameterByName("Resonance");
  resoParam->registerParameterObserver(this);

  vectorPad = new rsVectorPad;
  vectorPad->assignParameterX(cutoffParam);
  vectorPad->assignParameterY(resoParam);
  vectorPad->setPaintBackground(false);
  vectorPad->setDotSize(8.f);
  //vectorPad->setMargins(20, 40, 20, 40); // preliminary
  addWidget(vectorPad);
}

rsLadderPlotEditor::~rsLadderPlotEditor()
{
  cutoffParam->deRegisterParameterObserver(this);
  resoParam->deRegisterParameterObserver(this);
  delete freqRespPlot;
  delete vectorPad;
}

void rsLadderPlotEditor::parameterChanged(Parameter* p)
{
  double resoFreq = cutoffParam->getValue();
  double reso     = resoParam->getValue();
  double peakFreq = resoFreq * pow(reso, 0.25);  // formula found by trial and error
  peakFreq = jmax(peakFreq, 20.0);               // avoid NaN in drawing code
  freqRespPlot->setSpecialEvaluationPoint(0, 0, peakFreq);
}

void rsLadderPlotEditor::resized()
{
  freqRespPlot->setBounds(0, 0, getWidth(), getHeight());
  vectorPad->setBounds(   0, 0, getWidth(), getHeight());
  vectorPad->adjustMarginsToPlotX(freqRespPlot);
}

//=================================================================================================
// LadderEditor:

LadderEditor::LadderEditor(jura::Ladder *newLadderToEdit) : AudioModuleEditor(newLadderToEdit)
{
  ScopedLock scopedLock(*lock);

  // assign the pointer to the edited object:
  jassert(newLadderToEdit != nullptr ); // you must pass a valid module here
  ladderToEdit = newLadderToEdit;

  plotEditor = new rsLadderPlotEditor(ladderToEdit);
  addChildColourSchemeComponent(plotEditor);

  addWidget( cutoffSlider = new ModulatableSlider() );
  cutoffSlider->assignParameter( ladderToEdit->getParameterByName("Cutoff") );
  cutoffSlider->setSliderName("Cutoff");
  cutoffSlider->setDescription("Cutoff frequency in Hz");
  cutoffSlider->setDescriptionField(infoField);
  cutoffSlider->setStringConversionFunction(&hertzToStringWithUnitTotal5);

  addWidget( resonanceSlider = new ModulatableSlider() );
  resonanceSlider->assignParameter( ladderToEdit->getParameterByName("Resonance") );
  resonanceSlider->setSliderName("Resonance");
  resonanceSlider->setDescription("Amount of feedback");
  resonanceSlider->setDescriptionField(infoField);
  resonanceSlider->setStringConversionFunction(&valueToStringTotal5);

  addWidget( spreadSlider = new ModulatableSlider() );
  spreadSlider->assignParameter( ladderToEdit->getParameterByName("StereoSpread") );
  spreadSlider->setSliderName("Spread");
  spreadSlider->setDescription("Detunes cutoff frequencies of channels");
  spreadSlider->setDescriptionField(infoField);
  spreadSlider->setStringConversionFunction(&valueToStringTotal5);

  addWidget( modeComboBox = new AutomatableComboBox() );
  modeComboBox->assignParameter( ladderToEdit->getParameterByName("Mode") );
  modeComboBox->setDescription("Select frequency response type");
  modeComboBox->setDescriptionField(infoField);
  modeComboBox->registerComboBoxObserver(this); // to update plot when mode is switched
   // todo: pass the mode parameter to the frequency response plot, too - when use does 
   // right-click, the popup opens and the mode can be selected. then we don't need to observe
   // the combobox here...
   // maybe, we can allow also different modes for the two filters and different resonances

  // set up the widgets:
  updateWidgetsAccordingToState();

  setSize(430, 251);  
  //setSize(420, 240);
  // There's a bug somewhere - if we request a size of 400 x 300, the widgets become invisible.
  // I guess, 400 x 300 is the size, it already has by default, and when we set it again to this 
  // size, it will be recognized that nothing changed and some update function is not being called
  // ...or something. ...need to figure that out
  // ..i think, resized() is not being called and thus, the widget positions won't be set up.
  // Maybe, the default size should not be 400 x 300 but something like 0 x 0 or (1 x 1, if that's
  // not allowed) - something that is assured to be never requested as actual size of an editor.
}

void LadderEditor::resized()
{
  ScopedLock scopedLock(*lock);
  AudioModuleEditor::resized();

  int x  = 0;
  int y  = 0;
  int w  = getWidth();
  int w2 = w/2;
  int h  = getHeight();

  y = getPresetSectionBottom();
  plotEditor->setBounds(x, y+4, w, h-y-2*20-8); // 2*20 for 2 widget-rows below it
  y = plotEditor->getBottom();

  cutoffSlider->setBounds(      4, y+4, w2-4, 16);
  resonanceSlider->setBounds(w2+4, y+4, w2-8, 16);

  y = cutoffSlider->getBottom();  

  modeComboBox->setBounds(   4, y+4, w2-4, 16);
  spreadSlider->setBounds(w2+4, y+4, w2-8, 16);
}

void LadderEditor::rComboBoxChanged(RComboBox* comboBoxThatHasChanged)
{
  plotEditor->repaint();
}