
Ladder::Ladder()
{
  ScopedLock scopedLock(*plugInLock);
  moduleName = juce::String("Ladder");
  setActiveDirectory(getApplicationDirectory() + juce::String("/LadderPresets") );

  createStaticParameters();
}

void Ladder::createStaticParameters()
{
  ScopedLock scopedLock(*plugInLock);

  juce::Array<double> defaultValues;
  AutomatableParameter* p;

  p = new AutomatableParameter(plugInLock, "Cutoff", 20.0, 20000.0, 0.0, 1000.0, 
    Parameter::EXPONENTIAL, 74);
  defaultValues.clear();
  defaultValues.add(125.0);
  defaultValues.add(250.0);
  defaultValues.add(500.0);
  defaultValues.add(1000.0);
  defaultValues.add(2000.0);
  defaultValues.add(4000.0);
  defaultValues.add(8000.0);
  defaultValues.add(16000.0);
  p->setDefaultValues(defaultValues);
  addObservedParameter(p);
  p->setValueChangeCallback<Ladder>(this, &Ladder::setCutoff);

  p = new AutomatableParameter(plugInLock, "Resonance", 0.0, 2.0, 0.0, 1.0, 
    Parameter::LINEAR, 71);
  addObservedParameter(p);
  p->setValueChangeCallback<Ladder>(this, &Ladder::setResonance);

  p = new AutomatableParameter(plugInLock, "StereoSpread", -24.0, 24.0, 0.0, 0.0, 
    Parameter::LINEAR_BIPOLAR);
  addObservedParameter(p);
  p->setValueChangeCallback<Ladder>(this, &Ladder::setStereoSpread);

  p = new AutomatableParameter(plugInLock, "Mode", 0.0, 14.0, 4.0, 0.0, Parameter::STRING);
  p->addStringValue(juce::String("Flat"));
  p->addStringValue(juce::String("Lowpass 6 dB/oct"));
  p->addStringValue(juce::String("Lowpass 12 dB/oct"));
  p->addStringValue(juce::String("Lowpass 18 dB/oct"));
  p->addStringValue(juce::String("Lowpass 24 dB/oct"));
  p->addStringValue(juce::String("Highpass 6 dB/oct"));
  p->addStringValue(juce::String("Highpass 12 dB/oct"));
  p->addStringValue(juce::String("Highpass 18 dB/oct"));
  p->addStringValue(juce::String("Highpass 24 dB/oct"));
  p->addStringValue(juce::String("Bandpass 6/6 dB/oct"));
  p->addStringValue(juce::String("Bandpass 6/12 dB/oct"));
  p->addStringValue(juce::String("Bandpass 6/18 dB/oct"));
  p->addStringValue(juce::String("Bandpass 12/6 dB/oct"));
  p->addStringValue(juce::String("Bandpass 12/12 dB/oct"));
  p->addStringValue(juce::String("Bandpass 18/6 dB/oct"));
  addObservedParameter(p);
  p->setValueChangeCallback<Ladder>(this, &Ladder::setMode);

  // make sure that the parameters are initially in sync with the audio engine:
  for(int i = 0; i < (int)observedParameters.size(); i++)
    observedParameters[i]->resetToDefaultValue(true, true);

  // there seems to be bug - when saving a preset, it saves a controller-mapping of (nonexistent)
  // controller number -1 for the mode parameter -> check this...
}

//-------------------------------------------------------------------------------------------------
// Editor creation:

template<class ProcessorType, class EditorType>
AudioProcessorEditor* createAudioProcessorEditor(ProcessorType processor, EditorType *dummy)
{
  jura::ParameterObserver::guiAutomationSwitch = false;  // don't automate widgets during creation
  EditorType* editor = new EditorType(processor);
  AudioProcessorEditor *wrapper = new AudioModuleEditorWrapper(editor, processor);
  jura::ParameterObserver::guiAutomationSwitch = true;   // now, widgets can be automated again
  return wrapper;
}
AudioProcessorEditor* Ladder::createEditor()
{
  LadderEditor *dummy = nullptr;
  return createAudioProcessorEditor(this, dummy);
}
// The template function createAudioProcessorEditor and the createEditor function that makes use of 
// the template are a little trick to let the compiler generate the equivalent of the following 
// code:
// AudioProcessorEditor* Ladder::createEditor()
// {
//   jura::ParameterObserver::guiAutomationSwitch = false;  // don't automate widgets during creation
//   LadderEditor* editor = new LadderEditor(this);
//   AudioProcessorEditor *wrapper = new AudioModuleEditorWrapper(editor, this);
//   jura::ParameterObserver::guiAutomationSwitch = true;   // now, widgets can be automated again
//   return wrapper;
// }
// The advantage of doing it that way is, that the createAudioProcessorEditor template can be used 
// for any other pairs of processor/editor classes such that the body of the function (this 
// guiAutomation deactivation stuff etc.) does not need to be rewritten. We get rid of a lot of 
// identical boilerplate code which would otherwise have to be written out. Now, we only need to 
// write a two-liner for each processor/editor pair of classes instead of a five-liner. That 
// makes it also easier to change the content of the editor creation function for all classes at 
// once in one single location and thereby makes the code more maintenance friendly.

// However, in order to make the createAudioProcessorEditor template available for the creation of 
// other editors, we need to move it to another file which the other files that want to use it
// will include.....later...

//-------------------------------------------------------------------------------------------------
// the audio processing callback:

void Ladder::processBlock(double **inOutBuffer, int numChannels, int numSamples)
{
  jassert(numChannels == 2);
  for(int n = 0; n < numSamples; n++)
  {
    inOutBuffer[0][n] = ladderL.getSample(inOutBuffer[0][n]);
    inOutBuffer[1][n] = ladderR.getSample(inOutBuffer[1][n]);
  }
}

//-------------------------------------------------------------------------------------------------
// parameter setters (callback targets for the Parameter objects):

void Ladder::setSampleRate(double newSampleRate){
  ladderL.setSampleRate(newSampleRate); 
  ladderR.setSampleRate(newSampleRate);
}
void Ladder::setCutoff(double newCutoff){
  ladderL.setCutoff(freqFactorL * newCutoff); 
  ladderR.setCutoff(freqFactorR * newCutoff);
}
void Ladder::setResonance(double newResonance){
  ladderL.setResonance(newResonance);
  ladderR.setResonance(newResonance);
}
void Ladder::setMode(int newMode){
  ladderL.setMode(newMode);
  ladderR.setMode(newMode);
}
void Ladder::reset(){
  ladderL.reset();
  ladderR.reset();
}

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

double Ladder::getMagnitudeAt(double frequency)
{
  return 1.0; // not yet implemented

  /*
  computing the magnitude response should go something like:
  double w = 2*PI*frequency/sampleRate;
  complex<double> j(0, 1);            // imaginary unit
  complex<double> z = exp(j*w);
  complex<double> G0, G1, G2, G3, G4; // transfer functions of n-th stage output, n = 0..4
  complex<double> H, H2;              // transfer function with resonance and its magnitude-squared
  G0 = 1;
  G1 = b / (1 + a/z);
  G2 = G1*G1;
  G3 = G1*G2;
  G4 = G2*G2;
  H  = (c0*G0 + c1*G1 + c2*G2 + c3*G3 + c4*G4) / (1 + k * G4 / z);
  H2 = H * conj(H);          // magnitude squared
  return sqrt( H2.real() );  // magnitude

  // ..i think - i need to check the math..
  //ladderL.getMagnitudeAt(z);
  // we should do that in in RAPT::Ladder
  // i wonder, if a simpler formula is possible which avoids going through the complex transfer 
  // function -> computer algebra
  */
}

void Ladder::getMagnitudeResponse(const double *frequencies, double *magnitudes, int numBins,
  bool inDecibels)
{
  int i;
  for(i = 0; i < numBins; i++)
    magnitudes[i] = getMagnitudeAt(frequencies[i]);
  if(inDecibels)
    for(i = 0; i < numBins; i++)
      magnitudes[i] = amp2dBWithCheck(magnitudes[i], 1.e-6);
}

//=================================================================================================
// LadderSpectrumEditor:

// construction/destruction:

LadderSpectrumEditor::LadderSpectrumEditor(const juce::String& name) 
  : SpectrumDisplayOld(name)
{
  // nulling these pointers must occur before setting up the plot range:
  freqParameter = nullptr;
  resoParameter = nullptr;
  filterToEdit  = nullptr;

  // this stuff will be (re-) allocated in resized():
  numBins     = 0;
  frequencies = NULL;
  magnitudes  = NULL;


  setDescription("Drag around the node to adjust the filter's frequency and resonance");
  ParameterObserver::isGuiElement = true;

  // set up the plot range:
  setMaximumRange(15.625, 32000.0, -60.0, 30.0);
  setCurrentRange(15.625, 32000.0, -60.0, 30.0);
  setHorizontalCoarseGrid(12.0, false);
  setHorizontalFineGrid(   3.0, false);
  setVerticalCoarseGridVisible( false);
  setVerticalFineGridVisible(   false);

  plotColourScheme.setCurveColouringStrategy(PlotColourScheme::UNIFORM);

  currentMouseCursor = MouseCursor(MouseCursor::NormalCursor);
  setMouseCursor(currentMouseCursor);

  dotRadius = 5.f;

  // activate automation for this ParameterObserver:
  ParameterObserver::localAutomationSwitch = true;
}

LadderSpectrumEditor::~LadderSpectrumEditor(void)
{
  // remove ourselves as listeners from the Parameter objects, such that they do not try to notify 
  // a nonexistent listener:
  ParameterObserver::localAutomationSwitch = false;
  if( freqParameter != NULL )
    freqParameter->deRegisterParameterObserver(this);
  if( resoParameter != NULL )
    resoParameter->deRegisterParameterObserver(this);

  deleteAndZero(frequencies);
  deleteAndZero(magnitudes);
}

//-------------------------------------------------------------------------------------------------
// parameter-settings:

void LadderSpectrumEditor::setFilterToEdit(jura::Ladder* newFilterToEdit)
{
  filterToEdit = newFilterToEdit;
}

void LadderSpectrumEditor::assignParameterFreq(Parameter *parameterToAssign)
{
  freqParameter = parameterToAssign;
  if( freqParameter != NULL )
    freqParameter->registerParameterObserver(this);
}

void LadderSpectrumEditor::assignParameterReso(Parameter *parameterToAssign)
{
  resoParameter = parameterToAssign;
  if( resoParameter != NULL )
    resoParameter->registerParameterObserver(this);
}

void LadderSpectrumEditor::unAssignParameterFreq()
{
  if( freqParameter != NULL )
    freqParameter->deRegisterParameterObserver(this);
  freqParameter = NULL;
}

void LadderSpectrumEditor::unAssignParameterReso()
{
  if( resoParameter != NULL )
    resoParameter->deRegisterParameterObserver(this);
  resoParameter = NULL;
}

void LadderSpectrumEditor::parameterChanged(Parameter* parameterThatHasChanged)
{
  //// old code:
  //// send out a change-message: 
  //sendChangeMessage(); // ? why that?
  ////updatePlot();

  // new:
  updatePlot();
}

void LadderSpectrumEditor::parameterIsGoingToBeDeleted(Parameter* parameterThatWillBeDeleted)
{
  // clear reference to parameter that will be deleted
}

void LadderSpectrumEditor::updateWidgetFromAssignedParameter(bool sendMessage)
{
  updatePlot();
  if( sendMessage == true )
    sendChangeMessage();
}

//-------------------------------------------------------------------------------------------------
// callbacks:

void LadderSpectrumEditor::changeListenerCallback(ChangeBroadcaster *objectThatHasChanged)
{
  // temporarily switch the wantsAutomationNotification flag from the ParameterObserver base 
  // class off to avoid circular notifications and updates:
  localAutomationSwitch = false;

    // ?? what is this good for ?? we don't send out any messages here... fihure out in debugger

  // call the method which updates the widget:
  updatePlot();
  //updateWidgetFromAssignedParameter(false);

  // switch the wantsAutomationNotification flag on again:  
  localAutomationSwitch = true;
}

void LadderSpectrumEditor::mouseDown(const MouseEvent &e)
{
  if( filterToEdit == NULL )
    return;

  // preliminary: do not open the MIDI-learn menu on right-button - show the image export menu 
  // instead (inherited behaviour from CoordinateSytem):
  if( e.mods.isRightButtonDown() )
    CoordinateSystemOld::mouseDown(e);
  else
  {
    // get the position of the event in components coordinates
    double x = e.getMouseDownX();
    double y = e.getMouseDownY();

    setupFilterAccordingToMousePosition(x, y);

    // send out a change-message:
    sendChangeMessage();
  }
}

void LadderSpectrumEditor::mouseDrag(const juce::MouseEvent &e)
{
  if( filterToEdit == NULL )
    return;

  /*
  if( e.mods.isRightButtonDown() && xParameter != NULL && yParameter != NULL )
  {
  // ignore mouse drags whne the right button is down and we have assigned automatable 
  // parameters because in that case, the right click is used for opening the MIDI-learn popup
  }
  else...
  */

  // get the position of the event in components coordinates
  double x = e.getMouseDownX() + e.getDistanceFromDragStartX();
  double y = e.getMouseDownY() + e.getDistanceFromDragStartY();

  setupFilterAccordingToMousePosition(x, y);

  // send out a change-message:
  sendChangeMessage();
}

void LadderSpectrumEditor::setupFilterAccordingToMousePosition(double mouseX, double mouseY)
{
  if( filterToEdit == NULL )
    return;

  // get the position of the event in components coordinates
  double x = mouseX;
  double y = mouseY;

  // convert them into a frequency and resonance/q/gain value:
  double gain = y;
  transformFromComponentsCoordinates(x, gain);
  gain = clip(gain, -60.0, 30.0);
  double freq = x;
  double reso = yToReso(y);

  // restrict ranges (ToDo: actually the filter should take care of the itself....):
  //freq = clip(freq, 20.0, 20000.0);
  //reso = clip(reso,  0.0,     1.0);   // change the upper limit later...

  // set up the filter and raise automation events to update other widgets that represent the
  // parameters:
  if(freqParameter != NULL)
    freqParameter->setValue(freq, true, true);
  if(resoParameter != NULL)
    resoParameter->setValue(reso, true, true);
}

//-----------------------------------------------------------------------------------------------------------------------------------------
// drawing:

void LadderSpectrumEditor::resized()
{
  SpectrumDisplayOld::resized();

  // (re) allocate and fill the arrays for the magnitude plot
  numBins = getWidth();
  if( frequencies == NULL )
    delete[] frequencies;
  if( magnitudes == NULL )
    delete[] magnitudes;
  frequencies = new double[numBins];
  magnitudes  = new double[numBins];
  getDisplayedFrequencies(frequencies, numBins);
  for(int i = 0; i < numBins; i++)
    magnitudes[i]  = 0.0;
  updatePlot();
}

inline void clipBuffer(double *buffer, int length, double min, double max)
{
  for(int i = 0; i < length; i++)
    buffer[i] = clip(buffer[i], min, max);
}
void LadderSpectrumEditor::updatePlot()
{
  if( filterToEdit == NULL )
    return;

  // fill the magnitude array with the magnitudes:
  filterToEdit->getMagnitudeResponse(frequencies, magnitudes, numBins, true);
  clipBuffer(magnitudes, numBins, -120.0, 120.0);

  // overwrite the magnitude value at the bin closest to the cutoff-frequency with the magnitude at 
  // the exact cutoff frequency:
  double freq  = filterToEdit->getCutoff();
  double level = amp2dBWithCheck(filterToEdit->getMagnitudeAt(freq), 0.000001);
  level        = clip(level, -120.0, +120.0);

  for(int bin=0; bin < (numBins-1); bin++)
  {
    if( frequencies[bin] < freq && frequencies[bin+1] > freq )
    {
      if( fabs(frequencies[bin]-freq) <= fabs(frequencies[bin+1]-freq) ) // lower bin is closer
        magnitudes[bin]   = level;
      else                                                               // upper bin is closer
        magnitudes[bin+1] = level;
    }
  }

  setSpectrum(numBins, frequencies, magnitudes);
  //repaint();
}

void LadderSpectrumEditor::plotCurveFamily(Graphics &g, Image* targetImage, 
  XmlElement *targetSVG)
{
  //if( filterToEdit == nullptr )
  //  return;
  if(freqParameter == nullptr || resoParameter == nullptr) // pointers need to be assigned
    return;

  CurveFamilyPlotOld::plotCurveFamily(g, targetImage, targetSVG);

  //Colour graphColour = colourScheme.curves; // preliminary
  //if( colourScheme.plotColours.size() > 0 )
  //  graphColour = colourScheme.plotColours[0];
  Colour graphColour = plotColourScheme.getCurveColour(0);  
  //g.setColour(graphColour); 

  double x = freqParameter->getValue();  // frequency in physical
  double y = 0;                          // dummy, for the moment

  // determine the coordinates of the handle into image component coordinates (for export) or 
  // components coordinates for plot:
  if( targetImage == NULL )
    transformToComponentsCoordinates(x, y);
  else
    transformToImageCoordinates(x, y, targetImage);

  y = resoToY(resoParameter->getValue());


  // draw the handle and a crosshair:
  g.fillEllipse((float) (x-dotRadius), (float) (y-dotRadius), 
    (float) (2*dotRadius), (float) (2*dotRadius) );
  g.setColour(graphColour.withAlpha(0.4f));
  float w = (float) getPlotWidth(targetImage);
  float h = (float) getPlotHeight(targetImage);
  g.drawLine(       0,(float)y,        w, (float)y, 1.f);  // horizontal
  g.drawLine((float)x,       0, (float)x,        h, 1.f);  // vertical
}

// internal functions:

double LadderSpectrumEditor::resoToY(double reso, Image *targetImage)
{
  double maxRes = resoParameter->getMaxValue();
  //double y = (maxRes - reso/maxRes) * getPlotHeight(targetImage);
  double y = (1 - reso/maxRes) * getPlotHeight(targetImage);
  return y;
}

double LadderSpectrumEditor::yToReso(double y, Image *targetImage)
{
  double maxRes = resoParameter->getMaxValue();
  double reso = maxRes * ( 1.0 - y/getPlotHeight(targetImage) );  
  //double reso = maxRes * ( maxRes - y/getPlotHeight(targetImage) );  
  return reso;
}

//=================================================================================================
// the GUI editor class for the Ladder:

LadderEditor::LadderEditor(jura::Ladder *newLadderToEdit) : AudioModuleEditor(newLadderToEdit)
{
  ScopedLock scopedLock(*plugInLock);
  // maybe we should avoid this lock here and instead have a function that connects the widgets 
  // with the parameters where we acquire the lock - but maybe not

  // set the plugIn-headline:
  setHeadlineText( juce::String("Ladder") );

  // assign the pointer to the edited object:
  jassert(newLadderToEdit != nullptr ); // you must pass a valid module here
  ladderToEdit = newLadderToEdit;

  // create the widgets and assign the automatable parameters to them:

  frequencyResponseDisplay = new LadderSpectrumEditor(juce::String("SpectrumEditor"));
  frequencyResponseDisplay->setFilterToEdit(ladderToEdit);
  frequencyResponseDisplay->addChangeListener(this); // do we need this?
  frequencyResponseDisplay->assignParameterFreq( moduleToEdit->getParameterByName("Cutoff"));
  frequencyResponseDisplay->assignParameterReso( moduleToEdit->getParameterByName("Resonance"));
  addPlot( frequencyResponseDisplay );

  addWidget( cutoffSlider = new RSlider ("CutoffSlider") );
  cutoffSlider->assignParameter( ladderToEdit->getParameterByName("Cutoff") );
  cutoffSlider->setSliderName(juce::String("Cutoff"));
  cutoffSlider->setDescription(juce::String("Cutoff frequency in Hz"));
  cutoffSlider->setDescriptionField(infoField);
  cutoffSlider->setStringConversionFunction(&hertzToStringWithUnitTotal5);

  addWidget( resonanceSlider = new RSlider ("ResoSlider") );
  resonanceSlider->assignParameter( ladderToEdit->getParameterByName("Resonance") );
  resonanceSlider->setSliderName(juce::String("Resonance"));
  resonanceSlider->setDescription(juce::String("Amount of feedback"));
  resonanceSlider->setDescriptionField(infoField);
  resonanceSlider->setStringConversionFunction(&valueToStringTotal5);

  addWidget( spreadSlider = new RSlider ("SpreadSlider") );
  spreadSlider->assignParameter( ladderToEdit->getParameterByName("StereoSpread") );
  spreadSlider->setSliderName(juce::String("Spread"));
  spreadSlider->setDescription(juce::String("Detunes cutoff frequencies of channels"));
  spreadSlider->setDescriptionField(infoField);
  spreadSlider->setStringConversionFunction(&valueToStringTotal5);

  addWidget( modeComboBox = new RComboBox(juce::String("ModeComboBox")) );
  modeComboBox->assignParameter( ladderToEdit->getParameterByName("Mode") );
  modeComboBox->setDescription("Select frequency response type");
  modeComboBox->setDescriptionField(infoField);

  // change to midSideButton
  //addWidget( invertButton = new RButton(juce::String(T("Invert"))) );
  //invertButton->assignParameter( pitchShifterModuleToEdit->getParameterByName(T("Invert")) );
  //invertButton->setDescription(juce::String(T("Invert polarity of wet (shifted) signal")));
  //invertButton->setDescriptionField(infoField);
  //invertButton->setClickingTogglesState(true);

  // set up the widgets:
  updateWidgetsAccordingToState();

  setSize(420, 240);  
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
  ScopedLock scopedLock(*plugInLock);
  AudioModuleEditor::resized();

  // arrange widgets as follows:
  // the plot, below it the cutoff slider over th full width, below the mode-box and
  // reso slider over the half width

  int x = 0;
  int y = 0;
  int w = getWidth();
  int h = getHeight();

  //y = getPresetSectionBottom()+8;

  y = getPresetSectionBottom();
  frequencyResponseDisplay->setBounds(x, y+4, w, h-92);
  y = frequencyResponseDisplay->getBottom();

  cutoffSlider->setBounds(x+4, y+4, w-8, 16);
  y = cutoffSlider->getBottom();  
  resonanceSlider->setBounds(x+4, y+4, w-8, 16);
  y = resonanceSlider->getBottom();  
  spreadSlider->setBounds(x+4, y+4, w-8, 16);
  y = spreadSlider->getBottom(); 
  modeComboBox->setBounds(x+4, y+4, w/2-8, 16);


  // maybe here, we somehow have to also resize our wrapper, if any
}
