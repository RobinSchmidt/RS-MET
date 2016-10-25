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

  p = new AutomatableParameter(plugInLock, "Resonance", 0.0, 4.0, 0.0, 1.0, 
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

//=================================================================================================
// the associated GUI editor class:

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

  addWidget( cutoffSlider = new RSlider ("CutoffSlider") );
  cutoffSlider->assignParameter( ladderToEdit->getParameterByName("Cutoff") );
  cutoffSlider->setSliderName(juce::String("Cutoff"));
  cutoffSlider->setDescription(juce::String("Cutoff frequency in Hz"));
  cutoffSlider->setDescriptionField(infoField);
  cutoffSlider->setStringConversionFunction(&hertzToStringWithUnitTotal5); // use a hzToString function

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

  setSize(440, 128);
}

void LadderEditor::resized()
{
  ScopedLock scopedLock(*plugInLock);
  AudioModuleEditor::resized();
  int x = 0;
  int y = 0;
  int w = getWidth();
  int h = getHeight();

  y = getPresetSectionBottom()+8;

  cutoffSlider->setBounds(x+4, y, w-8, 16);
  y = cutoffSlider->getBottom();  
  resonanceSlider->setBounds(x+4, y+4, w-8, 16);
  y = resonanceSlider->getBottom();  
  spreadSlider->setBounds(x+4, y+4, w-8, 16);
  y = spreadSlider->getBottom(); 
  modeComboBox->setBounds(x+4, y+4, w/2-8, 16);

  // maybe here, we somehow have to also resize our wrapper, if any
}
