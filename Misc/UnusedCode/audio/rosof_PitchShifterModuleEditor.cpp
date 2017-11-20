#include "rosof_PitchShifterModuleEditor.h"
using namespace rosof;

//-------------------------------------------------------------------------------------------------
// construction/destruction:

PitchShifterModuleEditor::PitchShifterModuleEditor(CriticalSection *newPlugInLock, PitchShifterAudioModule* newPitchShifterAudioModule) 
: AudioModuleEditor(newPlugInLock, newPitchShifterAudioModule)
{
  ScopedLock scopedLock(*plugInLock);
    // maybe we should avoid this lock here and instead have a function that connects the widgets with the parameters where we acquire
    // the lock - but maybe not

  // set the plugIn-headline:
  setHeadlineText( juce::String(T("PitchShifter")) );

  // assign the pointer to the rosic::PitchShifter object to be used as aduio engine:
  jassert(newPitchShifterAudioModule != NULL ); // you must pass a valid module here
  pitchShifterModuleToEdit = newPitchShifterAudioModule;

  // create the widgets and assign the automatable parameters to them:
  addWidget( coarseSlider = new RSlider (T("CoarseSlider")) );
  coarseSlider->assignParameter( pitchShifterModuleToEdit->getParameterByName(T("DetuneCoarse")) );
  coarseSlider->setSliderName(juce::String(T("Coarse")));
  coarseSlider->setDescription(juce::String(T("Coarse pitch shifting factor in semitones")));
  coarseSlider->setDescriptionField(infoField);
  coarseSlider->setStringConversionFunction(&semitonesToStringWithUnit2);

  addWidget( fineSlider = new RSlider (T("FineSlider")) );
  fineSlider->assignParameter( pitchShifterModuleToEdit->getParameterByName(T("DetuneFine")) );
  fineSlider->setSliderName(juce::String(T("Fine")));
  fineSlider->setDescription(juce::String(T("Fine pitch shifting factor in cents")));
  fineSlider->setDescriptionField(infoField);
  fineSlider->setStringConversionFunction(&centsToStringWithUnit2);

  addWidget( grainLengthInMillisecondsSlider = new RSlider (T("GrainLengthSlider")) );
  grainLengthInMillisecondsSlider->assignParameter( 
    pitchShifterModuleToEdit->getParameterByName(T("GrainLengthInMilliseconds")) );
  grainLengthInMillisecondsSlider->setSliderName(juce::String(T("Grain Length")));
  grainLengthInMillisecondsSlider->setDescription(juce::String(T("Length of the grains in milliseconds")));
  grainLengthInMillisecondsSlider->setDescriptionField(infoField);
  grainLengthInMillisecondsSlider->setStringConversionFunction(&valueToStringTotal5);

  addWidget( grainLengthInCyclesSlider = new RSlider (T("CyclesPerGrainSlider")) );
  grainLengthInCyclesSlider->assignParameter( 
    pitchShifterModuleToEdit->getParameterByName(T("GrainLengthInPitchCycles")) );
  grainLengthInCyclesSlider->setSliderName(juce::String(T("Grain Length")));
  grainLengthInCyclesSlider->setDescription(juce::String(T("Length of the grains in pitch cylces")));
  grainLengthInCyclesSlider->setDescriptionField(infoField);
  grainLengthInCyclesSlider->setStringConversionFunction(&valueToStringTotal5);

  addWidget( grainLengthInBeatsSlider = new RSlider (T("GrainLengthInBeatsSlider")) );
  grainLengthInBeatsSlider->assignParameter( 
    pitchShifterModuleToEdit->getParameterByName(T("GrainLengthInBeats")) );
  grainLengthInBeatsSlider->setSliderName(juce::String(T("Grain Length")));
  grainLengthInBeatsSlider->setDescription(juce::String(T("Length of the grains in beats")));
  grainLengthInBeatsSlider->setDescriptionField(infoField);
  grainLengthInBeatsSlider->setStringConversionFunction(&valueToStringTotal5);

  addWidget( grainLengthUnitComboBox = new RComboBox(juce::String(T("GrainLengthUnitComboBox"))) );
  grainLengthUnitComboBox->assignParameter( 
    pitchShifterModuleToEdit->getParameterByName(T("GrainLengthUnit")) );
  grainLengthUnitComboBox->setDescription(T("Choose the unit for the grain length"));
  grainLengthUnitComboBox->setDescriptionField(infoField);
  grainLengthUnitComboBox->registerComboBoxObserver(this); // to update visibility of the sliders

  addWidget( feedbackSlider = new RSlider (T("FeedbackSlider")) );
  feedbackSlider->assignParameter( pitchShifterModuleToEdit->getParameterByName(T("Feedback")) );
  feedbackSlider->setSliderName(juce::String(T("Feedback")));
  feedbackSlider->setDescription(juce::String(T("Feeds the pitch-shifted output back to the input")));
  feedbackSlider->setDescriptionField(infoField);
  feedbackSlider->setStringConversionFunction(&percentToStringWithUnit1);

  addWidget( dryWetSlider = new RSlider (T("DryWet")) );
  dryWetSlider->assignParameter( pitchShifterModuleToEdit->getParameterByName(T("DryWet")) );
  dryWetSlider->setSliderName(juce::String(T("Dry/Wet")));
  dryWetSlider->setDescription(juce::String(T("Ratio between dry and wet signal (in % wet)")));
  dryWetSlider->setDescriptionField(infoField);
  dryWetSlider->setStringConversionFunction(&percentToStringWithUnit1);

  addWidget( antiAliasButton = new RButton(juce::String(T("Anti-Alias"))) );
  antiAliasButton->assignParameter( pitchShifterModuleToEdit->getParameterByName(T("AntiAlias")) );
  antiAliasButton->setDescription(juce::String(T("Switch anti-alias filter (for up-shifting) on/off")));
  antiAliasButton->setDescriptionField(infoField);
  antiAliasButton->setClickingTogglesState(true);

  addWidget( reverseButton = new RButton(juce::String(T("Reverse"))) );
  reverseButton->assignParameter( pitchShifterModuleToEdit->getParameterByName(T("Reverse")) );
  reverseButton->setDescription(juce::String(T("Reverse playback of the grains")));
  reverseButton->setDescriptionField(infoField);
  reverseButton->setClickingTogglesState(true);

  addWidget( invertButton = new RButton(juce::String(T("Invert"))) );
  invertButton->assignParameter( pitchShifterModuleToEdit->getParameterByName(T("Invert")) );
  invertButton->setDescription(juce::String(T("Invert polarity of wet (shifted) signal")));
  invertButton->setDescriptionField(infoField);
  invertButton->setClickingTogglesState(true);

  /*
  addWidget( formantPreserveButton = new RButton(juce::String(T("Formant"))) );
  formantPreserveButton->assignParameter( 
    pitchShifterModuleToEdit->getParameterByName(T("FormantPreserve")) );
  formantPreserveButton->setDescription(juce::String(T("Preserve formants")));
  formantPreserveButton->setDescriptionField(infoField);
  formantPreserveButton->setClickingTogglesState(true);

  addWidget( monoButton = new RButton(juce::String(T("Mono"))) );
  monoButton->assignParameter( pitchShifterModuleToEdit->getParameterByName(T("Mono")) );
  //monoButton->addRButtonListener(this);
  monoButton->setDescription(juce::String(T("Save CPU for mono signals")));
  monoButton->setDescriptionField(infoField);
  monoButton->setClickingTogglesState(true);
  */

  // set up the widgets:
  updateWidgetsAccordingToState();
}

//-------------------------------------------------------------------------------------------------
// callbacks:

void PitchShifterModuleEditor::rComboBoxChanged(RComboBox *rComboBoxThatHasChanged)
{
  ScopedLock scopedLock(*plugInLock);
  updateWidgetVisibility();
}

void PitchShifterModuleEditor::updateWidgetsAccordingToState()
{
  ScopedLock scopedLock(*plugInLock);
  AudioModuleEditor::updateWidgetsAccordingToState();
  updateWidgetVisibility();
}

void PitchShifterModuleEditor::resized()
{
  ScopedLock scopedLock(*plugInLock);
  AudioModuleEditor::resized();
  int x = 0;
  int y = 0;
  int w = getWidth();
  int h = getHeight();

  x  = 0;
  w /= 2;
  y = getPresetSectionBottom()+8;

  coarseSlider->setBounds(x+4, y, w-8, 16);

  y = coarseSlider->getBottom();  
  fineSlider->setBounds(x+4, y+4, w-8, 16);

  y = fineSlider->getBottom();  
  feedbackSlider->setBounds(x+4, y+4, w-8, 16);

  y = feedbackSlider->getBottom();  
  dryWetSlider->setBounds(x+4, y+4, w-8, 16);

  x = w;
  y = coarseSlider->getY();
 
  grainLengthInBeatsSlider->setBounds(x+4, y, w-64, 16);
  grainLengthInMillisecondsSlider->setBounds(grainLengthInBeatsSlider->getBounds());
  grainLengthInCyclesSlider->setBounds(grainLengthInBeatsSlider->getBounds());
  grainLengthUnitComboBox->setBounds(grainLengthInBeatsSlider->getRight()+4, y, 
    w-grainLengthInBeatsSlider->getWidth()-12, 16);

  y = grainLengthInMillisecondsSlider->getBottom()+8; 

  reverseButton->setBounds(x+4,    y, w/2-8, 16);
  invertButton->setBounds(x+w/2+4, y, w/2-8, 16);
  y += 24;
  //formantPreserveButton->setBounds(x+4, y+4, w/2-8, 16);
  x = reverseButton->getX() + reverseButton->getWidth()/2;
  w = invertButton->getX()  + invertButton->getWidth()/2   - x;
  antiAliasButton->setBounds(x, y+4, w, 16);

  //infoLabel->setBounds(0, getHeight()-20, 40, 20);
  infoField->setBounds(4, getHeight()-20, getWidth()-4,20);
  webLink->setBounds(getWidth()-112, getHeight()-20, 112-4, 20);
}

void PitchShifterModuleEditor::updateWidgetVisibility()
{
  ScopedLock scopedLock(*plugInLock);
  if( pitchShifterModuleToEdit == NULL )
    return;
  if( pitchShifterModuleToEdit->wrappedPitchShifter == NULL )
    return;

  // update the visibility for the 3 grain-length sliders:
  grainLengthInMillisecondsSlider->setVisible(false);  
  grainLengthInCyclesSlider->setVisible(false);
  grainLengthInBeatsSlider->setVisible(false);
  switch( pitchShifterModuleToEdit->wrappedPitchShifter->getGrainLengthUnit() )
  {
  case PitchShifterGrainAdaptive::MILLISECONDS: grainLengthInMillisecondsSlider->setVisible(true);  break;
  case PitchShifterGrainAdaptive::PITCH_CYCLES: grainLengthInCyclesSlider->setVisible(true);        break;
  case PitchShifterGrainAdaptive::BEATS:        grainLengthInBeatsSlider->setVisible(true);         break;
  }
}