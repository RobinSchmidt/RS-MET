//#include "rosof_AcidSequencerAudioModule.h"
//using namespace rosof;

//-------------------------------------------------------------------------------------------------
// construction/destruction:

AcidSequencerAudioModule::AcidSequencerAudioModule(CriticalSection *newPlugInLock,
  rosic::AcidSequencer *acidSequencerToWrap) : AudioModule(newPlugInLock)
{
  jassert(acidSequencerToWrap != NULL); // you must pass a valid rosic-object to the constructor
  wrappedAcidSequencer = acidSequencerToWrap;
  setModuleTypeName("AcidSequencer");
  createParameters();
}

//-------------------------------------------------------------------------------------------------
// automation:

void AcidSequencerAudioModule::parameterChanged(Parameter* parameterThatHasChanged)
{
  if( wrappedAcidSequencer == NULL )
    return;

  int    index = getIndexOfParameter(parameterThatHasChanged);
  double value = parameterThatHasChanged->getValue();
  switch( index )
  {
  case   0: wrappedAcidSequencer->setStepLength(     value);   break;
  case   1: wrappedAcidSequencer->setMode(      (int)value);   break;
  default:
    {
      // do nothing
    }
  } // end of switch( parameterIndex )

}

void AcidSequencerAudioModule::setStateFromXml(const XmlElement &xmlState,
                                               const juce::String &stateName, bool markAsClean)
{
  if( wrappedAcidSequencer != NULL )
  {
    // set permissibilities...

    for(int p=0; p<wrappedAcidSequencer->getNumPatterns(); p++ )
    {
      AcidPattern *pattern = wrappedAcidSequencer->getPattern(p);
      pattern->clear();
      XmlElement *xmlPattern = xmlState.getChildByName(juce::String("Pattern") + juce::String(p));
      if( xmlPattern != NULL )
      {
        pattern->setStepLength(xmlPattern->getDoubleAttribute("StepLength", 0.5));
        for(int s=0; s<pattern->getMaxNumSteps(); s++)
        {
          AcidNote   *note    = pattern->getNote(s);
          XmlElement *xmlStep = xmlPattern->getChildByName(juce::String("Step") + juce::String(s));
          if( xmlStep != NULL )
          {
            note->gate   = xmlStep->getBoolAttribute("Gate",   false);
            note->accent = xmlStep->getBoolAttribute("Accent", false);
            note->slide  = xmlStep->getBoolAttribute("Slide",  false);
            note->key    = xmlStep->getIntAttribute( "Key",    0);
            note->octave = xmlStep->getIntAttribute( "Octave", 0);
          }
        }
      }
    }
  }

  AudioModule::setStateFromXml(xmlState, stateName, markAsClean);
}

XmlElement* AcidSequencerAudioModule::getStateAsXml(const juce::String &stateName, bool markAsClean)
{
  XmlElement *xmlState = AudioModule::getStateAsXml(stateName, markAsClean);

  if( wrappedAcidSequencer != NULL )
  {
    // retrieve permissibilities...


    for(int p=0; p<wrappedAcidSequencer->getNumPatterns(); p++ )
    {
      AcidPattern *pattern = wrappedAcidSequencer->getPattern(p);
      if( !pattern->isEmpty() )
      {
        XmlElement *xmlPattern = new XmlElement(juce::String("Pattern") + juce::String(p));
        xmlPattern->setAttribute("StepLength", pattern->getStepLength());
        for(int s=0; s<pattern->getMaxNumSteps(); s++)
        {
          AcidNote *note = pattern->getNote(s);
          if( !note->isInDefaultState() )
          {
            XmlElement *xmlStep = new XmlElement(juce::String("Step") + juce::String(s));
            xmlStep->setAttribute("Gate",   note->gate);
            xmlStep->setAttribute("Accent", note->accent);
            xmlStep->setAttribute("Slide",  note->slide);
            xmlStep->setAttribute("Key",    note->key);
            xmlStep->setAttribute("Octave", note->octave);
            xmlPattern->addChildElement(xmlStep);
          }
        }
        xmlState->addChildElement(xmlPattern);
      }
    }
  }

  return xmlState;
}

//-------------------------------------------------------------------------------------------------
// internal functions:

void AcidSequencerAudioModule::createParameters()
{
  AutomatableParameter* p;

  // #000:
  p = new AutomatableParameter(lock, "StepLength", 0.25, 1.0, 0.01, 0.5, Parameter::LINEAR);
  addObservedParameter(p);

  // #001:
  p = new AutomatableParameter(lock, juce::String("Mode"), 0.0, 1.0, 1.0, 1.0, Parameter::STRING);
  p->addStringValue(juce::String("Off"));
  p->addStringValue(juce::String("Key"));
  //p->addStringValue(juce::String(T("Host")));
  p->setValue(0.0, false, false);
  addObservedParameter(p);

  // make a call to parameterChanged for each parameter in order to set up the DSP-core to reflect
  // the values the automatable parameters:
  for(int i=0; i < (int) parameters.size(); i++ )
    parameterChanged(parameters[i]);
}

//=================================================================================================


AcidPatternEditor::AcidPatternEditor(rosic::AcidSequencer *sequencerToEdit)
  : Component(juce::String("AcidPatternEditor"))
{
  patternToEdit         = NULL;
  this->sequencerToEdit = sequencerToEdit;

  whiteKeyColour           = Colours::white;
  blackKeyColour           = Colours::black;
  backgroundColourWhiteKey = Colours::white;      // for white key lanes
  backgroundColourBlackKey = Colours::lightgrey;  // for black key lanes
  handleColor              = Colours::black;
  textColour               = Colours::black;
  lineColour               = Colours::black;

  rowHeight     = 12.f;
  columnWidth   = 20.f;
  keyLength     = 48.f;
  topLaneHeight = 16.f;
}

//-------------------------------------------------------------------------------------------------
// setup:

void AcidPatternEditor::setPatternToEdit(rosic::AcidPattern *newPatternToEdit)
{
  patternToEdit = newPatternToEdit;
}

/*
void AcidPatternEditor::setColourScheme(const WidgetColourScheme& newColourScheme)
{
  RWidget::setColourScheme(newColourScheme);

  // colors we can use:
  // getOutlineColour(), getHandleColour(), getTextColour(), getSpecialColour1()

  whiteKeyColour           = getHandleColour();
  blackKeyColour           = Colours::black;
  backgroundColourWhiteKey = Colours::white;      // for white key lanes
  backgroundColourBlackKey = Colours::lightgrey;  // for black key lanes
  handleColor              = getHandleColour();
  textColour               = getTextColour();
  lineColour               = getOutlineColour();
}
*/

//-------------------------------------------------------------------------------------------------
// inquiry:

int AcidPatternEditor::getStepAt(float x)
{
  if( x < keyLength )
    return -1;
  else
  {
    float x2 = x - keyLength;
    return floorInt(x2/columnWidth);
  }
}

int AcidPatternEditor::getNoteAt(float y)
{
  if( y < 4*topLaneHeight )
    return -1;
  else
  {
    float y2 = (float) getHeight() - y;
    return floorInt(y2/rowHeight);
  }
}

int AcidPatternEditor::getVirtualKeyAt(float x, float y)
{
  int key = getNoteAt(y);
  if( key != -1 && x <= keyLength )
    return key;
  else
    return -1;
}

bool AcidPatternEditor::isInGateRow(float y)
{
  if( y < topLaneHeight )
    return true;
  else
    return false;
}

bool AcidPatternEditor::isInAccentRow(float y)
{
  if( y < 2*topLaneHeight && y >= topLaneHeight )
    return true;
  else
    return false;
}

bool AcidPatternEditor::isInSlideRow(float y)
{
  if( y < 3*topLaneHeight && y >= 2*topLaneHeight )
    return true;
  else
    return false;
}

bool AcidPatternEditor::isInOctaveRow(float y)
{
  if( y < 4*topLaneHeight && y >= 3*topLaneHeight )
    return true;
  else
    return false;
}

bool AcidPatternEditor::isInKeyboardColumn(float x)
{
  if( x < keyLength )
    return true;
  else
    return false;
}

//-------------------------------------------------------------------------------------------------
// callbacks:

void AcidPatternEditor::mouseDown(const MouseEvent &e)
{
  if( patternToEdit == NULL )
    return;

  int step = getStepAt((float) e.x);
  int key  = getNoteAt((float) e.y);
  if( step != -1 && key != -1 )
  {
    if( patternToEdit->getKey(step) == key && patternToEdit->getGate(step) == true )
    {
      patternToEdit->setGate(step, false);
    }
    else
    {
      patternToEdit->setKey( step, key);
      patternToEdit->setGate(step, true);
    }
  }
  else if( isInGateRow((float) e.y) )
  {
    if( step != -1 )
      patternToEdit->setGate(step, !patternToEdit->getGate(step));
  }
  else if( isInAccentRow((float) e.y) )
  {
    if( step != -1 )
      patternToEdit->setAccent(step, !patternToEdit->getAccent(step));
  }
  else if( isInSlideRow((float) e.y) )
  {
    if( step != -1 )
      patternToEdit->setSlide(step, !patternToEdit->getSlide(step));
  }
  else if( isInOctaveRow((float) e.y) )
  {
    if( step != -1 )
    {
      if( e.mods.isLeftButtonDown() )
        patternToEdit->setOctave(step, patternToEdit->getOctave(step)-1);
      else if( e.mods.isRightButtonDown() )
        patternToEdit->setOctave(step, patternToEdit->getOctave(step)+1);
      else if( e.mods.isMiddleButtonDown() )
        patternToEdit->setOctave(step, 0);
    }
  }
  else if( key != -1 && isInKeyboardColumn((float) e.x) )
  {
    if( sequencerToEdit != NULL )
      sequencerToEdit->toggleKeyPermissibility(key);
  }

  repaint();
}

void AcidPatternEditor::paint(juce::Graphics &g)
{
  g.fillAll(Colours::white);

  float x = 0.f;
  float y = 0.f;
  float w = (float) getWidth();
  float h = (float) getHeight();
  //float s = (float) rowHeight;
  float thickness = 2.f;


  drawBitmapFontText(g, (int)x+3, (int)y+3, "Gate:", &BitmapFontRoundedBoldA10D0::instance, textColour);
  y += topLaneHeight;
  drawBitmapFontText(g, (int)x+3, (int)y+3, "Accent:", &BitmapFontRoundedBoldA10D0::instance, textColour);
  y += topLaneHeight;
  drawBitmapFontText(g, (int)x+3, (int)y+3, "Slide:", &BitmapFontRoundedBoldA10D0::instance, textColour);
  y += topLaneHeight;
  drawBitmapFontText(g, (int)x+3, (int)y+3, "Octave:", &BitmapFontRoundedBoldA10D0::instance, textColour);

  w = keyLength;

  // draw background for white keys and black keys on top of it:
  float keyboardY = 4*topLaneHeight;

  g.setColour(whiteKeyColour);
  g.fillRect(0.f, keyboardY, keyLength, (float)(13*rowHeight));

  g.setColour(backgroundColourWhiteKey);
  g.fillRect(keyLength, keyboardY, (float)getWidth()-keyLength, (float)(13*rowHeight));

  // draw black keys:
  g.setColour(blackKeyColour);
  x = 0;
  w = 2*keyLength/3;
  h = (float) rowHeight;
  y = keyboardY + 11*rowHeight;
  g.fillRect(x, y, w, h); y -= 2*h;
  g.fillRect(x, y, w, h); y -= 3*h;
  g.fillRect(x, y, w, h); y -= 2*h;
  g.fillRect(x, y, w, h); y -= 2*h;
  g.fillRect(x, y, w, h);

  // draw lanes for the black keys:
  x = keyLength;
  w = (float)getWidth()-keyLength;
  y = keyboardY + 11*rowHeight;
  g.setColour(backgroundColourBlackKey);
  g.fillRect(x, y, w, h);  y -= 2*h;
  g.fillRect(x, y, w, h);  y -= 3*h;
  g.fillRect(x, y, w, h);  y -= 2*h;
  g.fillRect(x, y, w, h);  y -= 2*h;
  g.fillRect(x, y, w, h);

  // draw horizontal lines for the 4 top-lanes:
  thickness = 2.f;
  x = 0.f;
  y = topLaneHeight;
  w = (float) getWidth();
  g.setColour(lineColour);
  g.drawLine(x, y, w, y, thickness);
  y += topLaneHeight;
  g.drawLine(x, y, w, y, thickness);
  y += topLaneHeight;
  g.drawLine(x, y, w, y, thickness);
  y += topLaneHeight;
  g.drawLine(x, y, w, y, thickness);
  g.drawLine(keyLength, 0, keyLength, (float)getHeight(), thickness); // vertical

  // draw the lines between the piano-roll rows:
  g.setColour(lineColour);
  x = keyLength;
  h = (float) rowHeight;
  y = keyboardY + h;
  thickness = 1.f;
  g.drawLine(x, y, w, y, thickness);  y += h;
  g.drawLine(x, y, w, y, thickness);  y += h;
  g.drawLine(x, y, w, y, thickness);  y += h;
  g.drawLine(x, y, w, y, thickness);  y += h;
  g.drawLine(x, y, w, y, thickness);  y += h;
  g.drawLine(x, y, w, y, thickness);  y += h;
  g.drawLine(x, y, w, y, thickness);  y += h;
  g.drawLine(x, y, w, y, thickness);  y += h;
  g.drawLine(x, y, w, y, thickness);  y += h;
  g.drawLine(x, y, w, y, thickness);  y += h;
  g.drawLine(x, y, w, y, thickness);  y += h;
  g.drawLine(x, y, w, y, thickness);

  // draw the lines between the piano-roll white keys:
  g.setColour(lineColour);
  x = 0.f;
  w = keyLength;
  h = (float) rowHeight;
  y = keyboardY + 13*rowHeight - (h+h/2);
  g.drawLine(x, y, w, y, thickness); y -= 2*h;
  g.drawLine(x, y, w, y, thickness); y -= h+h/2;
  g.drawLine(x, y, w, y, thickness); y -= h+h/2;
  g.drawLine(x, y, w, y, thickness); y -= 2*h;
  g.drawLine(x, y, w, y, thickness); y -= 2*h;
  g.drawLine(x, y, w, y, thickness); y -= h+h/2;
  g.drawLine(x, y, w, y, thickness);

  // draw the vertical lines between the steps:
  int numSteps = rosic::AcidPattern::getMaxNumSteps();
  x = keyLength;
  y = 0.f;
  h = (float) getHeight();
  for(int i=0; i<numSteps; i++)
  {
    x += columnWidth;
    if( (i+1)%4 == 0 )
      thickness = 2.f;
    else
      thickness = 1.f;
    g.drawLine(x, y, x, h, thickness);
  }

  // draw the pattern data:
  if( patternToEdit != NULL )
  {
    x            = keyLength;
    w            = columnWidth;
    h            = topLaneHeight;
    float dx     = columnWidth   / 2.f;
    float dy     = topLaneHeight / 2.f;
    int numSteps = patternToEdit->getNumSteps();
    g.setColour(handleColor);
    for(int i=0; i<patternToEdit->getMaxNumSteps(); i++)
    {
      bool slide = patternToEdit->getSlide(i) && patternToEdit->getGate((i+1)%numSteps);
      y = 0;
      if( patternToEdit->getGate(i) == true )
      {
        g.fillEllipse(x+5, y+3, w-10, h-6);
        if( slide )
          g.drawLine(x+dx, y+dy, x+3.f*dx, y+dy, 4.f);
      }

      y += topLaneHeight;
      if( patternToEdit->getAccent(i) == true )
        g.fillEllipse(x+5, y+3, w-10, h-6);

      y += topLaneHeight;
      if( patternToEdit->getSlide(i) == true )
        g.fillEllipse(x+5, y+3, w-10, h-6);

      y += topLaneHeight;
      juce::String octString = valueToStringWithSign0( patternToEdit->getOctave(i) );

      drawBitmapFontText(g, (int)(x+dx), (int)(y+dy), octString,
        &BitmapFontRoundedBoldA10D0::instance, textColour, -1, Justification::centred);

      y = keyboardY + 12*rowHeight;
      if( patternToEdit->getGate(i) == true )
      {
        y -= patternToEdit->getKey(i) * rowHeight;
        float w2 = (float) patternToEdit->getStepLength();
        if( slide )
        {
          float x2 = x + columnWidth;
          float y2 = keyboardY + 12*rowHeight - patternToEdit->getKey(i+1) * rowHeight;
          g.drawLine(x2, y+dy, x2, y2+dy, 3.f);
          w2 = 1.f;
        }
        g.fillRect(x, y, w2*columnWidth, rowHeight);
      }
      x += columnWidth;
    }
  }

  // draw enclosing rectangle:
  g.drawRect(0, 0, getWidth(), getHeight(), 2);
}

//=================================================================================================
// class AcidSequencerModuleEditor:

//-------------------------------------------------------------------------------------------------
// construction/destruction:

AcidSequencerModuleEditor::AcidSequencerModuleEditor(CriticalSection *newPlugInLock, 
  AcidSequencerAudioModule* newAcidSequencerAudioModule) : AudioModuleEditor(newAcidSequencerAudioModule)
{
  setHeadlineStyle(SUB_HEADLINE);
  setLinkPosition(INVISIBLE);

  // assign the pointer to the rosic::AcidSequencer object to be used as aduio engine:
  jassert(newAcidSequencerAudioModule != NULL ); // you must pass a valid module here
  acidSequencerModuleToEdit = newAcidSequencerAudioModule;

  isTopLevelEditor = false;

  patternEditor = new AcidPatternEditor(acidSequencerModuleToEdit->wrappedAcidSequencer);
  addAndMakeVisible(patternEditor);
  //addWidget(patternEditor);
  patternEditor->setPatternToEdit(acidSequencerModuleToEdit->wrappedAcidSequencer->getPattern(0));

  addWidget( modeLabel = new RTextField( juce::String("Mode:")) );
  modeLabel->setDescription("Chooses the sequencer mode");
  modeLabel->setDescriptionField(infoField);

  addWidget( modeBox = new RComboBox(juce::String("ModeComboBox")) );
  modeBox->assignParameter( moduleToEdit->getParameterByName("Mode") );
  modeBox->setDescription(modeLabel->getDescription());
  //modeBox->setNoBackgroundAndOutline(true);
  modeBox->setDescriptionField(infoField);

  addWidget( stepLengthSlider = new RSlider ("StepLengthSlider") );
  stepLengthSlider->assignParameter( acidSequencerModuleToEdit->getParameterByName("StepLength") );
  stepLengthSlider->setDescription(juce::String("Length of the steps in 16th notes"));
  stepLengthSlider->setStringConversionFunction(valueToString2);
  stepLengthSlider->setDescriptionField(infoField);
  stepLengthSlider->addListener(this);

  addWidget( shiftLabel = new RTextField( juce::String("Shift:")) );
  shiftLabel->setDescription("Shift the whole pattern for left or right (circularly)");
  shiftLabel->setDescriptionField(infoField);

  addWidget( shiftLeftButton = new RButton(juce::String("L")) );
  shiftLeftButton->setDescription(juce::String("Shift the whole pattern one postion to the left (circularly)"));
  shiftLeftButton->setDescriptionField(infoField);
  shiftLeftButton->setClickingTogglesState(false);
  shiftLeftButton->addRButtonListener(this);

  addWidget( shiftRightButton = new RButton(juce::String("R")) );
  shiftRightButton->setDescription(juce::String("Shift the whole pattern one postion to the right (circularly)"));
  shiftRightButton->setDescriptionField(infoField);
  shiftRightButton->setClickingTogglesState(false);
  shiftRightButton->addRButtonListener(this);

  // set up the widgets:
  updateWidgetsAccordingToState();
}

//-------------------------------------------------------------------------------------------------
// setup:



//-------------------------------------------------------------------------------------------------
// callbacks:

void AcidSequencerModuleEditor::rButtonClicked(RButton *buttonThatWasClicked)
{
  if( acidSequencerModuleToEdit==NULL || acidSequencerModuleToEdit->wrappedAcidSequencer==NULL )
    return;

  if( buttonThatWasClicked == shiftLeftButton )
    acidSequencerModuleToEdit->wrappedAcidSequencer->circularShift(-1);
  else if( buttonThatWasClicked == shiftRightButton )
    acidSequencerModuleToEdit->wrappedAcidSequencer->circularShift(+1);

  // \todo: randomization stuff....


  patternEditor->repaint();
}

void AcidSequencerModuleEditor::updateWidgetsAccordingToState()
{
  AudioModuleEditor::updateWidgetsAccordingToState();
  patternEditor->repaint();
}

void AcidSequencerModuleEditor::rSliderValueChanged(RSlider *rSliderThatHasChanged)
{
  if( rSliderThatHasChanged == stepLengthSlider )
    patternEditor->repaint();
}

void AcidSequencerModuleEditor::paint(Graphics &g)
{
  AudioModuleEditor::paint(g);
}

void AcidSequencerModuleEditor::resized()
{
  AudioModuleEditor::resized();
  int x = 0;
  int y = 0;
  int w = getWidth();
  //int h = getHeight();

  y = getPresetSectionBottom()+4;
  x = 4;
  patternEditor->setBounds(x, y, 368, 220);

  x = patternEditor->getRight();
  w = getWidth()-x;

  modeLabel->setBounds(x+4,    y+4, 40,     16);
  modeBox->setBounds(  x+40+4, y+4, w-40-8, 16);
  y = modeBox->getBottom();
  stepLengthSlider->setBounds(x+4, y+4, w-8, 16);

  y = stepLengthSlider->getBottom()+4;
  shiftLabel->setBounds(x+4, y+4, 40, 16);
  x = shiftLabel->getRight();
  shiftLeftButton->setBounds(                           x+4, y+4, 24, 16);
  shiftRightButton->setBounds(shiftLeftButton->getRight()+4, y+4, 24, 16);
}
