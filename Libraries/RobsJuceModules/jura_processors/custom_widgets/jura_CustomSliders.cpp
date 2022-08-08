
TuningSlider::TuningSlider(const juce::String& componentName) : rsAutomatableSlider()
{
  //jassertfalse; // it needs to be checked carefully, if the popup menu works as intended
}

void TuningSlider::rPopUpMenuChanged(RPopUpMenu* menuThatHasChanged)
{
  RTreeViewNode *selectedItem = rightClickPopUp->getSelectedItem();
  if( selectedItem == NULL )
    return;

  switch( selectedItem->getNodeIdentifier() )
  {
  case OCTAVE_UP:   applySemitoneShiftToValue( 12.0); break;
  case OCTAVE_DOWN: applySemitoneShiftToValue(-12.0); break;
  default:          rsAutomatableSlider::rPopUpMenuChanged(menuThatHasChanged);
  }
}

void TuningSlider::addPopUpMenuItems()
{
  //jassertfalse; // it needs to be checked carefully, if the popup menu works as intended

  addPopUpEnterValueItem();
  addPopUpOctaveUpDownItems();
  addPopUpDefaultValueItems();
  addPopUpMidiItems();
}

void TuningSlider::addPopUpOctaveUpDownItems()
{
  rightClickPopUp->addItem(OCTAVE_UP,   "Octave Up");
  rightClickPopUp->addItem(OCTAVE_DOWN, "Octave Down");
}

void TuningSlider::applySemitoneShiftToValue(double numSemitones)
{
  if( assignedParameter != NULL ) 
    assignedParameter->setValue(assignedParameter->getValue()+numSemitones, true, true);
  else
    setValue(getValue()+numSemitones, true);
}
