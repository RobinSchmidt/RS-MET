#ifndef jura_CustomSliders_h
#define jura_CustomSliders_h

/** This class overrides RSlider to cater for some requirements that are specific to sliders that 
set up a tuning parameter.  */

class TuningSlider : public AutomatableSlider
{

public:

  //---------------------------------------------------------------------------------------------
  // construction/destruction:

  /** Constructor. */
  TuningSlider(const juce::String& componentName);


protected:

  /** Enumeration of the additional popup item identifiers specific to the TuningSlider class. */
  enum tuningItemIdentifiers
  {
    OCTAVE_UP = META_DETACH+1,
    OCTAVE_DOWN
  };

  // overrides:
  virtual void rPopUpMenuChanged(RPopUpMenu* menuThatHasChanged);
  virtual void addPopUpMenuItems();

  // new member functions:
  virtual void addPopUpOctaveUpDownItems();
  virtual void applySemitoneShiftToValue(double numSemitones);

  juce_UseDebuggingNewOperator;
};

#endif  
