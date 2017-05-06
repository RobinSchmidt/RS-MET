
//=========================================================================================================================================
// class FourPoleFilterModeComboBox:

/*
void FourPoleFilterModeComboBox::openPopUp()
{
  if (! menuActive)
  {
    const int currentId = getSelectedItemIndex() + 1;
    ComponentDeletionWatcher deletionWatcher (this);

    int id = 1;

    PopupMenu menu;
    menu.setLookAndFeel(&getLookAndFeel());
    menu.addItem(id, juce::String(T("Bypass")), true, currentId == id );
    id++;

    PopupMenu lowpassMenu;
    lowpassMenu.addItem(id, juce::String(T("6 dB/oct")), true, currentId == id );
    id++;
    lowpassMenu.addItem(id, juce::String(T("12 dB/oct")), true, currentId == id );
    id++;
    lowpassMenu.addItem(id, juce::String(T("18 dB/oct")), true, currentId == id );
    id++;
    lowpassMenu.addItem(id, juce::String(T("24 dB/oct")), true, currentId == id );
    id++;
    menu.addSubMenu(juce::String(T("Lowpass")), lowpassMenu);

    PopupMenu highpassMenu;
    highpassMenu.addItem(id, juce::String(T("6 dB/oct")), true, currentId == id );
    id++;
    highpassMenu.addItem(id, juce::String(T("12 dB/oct")), true, currentId == id );
    id++;
    highpassMenu.addItem(id, juce::String(T("18 dB/oct")), true, currentId == id );
    id++;
    highpassMenu.addItem(id, juce::String(T("24 dB/oct")), true, currentId == id );
    id++;
    menu.addSubMenu(juce::String(T("Highpass")), highpassMenu);

    PopupMenu bandpassMenu;
    bandpassMenu.addItem(id, juce::String(T("6+6 dB/oct")), true, currentId == id );
    id++;
    bandpassMenu.addItem(id, juce::String(T("12+12 dB/oct")), true, currentId == id );
    id++;
    bandpassMenu.addItem(id, juce::String(T("6+12 dB/oct")), true, currentId == id );
    id++;
    bandpassMenu.addItem(id, juce::String(T("12+6 dB/oct")), true, currentId == id );
    id++;
    bandpassMenu.addItem(id, juce::String(T("6+18 dB/oct")), true, currentId == id );
    id++;
    bandpassMenu.addItem(id, juce::String(T("18+6 dB/oct")), true, currentId == id );
    id++;
    menu.addSubMenu(juce::String(T("Bandpass")), bandpassMenu);

    PopupMenu notchMenu;
    notchMenu.addItem(id, juce::String(T("2nd Order")), true, currentId == id );
    id++;
    notchMenu.addItem(id, juce::String(T("4th Order")), true, currentId == id );
    id++;
    notchMenu.addItem(id, juce::String(T("Two Notches")), true, currentId == id );
    id++;
    menu.addSubMenu(juce::String(T("Notch")), notchMenu);

    PopupMenu peakMenu;
    peakMenu.addItem(id, juce::String(T("2nd Order")), true, currentId == id );
    id++;
    peakMenu.addItem(id, juce::String(T("4th Order")), true, currentId == id );
    id++;
    peakMenu.addItem(id, juce::String(T("Flat Top")), true, currentId == id );
    id++;
    peakMenu.addItem(id, juce::String(T("Two Peaks")), true, currentId == id );
    id++;
    menu.addSubMenu(juce::String(T("Peak/Dip")), peakMenu);

    PopupMenu lowShelvingMenu;
    lowShelvingMenu.addItem(id, juce::String(T("1st Order")), true, currentId == id );
    id++;
    lowShelvingMenu.addItem(id, juce::String(T("2nd Order")), true, currentId == id );
    id++;
    lowShelvingMenu.addItem(id, juce::String(T("3rd Order")), true, currentId == id );
    id++;
    lowShelvingMenu.addItem(id, juce::String(T("4th Order")), true, currentId == id );
    id++;
    menu.addSubMenu(juce::String(T("LowShelving")), lowShelvingMenu);

    PopupMenu highShelvingMenu;
    highShelvingMenu.addItem(id, juce::String(T("1st Order")), true, currentId == id );
    id++;
    highShelvingMenu.addItem(id, juce::String(T("2nd Order")), true, currentId == id );
    id++;
    highShelvingMenu.addItem(id, juce::String(T("3rd Order")), true, currentId == id );
    id++;
    highShelvingMenu.addItem(id, juce::String(T("4th Order")), true, currentId == id );
    id++;
    menu.addSubMenu(juce::String(T("HighShelving")), highShelvingMenu);

    // allpass, special (notch-peak-notch, punctured low/high/bandpass etc....), ...

    const int itemHeight = jlimit (12, 24, getHeight());
    menuActive = true;

    const int resultIndex = menu.showAt(this, currentId, getWidth(), 1, itemHeight) - 1;

    if (deletionWatcher.hasBeenDeleted())
      return;

    menuActive = false;

    if( resultIndex > -1 )
      setSelectedItemIndex(resultIndex);
  }
}
*/

//=========================================================================================================================================
// class WaveformComboBox:

WaveformComboBox::WaveformComboBox(const juce::String& componentName) : RComboBox(componentName)
{
  addItem(rosic::SINE,                   "Sine"    );
  addItem(rosic::TRIANGLE,               "Triangle");
  addItem(rosic::SQUARE,                 "Square"  );
  addItem(rosic::SAW,                    "Saw"     );
  addItem(rosic::NUM_STANDARD_WAVEFORMS, "Load..." );
}

/*
void WaveformComboBox::openPopUp()
{
  RComboBox::openPopUp();
  // handle the 'Load...' option
}
*/