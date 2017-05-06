
EffectSelectionPopup::EffectSelectionPopup(Component *componentToAttachTo) 
: RPopUpMenu(componentToAttachTo)
{
  // populate the menu:

  addItem(rosic::Quadrifex::MUTE,   juce::String("Mute"));
  addItem(rosic::Quadrifex::BYPASS, juce::String("Bypass"));

  RTreeViewNode *tmpNode1;

  tmpNode1 = new RTreeViewNode("Delays", -1, "Delayline based effects", true, false, false, NULL);
  tmpNode1->addChildNode(new RTreeViewNode("SimpleDelay",  rosic::Quadrifex::SIMPLE_DELAY));
  tmpNode1->addChildNode(new RTreeViewNode("PingPongEcho", rosic::Quadrifex::PINGPONG_ECHO));
  tmpNode1->addChildNode(new RTreeViewNode("Reverb",       rosic::Quadrifex::REVERB));
  addTreeNodeItem(tmpNode1);

  tmpNode1 = new RTreeViewNode("Filters", -1, "Filters", true, false, false, NULL);
  tmpNode1->addChildNode(new RTreeViewNode("TwoPoleFilter",     rosic::Quadrifex::TWO_POLE_FILTER));
  tmpNode1->addChildNode(new RTreeViewNode("DualTwoPoleFilter", rosic::Quadrifex::DUAL_TWO_POLE_FILTER));
  tmpNode1->addChildNode(new RTreeViewNode("CombBank",          rosic::Quadrifex::COMB_BANK));
  tmpNode1->addChildNode(new RTreeViewNode("CombResonator",     rosic::Quadrifex::COMB_RESONATOR));
  tmpNode1->addChildNode(new RTreeViewNode("Equalizer",         rosic::Quadrifex::EQUALIZER));
  tmpNode1->addChildNode(new RTreeViewNode("FourPoleFilter",    rosic::Quadrifex::FOUR_POLE_FILTER));
  tmpNode1->addChildNode(new RTreeViewNode("LadderFilter",      rosic::Quadrifex::LADDER_FILTER));
  tmpNode1->addChildNode(new RTreeViewNode("SlopeFilter",       rosic::Quadrifex::SLOPE_FILTER));
  addTreeNodeItem(tmpNode1);

  tmpNode1 = new RTreeViewNode("Modulation", -1, "Modulation", true, false, false, NULL);
  tmpNode1->addChildNode(new RTreeViewNode("Chorus",  rosic::Quadrifex::CHORUS));
  tmpNode1->addChildNode(new RTreeViewNode("Flanger", rosic::Quadrifex::FLANGER));
  tmpNode1->addChildNode(new RTreeViewNode("Phaser",  rosic::Quadrifex::PHASER));
  tmpNode1->addChildNode(new RTreeViewNode("Tremolo", rosic::Quadrifex::TREMOLO));
  tmpNode1->addChildNode(new RTreeViewNode("Vibrato", rosic::Quadrifex::VIBRATO));
  tmpNode1->addChildNode(new RTreeViewNode("WahWah",  rosic::Quadrifex::WAH_WAH));
  addTreeNodeItem(tmpNode1);

  tmpNode1 = new RTreeViewNode("Dynamics", -1, "Dynamics", true, false, false, NULL);
  tmpNode1->addChildNode(new RTreeViewNode("Compressor", rosic::Quadrifex::COMPRESSOR));
  tmpNode1->addChildNode(new RTreeViewNode("Expander",   rosic::Quadrifex::EXPANDER));
  tmpNode1->addChildNode(new RTreeViewNode("NoiseGate",  rosic::Quadrifex::NOISE_GATE));
  tmpNode1->addChildNode(new RTreeViewNode("Limiter",    rosic::Quadrifex::LIMITER));
  addTreeNodeItem(tmpNode1);

  tmpNode1 = new RTreeViewNode("Distortion", -1, "Distortion", true, false, false, NULL);
  tmpNode1->addChildNode(new RTreeViewNode("WaveShaper",       rosic::Quadrifex::WAVESHAPER));
  tmpNode1->addChildNode(new RTreeViewNode("BitCrusher",       rosic::Quadrifex::BIT_CRUSHER));
  tmpNode1->addChildNode(new RTreeViewNode("Harmonics",        rosic::Quadrifex::HARMONICS));
  tmpNode1->addChildNode(new RTreeViewNode("ModulatedAllpass", rosic::Quadrifex::MODULATED_ALLPASS));
  addTreeNodeItem(tmpNode1);

  tmpNode1 = new RTreeViewNode("Stereo Tools", -1, "Stereo Tools", true, false, false, NULL);
  tmpNode1->addChildNode(new RTreeViewNode("StereoPan",       rosic::Quadrifex::STEREO_PAN));
  tmpNode1->addChildNode(new RTreeViewNode("StereoWidth",     rosic::Quadrifex::STEREO_WIDTH));
  tmpNode1->addChildNode(new RTreeViewNode("PhaseStereoizer", rosic::Quadrifex::PHASE_STEREOIZER));
  tmpNode1->addChildNode(new RTreeViewNode("CombStereoizer",  rosic::Quadrifex::COMB_STEREOIZER));
  addTreeNodeItem(tmpNode1);

  // insert Spectral: FormantShifter, etc.

  tmpNode1 = new RTreeViewNode("Generators", -1, "Generators", true, false, false, NULL);
  tmpNode1->addChildNode(new RTreeViewNode("SineOscillator", rosic::Quadrifex::SINE_OSCILLATOR));
  tmpNode1->addChildNode(new RTreeViewNode("Noisifier",      rosic::Quadrifex::NOISIFIER));
  addTreeNodeItem(tmpNode1);

  tmpNode1 = new RTreeViewNode("Miscelanneous", -1, "Miscelanneous", true, false, false, NULL);
  tmpNode1->addChildNode(new RTreeViewNode("FrequencyShifter",        rosic::Quadrifex::FREQUENCY_SHIFTER));
  tmpNode1->addChildNode(new RTreeViewNode("PitchShifter",            rosic::Quadrifex::PITCH_SHIFTER));
  tmpNode1->addChildNode(new RTreeViewNode("RingModulator",           rosic::Quadrifex::RINGMODULATOR));
  tmpNode1->addChildNode(new RTreeViewNode("SingleSidebandModulator", rosic::Quadrifex::SSB_MODULATOR));
  tmpNode1->addChildNode(new RTreeViewNode("SlewRateLimiter",         rosic::Quadrifex::SLEWRATE_LIMITER));
  addTreeNodeItem(tmpNode1);
}

/*
void EffectSelectionPopup::createRightClickPopupMenu(PopupMenu *&mainMenu, PopupMenu *&defaultValueSubMenu, 
                                        int &defaultValueIndicesMin, int &defaultValueIndicesMax)
{
  RSlider::createRightClickPopupMenu(mainMenu, defaultValueSubMenu, defaultValueIndicesMin, 
    defaultValueIndicesMax);

  // add some items specific to this kind of slider:  
  int index = defaultValueIndicesMax+1;
  mainMenu->addItem(index, juce::String(T("octave up")) );
  index++;
  mainMenu->addItem(index, juce::String(T("octave down")) );
}

void EffectSelectionPopup::handleRightClickPopupMenuResult(int result, int defaultValueIndicesMin, 
                                                   int defaultValueIndicesMax)
{
  if( assignedParameter == NULL ) 
    return;

  if( result > defaultValueIndicesMax )
  {
    result -= defaultValueIndicesMax;
    if( result == 1 )
      assignedParameter->setValue(assignedParameter->getValue()+12.0, true);
    else if( result == 2 )
      assignedParameter->setValue(assignedParameter->getValue()-12.0, true);
  }
  else
    RSlider::handleRightClickPopupMenuResult(result,defaultValueIndicesMin,defaultValueIndicesMax);
}
*/

/*
void RSlider::openRightClickPopupMenu()
{
  // create the menu and its submenu:
  PopupMenu* menu                 = NULL;
  PopupMenu* defaultValuesSubMenu = NULL;
  int defaultValueIndicesMin, defaultValueIndicesMax;
  createRightClickPopupMenu(menu, defaultValuesSubMenu, defaultValueIndicesMin, 
    defaultValueIndicesMax);

  int result = menu->show();

  // we have retrieved the result - we don't need the menus anymore:
  if( defaultValuesSubMenu != NULL )
    delete defaultValuesSubMenu;
  if( menu != NULL )
    delete menu;

  handleRightClickPopupMenuResult(result, defaultValueIndicesMin, defaultValueIndicesMax);
}
*/