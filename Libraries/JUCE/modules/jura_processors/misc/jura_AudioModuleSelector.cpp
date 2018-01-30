AudioModuleSelector::AudioModuleSelector(AudioModuleFactory* factoryToUse) : RComboBox("ModuleSelector")
{
  setAudioModuleFactory(factoryToUse);


  //setSize(300, 300); // has no effect
}

void AudioModuleSelector::setAudioModuleFactory(AudioModuleFactory* newFactory)
{
  moduleFactory = newFactory;

  // todo: update contents of dropdown menu...
  
  
  // populate the tree (preliminary - to do: use moduleFactory object to figure out the available 
  // modules (and their creation functions, categories, etc.)

  RTreeViewNode *node;
  int i = 1;           //  the index is actually not used, but we need it as dummy
  popUpMenu->addTreeNodeItem(new RTreeViewNode("None",     i++));
#if JUCE_DEBUG
  popUpMenu->addTreeNodeItem(new RTreeViewNode("DebugAudioModule",    i++));
#endif

  node = new RTreeViewNode("Instruments", -1, "Instruments");
  node->addChildNode(new RTreeViewNode("AcidDevil",       i++));
  node->addChildNode(new RTreeViewNode("Straightliner",  i++));
  node->setOpen(false);
  popUpMenu->addTreeNodeItem(node);

  node = new RTreeViewNode("Effects", -1, "Effects");
  //node->addChildNode(new RTreeViewNode("Enveloper",     i++));
  node->addChildNode(new RTreeViewNode("FuncShaper",    i++));
  node->addChildNode(new RTreeViewNode("NodeShaper",    i++));
  //node->addChildNode(new RTreeViewNode("StereoDelay",   i++)); // include in Quadrifex
  //node->addChildNode(new RTreeViewNode("PitchShifter",  i++)); // include in Quadrifex
  node->addChildNode(new RTreeViewNode("EchoLab",       i++));
  node->addChildNode(new RTreeViewNode("Quadrifex",     i++));
  node->addChildNode(new RTreeViewNode("PingPongEcho",  i++));
  //node->addChildNode(new RTreeViewNode("AlgoVerb",      i++));
  node->setOpen(false);
  popUpMenu->addTreeNodeItem(node);


  node = new RTreeViewNode("Sources", -1, "Sources");
  node->addChildNode(new RTreeViewNode("EllipseOscillator", i++));
  node->addChildNode(new RTreeViewNode("Oscillator3D",   i++));
  node->addChildNode(new RTreeViewNode("RayBouncer",     i++));
  node->addChildNode(new RTreeViewNode("WaveOscillator",  i++));  // 
  node->addChildNode(new RTreeViewNode("FourOscSection",  i++));
  //node->addChildNode(new RTreeViewNode("NoiseGenerator",  i++));
  //node->addChildNode(new RTreeViewNode("SamplePlayer",    i++));
  node->setOpen(false);
  popUpMenu->addTreeNodeItem(node);


  node = new RTreeViewNode("Filters", -1, "Filters");
  node->addChildNode(new RTreeViewNode("Ladder",          i++));
  node->addChildNode(new RTreeViewNode("Equalizer",       i++));
  node->addChildNode(new RTreeViewNode("EngineersFilter", i++));
  //node->addChildNode(new RTreeViewNode("PhasorFilter",    i++));
  //node->addChildNode(new RTreeViewNode("CrossOver",       i++));
  node->setOpen(false);
  popUpMenu->addTreeNodeItem(node);

  node = new RTreeViewNode("Dynamics", -1, "Dynamics");
  //node->addChildNode(new RTreeViewNode("MultiComp",     i++));
  node->addChildNode(new RTreeViewNode("Limiter",     i++));
  popUpMenu->addTreeNodeItem(node);

  node = new RTreeViewNode("Modulators", -1, "Modulators");
  node->addChildNode(new RTreeViewNode("BreakpointModulator",  i++));
  //node->addChildNode(new RTreeViewNode("LowFrequencyOscillator",  i++));
  node->setOpen(false);
  popUpMenu->addTreeNodeItem(node);

  node = new RTreeViewNode("Analyzers", -1, "Analyzers");
  node->addChildNode(new RTreeViewNode("Scope",    i++));
  //node->addChildNode(new RTreeViewNode("PhaseScope2",   i++));
  node->addChildNode(new RTreeViewNode("MultiAnalyzer", i++));
  //node->addChildNode(new RTreeViewNode("TrackMeter",    i++));
  node->addChildNode(new RTreeViewNode("MidiMonitor",   i++));
  node->setOpen(false);
  popUpMenu->addTreeNodeItem(node);

  bool showUnfinishedModules = true;
  if(showUnfinishedModules)
  {
    node = new RTreeViewNode("UnderConstruction", -1, "UnderConstruction");
#ifdef _MSC_VER
    node->addChildNode(new RTreeViewNode("Liberty", i++));     // not yet available on gcc
#endif
    node->addChildNode(new RTreeViewNode("NewSynth",  i++));
    //node->addChildNode(new RTreeViewNode("MagicCarpet",    i++));
    node->addChildNode(new RTreeViewNode("SimpleSampler",  i++));
    //node->addChildNode(new RTreeViewNode("KeyShot",        i++));
    //node->addChildNode(new RTreeViewNode("Quadriga",       i++));
    //node->addChildNode(new RTreeViewNode("Workhorse",      i++));
    node->setOpen(false);
    popUpMenu->addTreeNodeItem(node);
  }

}

// The current release version includes:
// Instruments:
//   AcidDevil
//   Straightliner
// Effects:
//   FuncShaper
//   EchoLab
// Filters:
//   Equalizer
//   EngineersFilter
// Dynamics:
//   Limiter
// Analyzers:
//   Scope
//   MultiAnalyzer
//   MidiMonitor

void AudioModuleSelector::drawHighlighted(bool shouldBeHighlighted)
{
  highlighted = shouldBeHighlighted;
  repaint();
}

void AudioModuleSelector::paint(Graphics& g)
{
  if(!highlighted)
    RComboBox::paint(g);
  else
  {
    g.fillAll(getHandleColour()); // only difference to baseclass version - maybe refactor
    g.setColour(getOutlineColour());
    g.drawRect(0, 0, getWidth(), getHeight(), 2);
    int x = 4;
    int y = getHeight()/2 - font->getFontAscent()/2;
    drawBitmapFontText(g, x, y, getSelectedItemText(), font, getTextColour());
  }
}
