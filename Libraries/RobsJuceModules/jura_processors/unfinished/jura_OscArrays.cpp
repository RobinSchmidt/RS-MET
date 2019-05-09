BlepOscArrayModule::BlepOscArrayModule(CriticalSection *lockToUse)
  : AudioModule(lockToUse), oscArrayCore(&ratioGenerator)
{
  ScopedLock scopedLock(*lock);
  setModuleTypeName("BlepOscArray");
  setModuleName("OscArray");
  createParameters();
}
/*
AudioModuleEditor* BlepOscArrayModule::createEditor(int type)
{
  return new jura::BlepOscArrayEditor(this);
}
*/
void BlepOscArrayModule::createParameters()
{
  typedef Parameter Param;
  Param* p;

  typedef rosic::rsOscArrayPolyBlep1 OA;
  OA* oa = &oscArrayCore;


  p = new Param("Density", 1.0, 7.0, 1.0, Parameter::INTEGER);
  p->setValueChangeCallback<OA>(oa, &OA::setNumOscillators);
  addObservedParameter(p);

  p = new Param("Detune", 0.0, 100.0, 0.0, Parameter::LINEAR);
  p->setValueChangeCallback<OA>(oa, &OA::setDetune);
  addObservedParameter(p);
  // maybe rename to spread



  //p = new ChoiceParameter("Distribution", 1.0, 7.0, 1.0, Parameter::STRING);
  ChoiceParameter* cp;
  cp = new ChoiceParameter("Distribution");
  //p->setValueChangeCallback<OA>(oa, &OA::setFrequencyDistribution);

  // we would like to write things like
  typedef RAPT::rsRatioGenerator<double>::RatioKind RK;
  cp->addStringValue("Range Split Odd",    (int)RK::rangeSplitOdd);
  cp->addStringValue("Range Split Even",   (int)RK::rangeSplitEven);
  cp->addStringValue("Range Split Skewed", (int)RK::rangeSplitOdd);
  // is this possible? the enum-class does not support implicit conversion to strings - on which
  // the string parameter handling (unfortunately) relies - in Parameter we would somehow have to 
  // make an association between a string and a value from an enum class that is not yet known in
  // Parameter and should be generic - i think, we should deprecate the way, Strings are handled
  // ...maybe we should have baseclass Parameter with subclasses NumericParameter and 
  // ChoiceParameter
  addObservedParameter(cp);

  // ...
}


//=================================================================================================

BlepOscArrayEditor::BlepOscArrayEditor(BlepOscArrayModule* oscArrayToEdit)
  : AudioModuleEditor(oscArrayToEdit->lock), oscArrayModule(oscArrayToEdit)
{
  ScopedLock scopedLock(*lock);

  // ...
}

void BlepOscArrayEditor::resized()
{
  ScopedLock scopedLock(*lock);
  AudioModuleEditor::resized();

  int x = 0;
  int y = getPresetSectionBottom()+4;
  int w = getWidth();
  int h = getHeight() - y;

  // ...
}