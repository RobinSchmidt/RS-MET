BlepOscArrayModule::BlepOscArrayModule(CriticalSection *lockToUse)
  : AudioModuleWithMidiIn(lockToUse) /*, oscArrayCore(&ratioGenerator)*/
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


  p = new Param("Density", 1.0, 32.0, 7.0, Parameter::INTEGER);
  p->setValueChangeCallback<OA>(oa, &OA::setDensity);
  addObservedParameter(p);

  p = new Param("Detune", 0.0, 100.0, 5.0, Parameter::LINEAR);
  p->setValueChangeCallback<OA>(oa, &OA::setDetunePercent);
  addObservedParameter(p);
  // maybe rename to spread - the default value of 5 seems ok - but we may want to change that
  // later - in order to not break presets, it should be saved, even when in default setting - but
  // Parameter only allows to switch save/recall on/off globally for all parameters - maybe allow
  // to switch this locally


  // uses the new ChoiceParameter class - needs testing:
  rsChoiceParameter* cp;
  cp = new rsChoiceParameter("Distribution");
  cp->setValueChangeCallback<OA>(oa, &OA::setFrequencyDistribution);
  typedef RAPT::rsRatioGenerator<double>::RatioKind RK;
  cp->addStringValue("Range Split Odd",        (int)RK::rangeSplitOdd);     // 4
  cp->addStringValue("Range Split Even",       (int)RK::rangeSplitEven);    // 5
  cp->addStringValue("Range Split Skewed",     (int)RK::rangeSplitSkewed);  // 3
  cp->addStringValue("Prime Square Root",      (int)RK::primeSqrt);         // 1
  cp->addStringValue("Prime Square Root Diff", (int)RK::primeSqrtDiff );    // 2
  cp->addStringValue("Metallic",               (int)RK::metallic);          // 0
  cp->setDefaultValue((int)RK::primeSqrtDiff, true);
  // try to get rid of the explicit conversions to int here by introducing a template function
  // in ChoiceParameter that does this - but only if this doesn't lead to code bloat - figure this
  // out first - if it does bloat the binary, keep it as is - tried it - got a linker error
  addObservedParameter(cp);
  // Prime Square Root Diff is most neutral in terms of the transient - even when initial coherence
  // is 0, there is some sort of transient with the other distributions



  p = new Param("DistParam", 0.0, 1.0, 0.5, Parameter::LINEAR);
  p->setValueChangeCallback<OA>(oa, &OA::setDistributionParameter);
  addObservedParameter(p);

  p = new Param("InitialPhaseCoherence", 0.0, 1.0, 0.0, Parameter::LINEAR);
  p->setDisplayName("Coherence");
  p->setValueChangeCallback<OA>(oa, &OA::setIntialPhaseCoherence);
  addObservedParameter(p);
  // maybe call the display-parameter "Transient"...hmm - but that will be confusing for users who
  // want to look at and/or edit the xml files




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