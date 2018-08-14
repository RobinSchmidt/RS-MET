EllipseOscillatorAudioModule::EllipseOscillatorAudioModule(CriticalSection *lockToUse, 
  MetaParameterManager* metaManagerToUse, ModulationManager* modManagerToUse)
  : AudioModuleWithMidiIn(lockToUse, metaManagerToUse, modManagerToUse)
{
  ScopedLock scopedLock(*lock);
  setModuleTypeName("EllipseOscillator");
  createParameters();
}

void EllipseOscillatorAudioModule::createParameters()
{
  ScopedLock scopedLock(*lock);

  typedef rosic::rsEllipseOscillator EO;
  EO* eo = &oscCore;

  typedef ModulatableParameter Param;
  Param* p;

  p = new Param("Amplitude", -1.0, +1.0, 1.0, Parameter::LINEAR);
  addObservedParameter(p);
  p->setValueChangeCallback<EO>(eo, &EO::setAmplitude);

  p = new Param("Tune", -60.0, +60.0, 0.0, Parameter::LINEAR);
  addObservedParameter(p);
  p->setValueChangeCallback<EO>(eo, &EO::setDetune);

  //p = new Param("A", -1.0, 1.0, 0.0, Parameter::LINEAR);
  p = new Param("LowHigh", -1.0, 1.0, 0.0, Parameter::LINEAR);
  addObservedParameter(p);
  p->setValueChangeCallback<EO>(eo, &EO::setA);
  // determines, how the speed varies over the full cycle: 0 constant speed, 
  // +-1: skips half a cycle (i think)

  p = new Param("InOut", -10.0, 10.0, 1.0, Parameter::LINEAR);
  addObservedParameter(p);
  p->setValueChangeCallback<EO>(eo, &EO::setC);
  // determines, how speed varies over a half-cycles?

  p = new Param("Phase", -180.0, +180.0, 0.0, Parameter::LINEAR);
  addObservedParameter(p);
  p->setValueChangeCallback<EO>(eo, &EO::setRotationDegrees);
  // adjusts horizontal vs vertical extent of ellipse - at 45�, its a circle


  // A is Upper/Lower Bias, B is Left/Right Bias, C is Upper/Lower Pinch

  // Renormalize, Mix, FreqScaleY, FreqShiftY

  // we really need some rescaling of LowHigh and InOut
}

void EllipseOscillatorAudioModule::processBlock(double **inOutBuffer, int numChannels, int numSamples)
{
  for(int n = 0; n < numSamples; n++)
    oscCore.getSampleFrameStereo(&inOutBuffer[0][n], &inOutBuffer[1][n]);
}
void EllipseOscillatorAudioModule::processStereoFrame(double *left, double *right)
{
  oscCore.getSampleFrameStereo(left, right);
}

void EllipseOscillatorAudioModule::setSampleRate(double newSampleRate)
{
  oscCore.setSampleRate(newSampleRate);
}

void EllipseOscillatorAudioModule::reset()
{
  oscCore.reset();
}

void EllipseOscillatorAudioModule::noteOn(int noteNumber, int velocity)
{
  oscCore.setFrequency(pitchToFreq(noteNumber)); // preliminary - use tuning table
  //oscCore.reset();
}

//=================================================================================================

TriSawOscModule::TriSawOscModule(CriticalSection *lockToUse, 
  MetaParameterManager* metaManagerToUse, ModulationManager* modManagerToUse)
  : AudioModuleWithMidiIn(lockToUse, metaManagerToUse, modManagerToUse)
{
  ScopedLock scopedLock(*lock);
  setModuleTypeName("TriSawOscillator");
  createParameters();
}

void TriSawOscModule::createParameters()
{
  ScopedLock scopedLock(*lock);

  typedef RAPT::rsTriSawOscillator<double> TSO;
  TSO* tso = &oscCore;

  typedef ModulatableParameter Param;
  Param* p;

  std::vector<double> defaultValues; // used for figuring out map-points during development

  p = new Param("Asymmetry", -1.0, 1.0, 0.0, Parameter::LINEAR);
  //p->setMapper(new rsParameterMapperRationalBipolar(-1, +1, 0.80)); // test
  addObservedParameter(p);
  p->setValueChangeCallback<TSO>(tso, &TSO::setAsymmetry);

  /*
  // for later, to replace AttackBending and DecayBending parameters - but that doesn't work yet:
  p = new Param("Bending", -1.0, 1.0, 0.0, Parameter::LINEAR);
  //p->setMapper(new rsParameterMapperRationalBipolar(-1, +1, 0.9));
  addObservedParameter(p);
  p->setValueChangeCallback<TriSawOscModule>(this, &TriSawOscModule::setBend);

  p = new Param("BendAsym", -1.0, 1.0, 0.0, Parameter::LINEAR);
  //p->setMapper(new rsParameterMapperRationalBipolar(-1, +1, 0.9));
  addObservedParameter(p);
  p->setValueChangeCallback<TriSawOscModule>(this, &TriSawOscModule::setBendAsym);
  */
  
  p = new Param("AttackBending", -1.0, 1.0, 0.0, Parameter::LINEAR);
  //p->setMapper(new rsParameterMapperRationalBipolar(-1, +1, 0.9));
  //defaultValues.clear();
  //defaultValues.push_back(0.0);    // x = 0.0
  //defaultValues.push_back(0.5);    // x = 0.25
  //defaultValues.push_back(0.85);   // x = 0.5
  //defaultValues.push_back(0.975);  // x = 0.75
  //defaultValues.push_back(1.0);    // x = 1.0
  //p->setDefaultValues(defaultValues);
  addObservedParameter(p);
  p->setValueChangeCallback<TSO>(tso, &TSO::setAttackBending);

  p = new Param("DecayBending", -1.0, 1.0, 0.0, Parameter::LINEAR);
  //p->setMapper(new rsParameterMapperRationalBipolar(-1, +1, 0.9));
  addObservedParameter(p);
  p->setValueChangeCallback<TSO>(tso, &TSO::setDecayBending);


  p = new Param("AttackSigmoid", -1.0, 1.0, 0.0, Parameter::LINEAR);
  addObservedParameter(p);
  p->setValueChangeCallback<TSO>(tso, &TSO::setAttackSigmoid);

  p = new Param("DecaySigmoid", -1.0, 1.0, 0.0, Parameter::LINEAR);
  addObservedParameter(p);
  p->setValueChangeCallback<TSO>(tso, &TSO::setDecaySigmoid);

  // maybe later allow the (raw) sigmoid parameter go from -2..+1 because at -2, we get a flat 
  // section in the middle. i think a mapped (not raw anymore) sigmoid parameter should use a 
  // mapping curve that goes through (-1,-2),(0,0),(1,1) -> find a function through these points 
  // that feels natural


  // add AttackDrive/DecayDrive

  // for the envelope, add floor/ceil, set it up in terms of AttackTime/DecayTime

  // asymmetry and bedning need finer resolution at the ends / less in the center
  // maybe try tanh mapper or rational s-curve
  // ok, the mappers are currently commented - maybe thei shape paremeters should themselves be
  // availabe as user parameters - at least during development to figure out best (most musical) 
  // settings for them, and later maybe hard-coded - maybe to figure out the best curve, use a 
  // sinusoidal LFO and modulate the parameter with it and choose a setting in which the LFO sounds
  // "most sinusoidal" ...whatever that means


  // instead of AttackBending (ab) and DecayBending (db) provide
  // Bending = b = (ab+db)/2 and BendAsym = ba = (ab-db)/2
  // and similar for the two sigmoid parameters
}

void TriSawOscModule::processStereoFrame(double *left, double *right)
{
  *left = *right = oscCore.getSample();
}

void TriSawOscModule::setSampleRate(double newSampleRate)
{
  sampleRate = newSampleRate;
  oscCore.setPhaseIncrement(freq/sampleRate);
}

void TriSawOscModule::reset()
{
  oscCore.reset();
}

void TriSawOscModule::noteOn(int noteNumber, int velocity)
{
  freq = pitchToFreq(noteNumber);
  oscCore.setPhaseIncrement(freq/sampleRate); // maybe times 2 because it ranegs from -1..+1?
}

void TriSawOscModule::setBend(double newBend)
{
  bend = newBend;
  updateBending();
}

void TriSawOscModule::setBendAsym(double newAsym)
{
  bendAsym = newAsym;
  updateBending();
}

void TriSawOscModule::updateBending()
{
  // this does not work yet - when asym is too much, it messes up totally
  //oscCore.setAttackBending(bend + bendAsym);
  //oscCore.setDecayBending( bend - bendAsym);

  // Bend looks good but when asym ist introduced, it gets messed up...maybe i need some 
  // renormalization constant, i.e. scale attack/decay bend values according to the bend-asym?

  //double s = 1. / (1. - fabs(bendAsym)); // ad-hoc ...why....derive it
  //oscCore.setAttackBending(bend + s*bendAsym);
  //oscCore.setDecayBending( bend - s*bendAsym);

  // hmm..well when bend b=1 and bendAsym a=1 then
  // ba = b+a = 2, bd = b-a = 0 but we actually want ba = 1, bd = 0

  // that's not yet very nice:
  oscCore.setAttackBending(rosic::clip(bend + bendAsym, -1., +1.) );
  oscCore.setDecayBending( rosic::clip(bend - bendAsym, -1., +1.) );
  // but maybe it doesn't get any better

  //// this?
  //double s = 1 / (1 + abs(bendAsym));  // -1..1 -> 1..0..1
  //double p = 0.5 * (bendAsym + 1);     // -1..1 -> 0..1
  //double sb2 = 2*s*bend;
  //oscCore.setAttackBending(sb2 *   p  );
  //oscCore.setDecayBending( sb2 * (1-p));
  //// has a nice behavior but makes certain settings inacessible (for example 
  //// attack bend = 0.9, decay bend = -0.9

  //double s = 1 / sqrt(1 + bend*bend + bendAsym*bendAsym);
  //oscCore.setAttackBending(s * (bend+bendAsym));
  //oscCore.setDecayBending( s * (bend-bendAsym));



  // ...under construction - 
  // idea: 
  // -collect data for perceptually equidistant values for attackBend
  // -for example, finde the y value for x = 0, 0.25, 0.5, 0.75, 1.0, like this:
  //  (0,0),(0.25,0.5),(0.5,0.8),(0.75,0.9),(1,1) (values plucked out of thin air)
  // -to these values, fit a function of a given form, for example
  //  f(x) = (a1*x + a3*x^3) / (1 + b2*x^2 + b4*x^4)  -> odd symmetry
  //  https://www.desmos.com/calculator/6n6va65ige maybe make it 3-parameteric by imposing
  //  f(1) = 1  ...f(-1) = -1 is no additional constraint since the function is already symmetric 
  //  by construction, also we already have f(0) = 0 by construction 
  //  -constraint: 1 = (a+b)/(1+c+d) -> solve for d, leaving only a,b,c as parameters
  //   https://www.desmos.com/calculator/adns6zwdm7
  //  -maybe another constraint could be monotonicity in 0..1 -> f'(x) > 0 in 0..1
  //  -maybe try to define it via function values at 1/4, 1/2, 3/4
  // 
  // -write a pade approximation function for that (useful for the toolkit anyway) or use
  //  numpy/scipy curve fitting
  // -maybe this should allow also constraints for fixing some coeffs at given values and enforce
  //  to go through specific points...maybe pass arrays for the fixed coeffs where NaN encodes a 
  //  free coeff and an array of weights for data points where inf encodes to enforce passing
  //  through the value
  // -this procedure can be generally used to find good parameter mapping functions

  // -maybe it makes not much sense to figure out a each one-dimensional parameter mapping such as
  //  the one for attackBend
  // -for example, if only the raw attack-bend parameter is used, the mapping curve would contain
  //  (0,0),(0.25,5),(0.5,0.85),(0.75,0.975),(1,1) where in the first half of the curve, the 
  //  spectrum would fill up the even harmonics and in the second half raise the overtones, making
  //  it more sharp - but if we would use a symmetric bend parameter that sets attack and decay 
  //  bend at the same time, no such spectrum-filling first half would occur - no qualitative 
  //  change of the mapping function midway the parameter range
  // -maybe use a neural network to map from N user-parameters to M algo-parameters?
  //  -how to acquire the data?
  //  -how to impose constraints such as that the raw attack-bend and decay-bend (as network 
  //   outputs) should behave symmetric (in some sense)?
}