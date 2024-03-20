SineOscAudioModule::SineOscAudioModule(CriticalSection* lockToUse,
  MetaParameterManager* metaManagerToUse, ModulationManager* modManagerToUse)
  : AudioModuleWithMidiIn(lockToUse, metaManagerToUse, modManagerToUse)
{
  ScopedLock scopedLock(*lock);
  setModuleTypeName("SineOscillator");
  createParameters();
}

void SineOscAudioModule::createParameters()
{
  ScopedLock scopedLock(*lock);

  typedef SineOscCore SO;
  SO* so = &core;

  typedef ModulatableParameter Param;
  Param* p;

  p = new Param("Amplitude", -1.0, +1.0, 1.0, Parameter::LINEAR);
  addObservedParameter(p);
  p->setValueChangeCallback<SO>(so, &SO::setAmplitude);

  p = new Param("Detune", -60.0, +60.0, 0.0, Parameter::LINEAR);
  addObservedParameter(p);
  p->setValueChangeCallback<SO>(so, &SO::setDetune);
  // trying to modulate the freq doesn't work because in each call to setFrequency, the osc resets
  // its phase? ...maybe use a more basic osc class such as RAPT::rsSineIterator - but maybe a 
  // phasor-based implementation is better for freq modulation anyway - for the time being, just
  // don't modulate detune
  
  // maybe  make a start-phase parameter - or better: just a phase parameter to allow also for
  // phase-modulation
}

//-------------------------------------------------------------------------------------------------

SineOscAudioModulePoly::SineOscAudioModulePoly(CriticalSection* lockToUse,
  MetaParameterManager* metaManagerToUse, ModulationManager* modManagerToUse)
  : AudioModulePoly(lockToUse, metaManagerToUse, modManagerToUse) 
{
  ScopedLock scopedLock(*lock);
  setModuleTypeName("SineOscillatorPoly");  // change to SineOscillator later
  createParameters();
}

void SineOscAudioModulePoly::createParameters()
{
  ScopedLock scopedLock(*lock);

  typedef SineOscAudioModulePoly SOM;
  typedef ModulatableParameterPoly Param;
  Param* p;

  p = new Param("Frequency", 20.0, 20000.0, 1000.0, Parameter::EXPONENTIAL);
  addObservedParameter(p);
  p->setValueChangeCallbackPoly([this](double v, int i) { setFrequency(v, i); });
  // later we want a sort of exponential mapping but one that allows for zero values - maybe
  // something based on sinh that is matched to the actual desired exp-shape at the max-value
  // and at the mid-value (or some other selectable value) - or the user prescribes the mid 
  // value - what about frequencies < 0 - that shoudl reverse the direction of the oscillator
  // in the sine osc, we may realizes this by swapping y1,y2 (i think)

  p = new Param("Amplitude", -1.0, +1.0, 1.0, Parameter::LINEAR);
  addObservedParameter(p);
  p->setValueChangeCallbackPoly([this](double v, int i) { setAmplitude(v, i); });

  /*
  p = new Param("Detune", -24.0, +24.0, 0.0, Parameter::LINEAR);
  addObservedParameter(p);
  p->setValueChangeCallbackPoly([this](double v, int i) { setDetune(v, i); });
  */

  // now the parameters don't seem to do anything - that's not a good behavior! maybe we need both
  // monophonic and polyphonic callbacks for the parameters? but nah - that doesn't seem right. 
  // What we want is:
  // -when no modulator is connected to a slider, each voice uses the same value for frequency
 //   and detune and that value is determined by the slider as is
  //  ...but i think, it shoould already behave this way, if the modulated value is correctly
  //  initialized from the unmodulated value - maybe there's something wrong with that
  // -hmm..setFrequency is not called - why? ..ah - because we have no modulator assigned
  //  ...but that's not all - what about doModulationUpdate? do we need to override this in
  //  AudioModulePoly to also call do the polyphonic callbacks?
}

void SineOscAudioModulePoly::allocateVoiceResources() 
{
  if(voiceManager)
    voices.resize(voiceManager->getMaxNumVoices());
  else
    voices.resize(1); // monophonic in absence of a voice manager
}

void SineOscAudioModulePoly::setFrequency(double newFrequency, int voice)
{
  jassert(voice < voices.size());
  double omega = 2*PI*newFrequency / sampleRate;
  voices[voice].setOmega(omega);
}

void SineOscAudioModulePoly::setAmplitude(double newAmplitude, int voice)
{
  jassert(voice < voices.size());
  voices[voice].setAmplitude(newAmplitude);
}

void SineOscAudioModulePoly::setDetune(double newDetune, int voice)
{
  jassert(voice < voices.size());
  RAPT::rsError("Not yet implemented");
}



//=================================================================================================

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
  // adjusts horizontal vs vertical extent of ellipse - at 45°, its a circle


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
  oscCore.setFrequency(RAPT::rsPitchToFreq(noteNumber)); // preliminary - use tuning table
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

  p = new Param("Amplitude", -2.0, 2.0, 1.0, Parameter::LINEAR);
  addObservedParameter(p);
  p->setValueChangeCallback<TriSawOscModule>(this, &TriSawOscModule::setAmplitude);

  p = new Param("Asymmetry", -1.0, 1.0, 0.0, Parameter::LINEAR);
  //p->setMapper(new rsParameterMapperRationalBipolar(-1, +1, 0.80)); // test
  addObservedParameter(p);
  p->setValueChangeCallback<TSO>(tso, &TSO::setAsymmetry);

  // for later, to replace AttackBending and DecayBending parameters - but that doesn't work yet:
  p = new Param("Bending", -1.0, 1.0, 0.0, Parameter::LINEAR);
  //p->setMapper(new rsParameterMapperRationalBipolar(-1, +1, 0.9));
  addObservedParameter(p);
  p->setValueChangeCallback<TriSawOscModule>(this, &TriSawOscModule::setBend);

  p = new Param("BendAsym", -1.0, 1.0, 0.0, Parameter::LINEAR);
  //p->setMapper(new rsParameterMapperRationalBipolar(-1, +1, 0.9));
  addObservedParameter(p);
  p->setValueChangeCallback<TriSawOscModule>(this, &TriSawOscModule::setBendAsym);

  /*
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
  */

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
  *left = *right = amplitude * oscCore.getSample();
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
  freq = RAPT::rsPitchToFreq(noteNumber);
  oscCore.setPhaseIncrement(freq/sampleRate); // maybe times 2 because it ranges from -1..+1?
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


  //// that's not yet very nice:
  //oscCore.setAttackBending(rosic::clip(bend + bendAsym, -1., +1.) );
  //oscCore.setDecayBending( rosic::clip(bend - bendAsym, -1., +1.) );
  //// but maybe it doesn't get any better


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


  // from here: https://github.com/RobinSchmidt/RS-MET/issues/235#issuecomment-424092660
  double bendAbs = bend;
  double target  = 1;
  if(bend < 0) {
    bendAbs = -bend;
    target  = -target;
  }
  oscCore.setAttackBending(juce::jmap(bendAbs,  bendAsym, target));
  oscCore.setDecayBending( juce::jmap(bendAbs, -bendAsym, target));
  // jmap Remaps a normalised value (between 0 and 1) to a target range.
  // This effectively returns (targetRangeMin + value0To1 * (targetRangeMax - targetRangeMin)).
  // why does this work? is really every possible combination of attack-bend and decay-bend 
  // reachable? ...plot attack-bending as funtion of bend and bendAsym

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

//=================================================================================================

FlatZapperModule::FlatZapperModule(CriticalSection *lockToUse, 
  MetaParameterManager* metaManagerToUse, ModulationManager* modManagerToUse)
  : AudioModuleWithMidiIn(lockToUse, metaManagerToUse, modManagerToUse)
{
  ScopedLock scopedLock(*lock);
  setModuleTypeName("FlatZapper");
  createParameters();
}

void FlatZapperModule::createParameters()
{
  ScopedLock scopedLock(*lock);

  typedef rosic::rsFlatZapper FZ;
  FZ* fz = &zapperCore;

  typedef Parameter Param;
  Param* p;

  // Input/output settings
  p = new Param("Level", -48, +96, 0.0, Parameter::LINEAR, 0.01);
  addObservedParameter(p);
  p->setValueChangeCallback<FlatZapperModule>(this, &FlatZapperModule::setLevel);
  // Maybe rename to outputLevel or OutVolume, OutGain

  p = new Param("Exciter", -100, +100, 100.0, Parameter::LINEAR, 0.01); // %
  addObservedParameter(p);
  p->setValueChangeCallback<FlatZapperModule>(this, &FlatZapperModule::setExciterPercent);

  p = new Param("Input", -100, +100, 0.0, Parameter::LINEAR, 0.01); // %
  addObservedParameter(p);
  p->setValueChangeCallback<FlatZapperModule>(this, &FlatZapperModule::setInputPercent);


  // Filter mode settings:
  p = new Param("NumStages", 0.0, 256, 128.0, Parameter::LINEAR, 1.0);
  addObservedParameter(p);
  p->setValueChangeCallback<FZ>(fz, &FZ::setNumStages);
  // Maybe use exp-scaling...or maybe not


  p = new Param("AllpassMode", 0.0, 1.0, 1.0, Parameter::STRING);
  p->addStringValue("OnePole");
  p->addStringValue("Biquad");
  addObservedParameter(p);
  p->setValueChangeCallback<FlatZapperModule>(this, &FlatZapperModule::setAllpassMode);
  // The order in which we add the string-values for the modes must match the order in the switch
  // statement in setAllpassMode()


  // Freq settings:
  p = new Param("FreqLow", 1.0, 20000.0, 15.0, Parameter::EXPONENTIAL);
  addObservedParameter(p);
  p->setValueChangeCallback<FZ>(fz, &FZ::setLowFreq);

  p = new Param("FreqHigh", 10.0, 20000.0, 8000.0, Parameter::EXPONENTIAL);
  addObservedParameter(p);
  p->setValueChangeCallback<FZ>(fz, &FZ::setHighFreq);

  p = new Param("FreqShape", -5, +5, 0.0, Parameter::LINEAR, 0.01);
  addObservedParameter(p);
  p->setValueChangeCallback<FZ>(fz, &FZ::setFreqShape);


  // Q settings:
  p = new Param("QLow", 0.125, 8.0, 1.0, Parameter::EXPONENTIAL, 0.0);
  addObservedParameter(p);
  p->setValueChangeCallback<FZ>(fz, &FZ::setLowQ);

  p = new Param("QHigh", 0.125, 8.0, 1.0, Parameter::EXPONENTIAL, 0.0);
  addObservedParameter(p);
  p->setValueChangeCallback<FZ>(fz, &FZ::setHighQ);

  p = new Param("QShape", -5, +5, 0.0, Parameter::LINEAR, 0.01);
  addObservedParameter(p);
  p->setValueChangeCallback<FZ>(fz, &FZ::setQShape);



  int dummy = 0;
}

void FlatZapperModule::processStereoFrame(double *left, double *right)
{
  // Produce the exitation signal. We use the input and add unit impulses at note events:
  double inL = inputAmp * (*left);
  double inR = inputAmp * (*right);
  if(receivedKey != -1)
  {
    if(receivedVel != 0)
    {
      inL += exciterAmp;
      inR += exciterAmp;
    }
    else
    {
      // todo: maybe trigger a release excitation
    }
  }
  

  //*left = *right = inL;   // for test
  
  *left = *right = masterAmp * zapperCore.getSample(inL); 
  // Preliminary - todo: implement stereo processing by moving the DSP core class to RAPT, 
  // templatizing it and using an instance for <double, rsFloat64x2> here. See implementation of
  // EngineersFilter for how that works.


  // The last noteOn event has been consumed - so we reset the members that store it in noteOn:
  receivedKey = -1;
  receivedVel = -1;
}

void FlatZapperModule::setSampleRate(double newSampleRate)
{
  zapperCore.setSampleRate(newSampleRate);
}

void FlatZapperModule::reset()
{
  zapperCore.reset();
}

void FlatZapperModule::noteOn(int noteNumber, int velocity)
{
  receivedKey = noteNumber;
  receivedVel = velocity;
}

void FlatZapperModule::setAllpassMode(double newMode)
{
  using MD = rosic::rsFlatZapper::Mode;
  int mode = round(newMode);
  switch(mode)
  {
  case 0: zapperCore.setMode(MD::onePole); break;
  case 1: zapperCore.setMode(MD::biquad);  break;
  }
  // The order in which the cases occur here must match the order in which we add the corresponding
  // string values to the "AllpassMode" parameter in createParameters()
}

// Observations:
// -This algorithm is good for zap, kick and tom sounds
// -When retriggering in fast succession like in a machine gun effect, it behaves nicely natural. 
//  Hmm...or well - the signal seem to add up which makes sense. maybe optionally reset the filter
//  states on note on - maybe using a soft-reset parameter ranging from 0 to 1. Hmm - may it's not 
//  so nice after all - try hitting 2 or 3 keys simultaneously - sounds bad. But a soft-reset will
//  behave differently for different filter implementations (DF1 vs DF2 vs SVF vs ...). That may be
//  a good thing, if we let the user also choose the filter structure explicitly but a bad thing
//  if we just hardcode it and want to change it later.
// -Using it as effect on a sawtooth output of a synth is also interesting
//  -the waveform looks different for every pitch
//  -we get a transient before the waveform settles to steady state
//  -putting a saturator after it is also nice - the saturation sound different for every note due 
//   to the different waveshape. A stronger fuzz is also cool
//  -sounds all goa-esque
//  -When the saw-freq is low enough with respct to the length of the zap, one zap per cycle 
//   happens - but when the zap bleeds into the next cycle, things get crazy and the waveform 
//   becomes pseudo-random
// -It can also serve a simple impulse generator when setting numStages to 0.
//
//
// ToDo:
// -Maybe have a 1-pole lowpass, a DC-blocker (Butterworth, adjustable order) and a tilt-filter 
//  built in. Especially the 1-pole lowpass is important because without it, all presets will rely
//  on a post-processing Equalizer which is not good. But the internal filters should all be 
//  optional
// -Add parameter for the allpass mode (biquad, 1-pole)
// -Make a custom GUI that shows a preview of the output and later maybe also the phase-response,
//  phase-delay and group delay response

//=================================================================================================

SweepKickerModule::SweepKickerModule(CriticalSection *lockToUse, 
  MetaParameterManager* metaManagerToUse, ModulationManager* modManagerToUse)
  : AudioModuleWithMidiIn(lockToUse, metaManagerToUse, modManagerToUse)
{
  ScopedLock scopedLock(*lock);
  setModuleTypeName("SweepKicker");
  createParameters();
}

void SweepKickerModule::createParameters()
{
  ScopedLock scopedLock(*lock);

  using SK     = rosic::rsSweepKicker;
  using SKM    = jura::SweepKickerModule;
  using FixPar = jura::Parameter;
  //using AutPar = jura::AutomatableParameter;
  using ModPar = jura::ModulatableParameter;
  // In the framework, we have 3 types of parameters with an incremental feature set:
  // fixed, automatable, modulatable. 

  FixPar* fp;
  //AutPar* ap;
  ModPar* mp;

  // Amplitude parameters:
  mp = new ModPar("Amplitude", -1.0, +1.0, 1.0, Parameter::LINEAR);
  addObservedParameter(mp);
  mp->setValueChangeCallback<SKM>(this, &SKM::setAmplitude);
  // Modulatable to get away without having a built-in amp-env.

  fp = new FixPar("FadeOut", 0.0, 500.0, 100.0, Parameter::LINEAR);
  addObservedParameter(fp);
  fp->setValueChangeCallback<SK>(&core, &SK::setFadeOutTimeMs);

  fp = new FixPar("PassThrough", -1.0, +1.0, 0.0, Parameter::LINEAR);
  addObservedParameter(fp);
  fp->setValueChangeCallback<SKM>(this, &SKM::setPassThroughAmplitude);


  // Frequency parameters:
  // Interpretation of the frequency parameters: we use them unchanged when the incoming note is on
  // the reference key (which is 69). But they respond both to the midi-key with a 100% (i.e. 1 
  // oct/oct) keytrack characteristic:

  fp = new FixPar("FreqHigh", 500.0, 20000.0, 10000.0, Parameter::EXPONENTIAL);
  addObservedParameter(fp);
  fp->setValueChangeCallback<SK>(&core, &SK::setHighFreq);

  fp = new FixPar("FreqLow", 0.0, 400.0, 0.0, Parameter::LINEAR);
  addObservedParameter(fp);
  fp->setValueChangeCallback<SK>(&core, &SK::setLowFreq);


  // Freq envelope parameters:

  fp = new FixPar("SweepTime", 50.0, 500.0, 200.0, Parameter::EXPONENTIAL); // In seconds
  addObservedParameter(fp);
  fp->setValueChangeCallback<SK>(&core, &SK::setSweepTimeInMs);

  fp = new FixPar("Chirp", -1.0, +1.0, 0.0, Parameter::LINEAR);
  addObservedParameter(fp);
  fp->setValueChangeCallback<SK>(&core, &SK::setChirpAmount);

  fp = new FixPar("ChirpShape", -1.0, +1.0, 0.0, Parameter::LINEAR);
  addObservedParameter(fp);
  fp->setValueChangeCallback<SK>(&core, &SK::setChirpShape);


  // Oscillator parameters:

  mp = new ModPar("Phase", -180.0, +180.0, 0.0, Parameter::LINEAR);
  addObservedParameter(mp);
  mp->setValueChangeCallback<SK>(&core, &SK::setStartPhase);

  mp = new ModPar("PhaseStereoShift", -180.0, +180.0, 0.0, Parameter::LINEAR);
  addObservedParameter(mp);
  mp->setValueChangeCallback<SK>(&core, &SK::setStereoPhaseShift);
  // It is interesting to compare the stereo phase shift dialed in here to the same value with the
  // PhaseStereoizer as effect. It does indeed sound different


  // ToDo:
  // -Add FreqHighByVel, FreqScale (modulatable), ChirpByVel, Formula (maybe - let the use choose 
  //  different formulas for freq-env), SinToSaw/PhaseShape/WarpHalf (see Straighliner's oscs),
  //  Detune (maybe - but would be somewhat redundant with FreqScale - but may respond differently
  //  to modulation - linear vs exponential FM - so it might be worthwhile - but we could also 
  //  modulate a single FreqScale parameter with different modes, I guess - and also: for linear FM
  //  it would have to be an additive offset parameter, I think)
  // -Add descriptions for the parameters
}

void SweepKickerModule::processStereoFrame(double* left, double* right)
{
  double tmpL, tmpR;
  core.getSampleFrameStereo(&tmpL, &tmpR);
  *left  = passThroughAmp * *left   +  amplitude * tmpL;
  *right = passThroughAmp * *right  +  amplitude * tmpR;
}