namespace rosic {
namespace Sampler {

void Processor::addParameter(Opcode opcode)
{
  params.push_back(Parameter(opcode));
}

void Processor::removeParameter(Opcode opcode)
{
  int i = findParameter(opcode);
  if(i != -1)
    RAPT::rsRemove(params, (size_t)i);
}

int Processor::replaceOpcode(Opcode oldOpcode, Opcode newOpcode)
{
  int foundAt = findParameter(oldOpcode);
  RAPT::rsAssert(foundAt != -1, "No such opcode found");
  if(foundAt != -1)
    params[foundAt].setOpcode(newOpcode);
  return foundAt;
}

void Processor::setParameter(Opcode opcode, float value)
{
  size_t i = 0;
  for(i = 0; i < params.size(); i++) {
    if(params[i].getOpcode() == opcode) {
      params[i].setValue(value);
      return; }}
  RAPT::rsError("Parameter not found in Processor::setParameter");
}

Parameter* Processor::getParameter(Opcode op)
{
  for(size_t i = 0; i < params.size(); ++i) {
    if(params[i].getOpcode() == op)
      return &params[i]; }
  return nullptr;
}

void Processor::setParametersToDefaults(int index)
{
  SfzCodeBook* cb = SfzCodeBook::getInstance();
  for(size_t i = 0; i < params.size(); i++) {
    Opcode op = params[i].getOpcode();
    float defVal = cb->opcodeDefaultValue(op, index);
    params[i].setValue(defVal); }
}

void Processor::prepareToPlay(uchar key, uchar vel, double sampleRate)
{
  this->key = key;
  this->vel = vel;
  dirty = true;

  // Maybe we need to init modulatedValues to the unmodulated ones:
  for(size_t i = 0; i < params.size(); i++)  // maybe factor out into function and call it
    params[i].initModulatedValue();          // also from  setParametersToDefaults

  updateCoeffs(sampleRate);
  resetState();
}

int Processor::findParameter(Opcode opcode)
{
  for(size_t i = 0; i < params.size(); ++i) 
    if(params[i].getOpcode() == opcode) 
      return (int) i;
  return -1;
}

//=================================================================================================

LowFreqOsc::LowFreqOsc()
{
  type = OpcodeType::FreeLfo;
  params.reserve(3);                      // index
  addParameter(Opcode::lfoN_freq);        //   0
  addParameter(Opcode::lfoN_amp);         //   1
  //addParameter(Opcode::lfoN_fade);
  //addParameter(Opcode::lfoN_delay);
  // ToDo: phase, wave, sync, ...
  // allow the "wave" parameter to be something like "Samples/SingleCyle/Sine.wav". When the 
  // implementation works, it can used as model for the SamplePlayer object which should replace
  // the direct handling of samples in the RegionPlayer. Maybe the "wave" parameter could have a 
  // setting "sample" in which case we also store a pointer to an AudioStream object containing
  // the waveform
}

void LowFreqOsc::updateCoeffs(double fs)
{
  float k = 0.01f;     // scales percents to raw factors
  core.setup(params[0].mv(), k*params[1].mv(), 0.f, 0.f, 0.f, (float)fs);
  //core.setup(params[0].mv(), 1.f, 0.f, 0.f, 0.f, (float)fs);

  dirty = false;
}

LowFreqOscAmp::LowFreqOscAmp()
{
  type = OpcodeType::AmpLfo;
  replaceOpcode(Opcode::lfoN_freq, Opcode::amplfo_freq);
  //replaceOpcode(Opcode::lfoN_amp,  Opcode::amplfo_amp);


  //removeParameter(Opcode::lfoN_amp);
  // It would then be better to not create it in the first place, memory-wise...not a big deal,
  // but still. But no - we shouldn't remove it - if we do, we also need to re-implement 
  // updateCoeffs because it accesses the params array and assumed a certain content
}

LowFreqOscFil::LowFreqOscFil()
{
  type = OpcodeType::FilterLfo;
  replaceOpcode(Opcode::lfoN_freq, Opcode::fillfo_freq);
  //replaceOpcode(Opcode::lfoN_amp,  Opcode::fillfo_amp);
}




//=================================================================================================

EnvGen::EnvGen()
{
  type = OpcodeType::FreeEnv;
  using OC = Opcode;
  params.reserve(12);                      // index
  addParameter(OC::adsrN_start);           //   0
  addParameter(OC::adsrN_delay);           //   1
  addParameter(OC::adsrN_attack);          //   2
  addParameter(OC::adsrN_peak);            //   3
  addParameter(OC::adsrN_hold);            //   4
  addParameter(OC::adsrN_decay);           //   5
  addParameter(OC::adsrN_sustain);         //   6
  addParameter(OC::adsrN_release);         //   7
  addParameter(OC::adsrN_end);             //   8
  addParameter(OC::adsrN_attack_shape);    //   9
  addParameter(OC::adsrN_decay_shape);     //  10
  addParameter(OC::adsrN_release_shape);   //  11

  // ToDo: 
  // -Add the velocity-scaling params: vel2delay, vel2attack, vel2hold, etc. defined in sfz for
  //  ampeg, fileg, etc.
  // -Add key2attack, key2decay, etc. and/or maybe have a general key2length or key2timescale as we
  //  have in Straightliner. These are not defined in sfz but make a lot of sense musically.
  // -peak is actually not defined in sfz - so maybe we should remove it but maybe it's useful 
  //  enough to keep it
}

void EnvGen::updateCoeffs(double sampleRate)
{
  const std::vector<Parameter>& p = params;
  float fs = (float) sampleRate;  // scales seconds to samples
  float k  = 0.01f;               // scales percents to raw factors
  core.setup(p[0].mv()*k, p[1].mv()*fs, p[2].mv()*fs, p[3].mv()*k, p[4].mv()*fs, p[5].mv()*fs,
    p[6].mv()*k, p[7].mv()*fs, p[8].mv()*k, p[9].mv(), p[10].mv(), p[11].mv());
  dirty = false;
}

EnvGenAmp::EnvGenAmp()
{
  // We just replace the opcodes to which our parameters listen. Functionally, the ampeg opcodes 
  // work exactly the same as the corresponding general adsrN, so we can re-use the code. Just the
  // opcodes have different names reflecting the specialization for the specific modulation target.

  type = OpcodeType::AmpEnv;
  using OC = Opcode;
  replaceOpcode(OC::adsrN_start,         OC::ampeg_start);
  replaceOpcode(OC::adsrN_delay,         OC::ampeg_delay);
  replaceOpcode(OC::adsrN_attack,        OC::ampeg_attack);
  replaceOpcode(OC::adsrN_peak,          OC::ampeg_peak);
  replaceOpcode(OC::adsrN_hold,          OC::ampeg_hold);
  replaceOpcode(OC::adsrN_decay,         OC::ampeg_decay);
  replaceOpcode(OC::adsrN_sustain,       OC::ampeg_sustain);
  replaceOpcode(OC::adsrN_release,       OC::ampeg_release);
  replaceOpcode(OC::adsrN_end,           OC::ampeg_end);
  replaceOpcode(OC::adsrN_attack_shape,  OC::ampeg_attack_shape);
  replaceOpcode(OC::adsrN_decay_shape,   OC::ampeg_decay_shape);
  replaceOpcode(OC::adsrN_release_shape, OC::ampeg_release_shape);


  // ToDo:
  // -Do the same for pitcheg 
  // -We also need to do soemthing similar for the LFO
  // -For the pre-allocation of env objects, it would actually be much more desirable, if the
  //  different types of specialized envelopes would be interchangable at runtime, i.e. EnvGenAmp
  //  is not a subclass of EnveGen but instead, we can somehow at runtime switch a general EG
  //  into an AmpEG (and back). Then we can just pre-allocate an array of general EGs and use them
  //  for all purposes. ...not yet sure, how to implement this, though - might be tricky.
}

EnvGenFil::EnvGenFil()
{
  type = OpcodeType::FilterEnv;
  using OC = Opcode;
  replaceOpcode(OC::adsrN_start,         OC::fileg_start);
  replaceOpcode(OC::adsrN_delay,         OC::fileg_delay);
  replaceOpcode(OC::adsrN_attack,        OC::fileg_attack);
  replaceOpcode(OC::adsrN_peak,          OC::fileg_peak);
  replaceOpcode(OC::adsrN_hold,          OC::fileg_hold);
  replaceOpcode(OC::adsrN_decay,         OC::fileg_decay);
  replaceOpcode(OC::adsrN_sustain,       OC::fileg_sustain);
  replaceOpcode(OC::adsrN_release,       OC::fileg_release);
  replaceOpcode(OC::adsrN_end,           OC::fileg_end);
  replaceOpcode(OC::adsrN_attack_shape,  OC::fileg_attack_shape);
  replaceOpcode(OC::adsrN_decay_shape,   OC::fileg_decay_shape);
  replaceOpcode(OC::adsrN_release_shape, OC::fileg_release_shape);
}

//=================================================================================================

MidiController::MidiController()
{
  using OC = Opcode;
  type = OpcodeType::MidiCtrl;
  params.reserve(1);                 // index
  addParameter(OC::controlN_index);  //   0

  // Maybe add the following parameters:
  // -controlN_amount or _depth or _scale: scales the normalized output by a factor, maybe unit
  //  should be percent and default 100
  // -controlN_shift or _offset: shifts the controller output by a constant, maybe it should also
  //  be expressed in %, default is zeor, of course
  // -controlN_neutral: defines the neutral value as integer midi number in 0..127. Is actually 
  //  also am offset, but applied before normalizing to 0..1 or 0..100%
  // -controlN_smooth: smoothing time in seconds. We can have a linear smoother in the DSP core
  //  object. Default is zero.
  // -controlN_quantize: an optional quantization, defaults to zero
}

void MidiController::processFrame(float* L, float* R)
{ 
  // Outputs the normalized value of the controller in the range 0..1 where a midi value of 0 maps 
  // to 0.f and a midi value of 127 maps to 1.f.

  RAPT::rsAssert(playStatus != nullptr);
  if(playStatus) {
    RAPT::rsUint8 rawVal = playStatus->getMidiControllerCurrentValue(ctrlIndex);
    *L = *R = (1.f/127.f) * (float)rawVal; }
  else {
    *L = *R = 0.f; }
    // We are playing safe here with the if-conditional in the spirit of defensive programming. 
    // Later when the code stabilizes, maybe we can just assume that playStatus never is a nullptr
    // and optimize this branch away.
}

void MidiController::updateCoeffs(double sampleRate)
{
  ctrlIndex = (int) params[0].mv();
}

//=================================================================================================

Amplifier::Amplifier()
{
  type = OpcodeType::Amplifier;
  params.reserve(8);                      // index
  addParameter(Opcode::amplitudeN);       //   0
  addParameter(Opcode::volumeN);          //   1
  addParameter(Opcode::panN);             //   2
  addParameter(Opcode::widthN);           //   3
  addParameter(Opcode::positionN);        //   4
  addParameter(Opcode::ampN_veltrack);    //   5
  addParameter(Opcode::ampN_keytrack);    //   6
  addParameter(Opcode::ampN_keycenter);   //   7

  // Having to pass a magic number to reserve() is bad and error-prone -> try to find a better
  // way. The number of parameters is actually known at compile time. Maybe use std::array 
  // instead of std::vector...hmm...but the number varies between the subclasses, so the array
  // could not be a baseclass member then...hmm...Maybe SignalProcessor could have an int as
  // template parameter - but no: then i can't treat them uniformly in the DSP chain...ok, then 
  // maybe try to set up a unit test that fires when we reserve the wrong size here. Maybe we can 
  // use an initializer list. But it's actually quite useful to have a more verbose textual "list" 
  // here to refer to to figure out the indices of the parameters. Storing the parameters on the 
  // heap may actually be beneficial for the memory layout of the DSPs because then they don't 
  // share the same memory with the algo-parameters. The params are only accessed occasionally 
  // during processing (in case of handling midi controllers) but most of the time, it's good that 
  // they don't intefere with the algo params (stored directly in the DSP objects).

  // Maybe add these:
  // https://sfzformat.com/opcodes/amplitude
}

void Amplifier::processFrame(float* L, float* R) 
{
  core.processFrame(L, R);
  // ToDo:
  // -Let the baseclass maintain a pointer to some sort of MidiStatus object which contains a 
  //  "dirty" flag which is set whenever the status changes, for example because a 
  //  control-change was received.
  // -Inspect the flag here, if it is dirty, we have two options:
  //  -Recalculate out coeffs immediately, taking into account the new settings of controllers
  //   etc, or:
  // -Set up this object for a parameter transition (smoothing) by initializing a sampleCounter
  //  to numSmoothingSamples and checking that counter in each sample and if it nonzero, do an 
  //  update step and count down - when zero is reached, we have reached the new target 
  //  parameters
}

void Amplifier::processBlock(float* L, float* R, int N) 
{
  for(int n = 0; n < N; n++)
    processFrame(&L[n], &R[n]);
}

void Amplifier::updateCoeffs(double sampleRate)
{

  // Extract nominal values from the parameters:

  float amplitude = params[0].mv();
  float volume    = params[1].mv();
  float pan       = params[2].mv();
  float width     = params[3].mv();
  float position  = params[4].mv();

  // Extract key/vel modifiers:
  float ampN_veltrack  = params[5].mv();
  float ampN_keytrack  = params[6].mv();
  float ampN_keycenter = params[7].mv();

  // Apply modifiers:
  volume += ampN_veltrack * 0.01f * 40 * log10f(127.f/(float)vel);
  volume += ampN_keytrack * ((float)key - ampN_keycenter);
  // Unit of veltrack is percent and the formula is dB = 20 log (127^2 / Velocity^2). Keytracking 
  // is adjusted in dB per key. See https://sfzformat.com/legacy/

  // Set up the core:
  core.setup(amplitude*0.01f, volume, pan, width, position);
  dirty = false;

  // ToDo:
  // -get rid of the getValue() calls by allowing the params to convert to float
  // -it's not ideal that this code depends on the order, how we add the params in the 
  //  constructor - try to avoid that - not sure, if that's possible
}

//=================================================================================================

Filter::Filter() 
{ 
  using OC = Opcode;
  type = OpcodeType::Filter;
  params.reserve(6);                 // index
  addParameter(OC::filN_type);       //   0
  addParameter(OC::cutoffN);         //   1
  addParameter(OC::resonanceN);      //   2
  addParameter(OC::filN_keytrack);   //   3
  addParameter(OC::filN_keycenter);  //   4
  addParameter(OC::filN_veltrack);   //   5
  // ToDo
  // -Maybe support peak filters, re-use the resonance parameter for their gain, introduce a 
  //  bandwidth parameter ...on the other hand, we have the eqN_... opcodes, so maybe supporting 
  //  peak would be redundant - we'll see
  // -Maybe support a reso_drive opcode to allow self oscillation in ladder filters
  // -More parameters: filN_bw: bandwidth in octaves (for peak and shelv filters), filN_design =
  //  default, moog, biquad, butter, cheby1, cheby2, ellip, ...
  //
}

void Filter::processFrame(float* L, float* R) 
{ 
  core.processFrame(L, R); 
}

void Filter::processBlock(float* L, float* R, int N) 
{
  for(int n = 0; n < N; n++)
    processFrame(&L[n], &R[n]);
}

void Filter::updateCoeffs(double fs)
{
  // This is still somewhat ugly:
  FilterType sfzType = (FilterType)(int)params[0].getValue();
  FilterCore::Type coreType = convertTypeEnum(sfzType);

  // Extract numeric parameters:
  float cutoff    = params[1].mv();
  float resonance = params[2].mv();
  float keytrack  = params[3].mv();
  float keycenter = params[4].mv();
  float veltrack  = params[5].mv();

  // Apply modifiers to cutoff:
  float pitchOffset = ((float)key - keycenter) * keytrack * 0.01f;  // pitch-offset from keytrack
  float scl = 1.f - ((vel-1)/126.f);                                // 0 at vel = 127, 1 at vel = 1
  pitchOffset += veltrack * 0.01f * scl;                            // pitch-offset from veltrack
  cutoff *= RAPT::rsPitchOffsetToFreqFactor(pitchOffset);

  // Set up core:
  core.setupCutRes(coreType, cutoff*float(2*PI/fs), resonance);
  dirty = false;

  // ToDo:
  // Verify the formula used for velocity tracking. It's just a guess based on what I think, the 
  // behavior should be. I think, at vel = 127, the cutoff should be unmodified and at vel=1, the 
  // cutoff should be modified by the given amount of veltrack in cents. This is achieved by the
  // formula with a percpetually sensible transition in between. But I don't know, if that's what
  // a reference implementation does. 
  // Maybe the implementation of key/vel tracking should be based on using the modulation system by
  // allowing midi key/vel as modulation sources? But that may make key/vel tracking more expensive 
  // because we would then need modulation connections which we otherwise don't need.
}

FilterCore::Type Filter::convertTypeEnum(FilterType sfzType)
{
  // Conversion of filter type enum values used in the sfz data and those used in the dsp core.
  // Maybe we should try to avoid the translation step between the core-enum and sfz-enum by 
  // using a single enum for both. I'm not yet sure, if that's practical - we'll see. It could
  // turn out to be problematic when we want to use bit-twiddling of the enum-values to switch 
  // between different filter topologies in the core while the sfz-type number must allow for
  // lossless roundtrip with a float.
  using TC = FilterCore::Type;      // enum used in the DSP core
  using TO = FilterType;            // enum used in the sfz opcode
  switch(sfzType)
  {
  case TO::lp_6:   return TC::FO_Lowpass;
  case TO::hp_6:   return TC::FO_Highpass;

  case TO::lp_12:  return TC::BQ_Lowpass;        // ToDo: Use SVF as default implementation
  case TO::hp_12:  return TC::BQ_Highpass;       // for 2nd order filters...maybe...
  case TO::bp_6_6: return TC::BQ_Bandpass_Skirt;
  case TO::br_6_6: return TC::BQ_Bandstop;


    //case TO::lp_12: return TC::SVF_Lowpass_12;
    //case TO::hp_12: return TC::SVF_Highpass_12;
  }
  //RAPT::rsError("Unknown filter type in convertTypeEnum.");
  //return TC::Unknown;
  // This may actually happen when the user only defines cutoff but not fil_type. In this 
  // case, sfz prescribes that the default mode is lpf_2p.

  return TC::BQ_Lowpass;
}
// todo: avoid this conversion - use the same enum in both, the sfz codebook and the 
// FilterCore just like we do with the waveshaper's distortion shapes

//=================================================================================================

Equalizer::Equalizer()
{
  type = OpcodeType::Equalizer;
  params.reserve(3);
  addParameter(Opcode::eqN_gain);
  addParameter(Opcode::eqN_freq);
  addParameter(Opcode::eqN_bw);

  // Maybe add these:
  // https://sfzformat.com/opcodes/eq_type
  // https://sfzformat.com/opcodes/eqN_type
}

void Equalizer::processFrame(float* L, float* R)
{ 
  core.processFrame(L, R); 
}

void Equalizer::processBlock(float* L, float* R, int N)
{
  for(int n = 0; n < N; n++)
    processFrame(&L[n], &R[n]);
}

void Equalizer::updateCoeffs(double fs)
{
  core.setupGainFreqBw(
    FilterCore::Type::BQ_Bell,
    params[0].mv(),
    params[1].mv() * float(2*PI/fs),
    params[2].mv()
  );
  dirty = false;
}

//=================================================================================================

WaveShaper::WaveShaper()
{
  type = OpcodeType::WaveShaper;
  params.reserve(3);
  addParameter(Opcode::distortN_shape);
  addParameter(Opcode::distortN_drive);  // maybe drive should be in dB
  addParameter(Opcode::distortN_dc);
}

void WaveShaper::processFrame(float* L, float* R)
{ 
  core.processFrame(L, R); 
}

void WaveShaper::processBlock(float* L, float* R, int N)
{
  for(int n = 0; n < N; n++)
    processFrame(&L[n], &R[n]);
}

void WaveShaper::updateCoeffs(double sampleRate)
{
  //core.setup((DistortShape)(int)params[0].getValue(), params[1].getValue(),
  //  params[2].getValue(), 1.f, 0.f, 0.f);

  core.setup((DistortShape)(int)params[0].getValue(), params[1].mv(),
    params[2].mv(), 1.f, 0.f, 0.f);

  // we use getValue for the unmodulatable values and mv() for the modulatable ones

  dirty = false;
  // ToDo: we actually need to retrieve the modulatedValue not the (nominal) value
}

//=================================================================================================

DspResourcePool::DspResourcePool(PlayStatus* playStatusToUse)
  : playStatus(playStatusToUse)
{
  allocateEffects();
  allocateModulators();
  allocateConnectors();
  // Maybe don't do this on construction. Maybe client code should explicitly request this when a
  // new sfz file ist loaded. The content of the file should determine, how many of each type we 
  // need to allocate. Maybe have an allocate function that takes a reference to the SfzInstrument
  // object. Or: the xml file should say how much to pre-allocate - maybe it can have an option 
  // that let's use allocate heuristically based on the sfz content. Or maybe define sfz opcodes
  // for pre-allocation: max_filters, max_amplifiers, max_connections, etc.
}

void DspResourcePool::allocateEffects()
{
  amplifiers.init(64);
  filters.init(64);
  equalizers.init(64); 
  waveShapers.init(8);
  // These numbers are ad hoc and preliminary. We need to do something more sensible here later. 
  // Perhaps, this function should be called when a new sfz is loaded and it should have arguments
  // for how many objects of each type are needed. The engine should analyze, how many filters, 
  // waveshapers, etc. could be needed in the worst case, fill a suitable data structure with that 
  // information and pass it in.

  // Maybe instead of manually resizing everything, do it programmatically, maybe using a
  // SignalProcessorFactory to create the objects. But then: how do we keep track of the different
  // kinds of modules?
}

Processor* DspResourcePool::grabEffect(OpcodeType type)
{
  using OT = OpcodeType;
  Processor* p = nullptr;
  switch(type)
  {
  case OT::Amplifier:  p = amplifiers.grabItem();  break;
  case OT::Filter:     p = filters.grabItem();     break;
  case OT::Equalizer:  p = equalizers.grabItem();  break;
  case OT::WaveShaper: p = waveShapers.grabItem(); break;
  };
  return p;
}

void DspResourcePool::repositEffect(Processor* p)
{
  using OT = OpcodeType;
  int i = -1;
  switch(p->getType())
  {
  case OT::Amplifier:  i = amplifiers.repositItem(p);  break;
  case OT::Filter:     i = filters.repositItem(p);     break;
  case OT::Equalizer:  i = equalizers.repositItem(p);  break;
  case OT::WaveShaper: i = waveShapers.repositItem(p); break;
  }
  RAPT::rsAssert(i != -1, "Reposited processor was not in pool");
}

void DspResourcePool::allocateModulators()
{
  int N = 64; // Preliminary. Maybe it should depend on the patch. But even if not, N=64 is too
              // small in practice.

  freeEnvGens.init(N);
  ampEnvGens.init(N);
  filEnvGens.init(N);

  freeLowFreqOscs.init(N);
  ampLowFreqOscs.init(N);
  filLowFreqOscs.init(N);

  midiControllers.init(N);
  for(int i = 0; i < N; i++)
    midiControllers.getItemPointer(i)->setPlayStatusToUse(playStatus);
}

Processor* DspResourcePool::grabModulator(OpcodeType type)
{
  using OT = OpcodeType;
  Processor* p = nullptr;
  switch(type)
  {
  case OT::FreeEnv:   p = freeEnvGens.grabItem(); break;
  case OT::AmpEnv:    p = ampEnvGens.grabItem();  break;
  case OT::FilterEnv: p = filEnvGens.grabItem();  break;
  //case OT::PitchEnv:  p = envGens.grabItem(); break;

  case OT::FreeLfo:   p = freeLowFreqOscs.grabItem(); break;
  case OT::AmpLfo:    p = ampLowFreqOscs.grabItem();  break;
  case OT::FilterLfo: p = filLowFreqOscs.grabItem();  break;
  //case OT::PitchLfo:  p = lowFreqOscs.grabItem(); break;

  case OT::MidiCtrl: p = midiControllers.grabItem(); break;

  default: { RAPT::rsError("Unknown modulator type"); }
  };
  return p;
  // ToDo: maybe consolidate the cases that return from envGens or lowFreqOscs into some sort of
  // if-statement that checks, if type is within some range, like
  //   if( isLfo(type)    ) return lowFreqOscs.grabItem();
  //   if( isEnvGen(type) ) return envGens.grabItem();
  // and only then enter the switch statement to switch between other kinds of modulators. Maybe
  // also have isMidiInput(type) ...the implementation of midi input controllers can be simple 
  // dsp-wise: it needs to hold a pointer to the PlayStatus object and just return the value 
  // stored there as constant (but changing when the user changes it - but these changes will be
  // seen in the PlayStatus anyway, so it should automatically work)
}

Processor* DspResourcePool::grabProcessor(OpcodeType type)
{
  if(SfzCodeBook::isEffectSetting(type))    return grabEffect(type);
  if(SfzCodeBook::isModSourceSetting(type)) return grabModulator(type);
  RAPT::rsError("type should be and effect or modulator setting");
  return nullptr;
}

void DspResourcePool::repositModulator(Processor* p)
{
  using OT = OpcodeType;
  int i = -1;
  switch(p->getType())
  {
  case OT::FreeEnv:   i = freeEnvGens.repositItem(p); break;
  case OT::AmpEnv:    i = ampEnvGens.repositItem(p);  break;
  case OT::FilterEnv: i = filEnvGens.repositItem(p);  break;
  //case OT::PitchEnv:  i = envGens.repositItem(p); break;

  case OT::FreeLfo:   i = freeLowFreqOscs.repositItem(p); break;
  case OT::AmpLfo:    i = ampLowFreqOscs.repositItem(p);  break;
  case OT::FilterLfo: i = filLowFreqOscs.repositItem(p);  break;
  //case OT::PitchLfo:  i = lowFreqOscs.repositItem(p);  break;

  case OT::MidiCtrl: i = midiControllers.repositItem(p); break;

  }
  RAPT::rsAssert(i != -1, "Reposited processor was not in pool");
}





}} // namespaces

//=================================================================================================
/*

ToDo:

-Make it possible to let the user adjust a parameter and have the effect to be immediately audible
 not only on the next noteOn. 
 -I think, we need to call Processor::updateCoeffs() on any such parameter change on all currently 
  running effects. We need to loop through the activePlayers and then through their effect-chain 
  and for each effect, figure out if it is affected by the parameter change and if so, call 
  updateCoeffs() on it.
 -This call will set dirty = false;  ...I guess that's OK but we need to make sure that this 
  doesn't interfere with modulation system - probably not.

Notes:

Modulators should be stereo. Stereo LFOs should have a parameter for the phase offset between 
modulator for left and right channel. For envelopes, maybe a delay would be more appropriate. The 
LFOs should support the sample opcode to load single cycle LFO waveforms. 

Maybe use rsFloat32x4 to pass around signals. We could use the 4 channels for L/R/M/S and process
all simultaneously. ...but for delayline based effects, that may imply to double the memory 
requirements. If a filter is modulated by a stereo LFO, it would also need stereo fields for its
coeffs. For the state vars, this is needed anyway. Or maybe all effects should have an opcode
stereo_mode with options: linked, left_right, mid_side. How would we realize a mid/side eq, for 
example? Maybe we need opcodes eqN_gain_left, eqN_gain_right, eqN_gain_mid, eqN_gain_side which
accumulate to the global eqN_gain opcode? But how would that work for a cutoff of a lowpass?
Maybe we should have stereo and/or 4 channels versions of the effects and select from those 
options according to which opcodes are defined.


rsSamplerFilter:
 -Controlled by opcodes: 
  -done: cutoffN, resonanceN, filN_type,
  -todo: filN_gain, filN_slope/tilt
 -Modes:
  -done: lpf_1p, lpf_2p, hpf_1p, hpf_2p, bpf_2p (preliminary)
  -todo: eq_peak, eq_loshelf, eq_hishelf, eq_tilt, brf_2p (look up what other engines do, be 
         compatible)
 -Maybe make also a struct for very basic 1-pole/1-zero filter. It can be realized by the biquad 
 (and by the other structures, too), but maybe it's more efficient to do it like that. I expect 
 that a simple 1st order lowpass is a quite common thing to use.
 -Use the filter also for the equalizer opcode. No need to define a different class for that. Maybe
  extend sfz to support 4 instead of 3 bands when we later can realize 2 bands per filter...
-Implement Martin Vicanek's design formulas based on the IIT for the poles and magnitude matching for
 the zeros. Maybe implement a general biquad s-to-z transform based on that idea but try to also 
 derive optimized cookbook formulas for the same types as the RBJ cookbook, just not using bilinear
 trafo but Martin's approach.
-For the filter realization, use a ZDF-SVF design. Maybe try to directly derive cookbook formulas
 for the SVF.

 -to figure out the correct design formulas, check source code of sfizz and linuxsampler
  https://github.com/sfztools/sfizz/blob/develop/src/sfizz/dsp/filters/filters_modulable.dsp
  https://github.com/sfztools/sfizz/blob/develop/src/sfizz/dsp/filters/rbj_filters.dsp

  https://github.com/sfztools/sfizz/blob/develop/src/sfizz/dsp/filters/sfz_filters.dsp
  sfzGetQFromSlope
  cutoff = hslider("[01] Cutoff [unit:Hz] [scale:log]", 440.0, 50.0, 10000.0, 1.0) : max(1.0) : min(20000.0);
  Q = vslider("[02] Resonance [unit:dB]", 0.0, 0.0, 40.0, 0.1) : max(-60.0) : min(60.0) : ba.db2linear;
  ...

  https://github.com/sfztools/sfizz/blob/develop/src/sfizz/gen/filters/sfzPeq.hxx
  https://github.com/sfztools/sfizz/blob/develop/src/sfizz/gen/filters/sfzLpf2p.hxx
  https://github.com/sfztools/sfizz/blob/develop/src/sfizz/gen/filters/sfzHpf2p.hxx

  https://github.com/sfztools/sfizz/blob/develop/src/sfizz/effects/Filter.cpp
  _filter.prepare(_desc.cutoff, _desc.resonance, _desc.gain);
  ..look for that function ...  _filter is of class sfz::Filter

  https://github.com/sfztools/sfizz/blob/acd866fd3d247d2fc659593cac96e88e801c29e2/src/sfizz/SfzFilter.h
  https://github.com/sfztools/sfizz/blob/acd866fd3d247d2fc659593cac96e88e801c29e2/src/sfizz/SfzFilter.cpp


  https://www.electronics-tutorials.ws/filter/second-order-filters.html


  https://github.com/linuxsampler/linuxsampler/blob/master/src/engines/common/BiquadFilter.h


rsSamplerWaveShaper:
-Try to find a similar parametrization for each of the available shapes.
-One parameter should be interpretable as "hardness". it should go from 0 to 100% where 100% maps 
 to a hard-clipper and 0% should perhaps be the default, suing the prototype sigmoid without 
 modification like using e.g. tanh(drive * x) directly without applying a "hardening" 
 transformation.
 -for tanh, maybe try functions like tanh(a*x^b), tanh(a*x+b*x^3), tanh(sinh(a*x)/x), 
  tanh((a*x+b*x^3)/(c*x+d*x^3))
 -try tanh(asinh(a*x)/a) or asinh(tanh(a*x)/a)  
  https://www.desmos.com/calculator/uiobrmwgyy, https://www.desmos.com/calculator/go7o7j6eil

rsSamplerAmplifier:
-Controlled by opcodes:
 -done (preliminary): amplitudeN (new), volumeN (sfz), panN (sfz), widthN (sfz), positionN (sfz)
-internal algo parameters should be gLL, gLR, gRL, gRR - the 4 gains for the channel-mix matrix
-user params should be: scale (linear overall scale factor), pan, width and as 3rd...dunno - maybe
 pos? whatever is left to determine the rest.




-rsSamplerSlopeFilter
 -should be based on 2 shelving filters using the Vicanek designs
 -numerically optimize the values of the free design parameters, i.e. the shelving freqs and Qs. 
  The optimization criterion should be the match between an ideal slope and what the filter 
  actually produces, taking only the freq-range 20-20k into account. We specifically don't care what
  happens below 20 Hz - the ideal 1/f^s slope would approach infinity at DC which we clearly don't 
  want. It's actually good that our filter deviates from the mathematically idealized function.
  ...hmm...but what sort of slope should we choose for the optimization? Maybe 6 dB/oct? Or maybe
  the freqs and Qs should depend on the gain setting? Maybe do the optimization for different 
  slopes (1,2,3,...,48) and fit simple functions to the results.
 -maybe it should not a module in its own right but rather yet another mode for the regular filter 
  and/or eq like lpf_2p, peak_2p, tilt_2p, tilt_4p, see:
  https://hofa-plugins.de/plugins/4u-dynamictilteq/
  https://www.elysia.com/de/plugins/niveau-filter/

-rsSamplerComb and/or rsSamplerDelay
 -maybe it should be a delay bank controlled by opcodes like:
  delayN_timeX or combN_freqX

ToDo:
-Maybe at some later stage, generalize the dspChain in the sampler engine to a more flexible 
 dspGraph where the routing is not necessarily fixed to x -> 1 -> 2 -> 3 -> y, say (where x,y stand
 for input and output and the numbers stand for the dsps). Maybe the routing could be controlled in 
 the sfz file via a new opcode dsp_graph_... where the ... could be a string like x>1_1>2_2>3_3>y 
 to build the simple chain. The string x>1_x>2_1>3_2>3_3>y could mean: input goes into 1,2 in 
 paralell, 1,2 get summed and go into 3, 3 goes to out. Maybe using > could be problematic when 
 embedding such strings in xml files? maybe use x-1_x-2_1-3_2-3_3-y instead? Or some other symbol?
 Maybe we need a constraint that the left number (source node) must always be less than the right 
 number (target node). This will avoid loops and presumably make parsing easier. Maybe the 
 restriction can be lifted later...
-We actually already need a similar syntax to specify the wiring of the modulations. Maybe 
 something like lfo2_to_eq5_freq, eg2_to_cutoff3, etc. ...so, maybe for the DSPs, we should have 
 one opcode per wire instead of specifying the complete wiring in a single opcode? Maybe something
 like in_to_fil1=1, fil1_to_fil2=1, fil2_to_out=1 specifies the chain in -> fil1 -> fil2 -> out and 
 that wiring is implicitly assumed if nothing else is said?
-Test the computational load incurred by the on-demand assembling of the DSP chain. Maybe 
 pre-assemble the chains for all notes, if it turns out to be a performance problem. That would 
 increase the RAM usage beause we would have to keep around all those pre-built DSP chains.
-Maybe for certain effects, we should have an "algo" parameter, for example, for a pitch-shifter,
 it could have values: granular, spectral, etc. If it makes sense to implement different 
 algorithms as different DSP objects (because different algos may differ vastly in resource
 requirements), maybe the sampler-pool should abstract this away in some way - maybe by letting
 grabProcessor have another integer parameter for selecting the algo. For repositing, we may 
 identify the object class via RTTI (i.e. dynamic_cast or typeid or whatever).


Ideas:
-Maybe use std::variant for the DSP chain, see here: https://www.youtube.com/watch?v=h-zy1hBqT74
 towards the end. The same thing should actually also be possible with unions and then maybe with
 a struct consisting of such a union and a function-pointer (or a bunch of them) for the process
 etc. functions.
 -This may let us do away with the pre-allocation of specific DSP devices. We could perhaps use 
  an array or vector of direct objects for the DSPs instead of pointers (the objects would be the
  variants or unions). That would also increase the locality of the data used in the DSP chain. 
  -> good!
 -But we would perhaps have to commit to a certain maximum length of the DSP chain and the memory 
  usage of a variant DSP would be dictated by the maximum size of all of them which may be 
  wasteful. It would also be wasteful because the actual length may always be equal to the maximum
  length.
  -> bad!
-We could also have a general pool of pre-allocated memory and let the DSP chain use a std::vector
 with a custom allocator tapping into that pool.
 -That would require to write a custom memory allocator. Such an allocator could pre-allocate a 
  certain amount of memory (in bytes), then maybe split that into buckets of different sizes

Notes:
-Different effects have different requirements on the resources. At least two relevant resources
 to consider are CPU usage and RAM usage mapping to the classical time- and space-requirements in
 computer science. Here is a table of what I would expect from theoretical considerations:

   DSP Process                        CPU   RAM

   simple filter (biquad, etc.)        L     L     L: low, H: high, M: medium
   gain, pan, width, tremolo, wah      L     L
   naive distortion (waveshaping)      L     L
   flanger, vibrato                    L     M
   echo, delay                         L     H

   anti-aliased distortion             M     L
   ringmod, freq-shift, dynamics       M     L
   chorus, pitch-shift                 M     M
   dynamics with look-ahead            M     M
   adaptive filters (LMS, etc)         M     M
   simple reverb                       M     H

   complex filter (VA, ZDF, etc.)      H     L
   spectral stuff (FFT)                H     M
   good reverb, convolution            H     H

-Another thing to consider is the possible introduction of latency and/or the spikeyness of 
 CPU-load distribution due to block-based processing. Anything with look-ahead will introduce 
 latency. Pitch-shifters will also. FFT-based stuff also will, except for convolution, if a 
 zero-delay algo is used.

*/
