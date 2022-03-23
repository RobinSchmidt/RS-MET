namespace rosic {
namespace Sampler {

void Processor::addParameter(Opcode opcode)
{
  params.push_back(Parameter(opcode));
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

//=================================================================================================

LowFreqOsc::LowFreqOsc()
{
  type = OpcodeType::FreeLfo;
  params.reserve(3);                      // index
  addParameter(Opcode::lfoN_freq);        //   0
  //addParameter(Opcode::lfoN_amp);         //   1
  //addParameter(Opcode::lfoN_fade);
  //addParameter(Opcode::lfoN_delay);
  // ToDo: phase, wave, sync, ...
}

void LowFreqOsc::updateCoeffs(double fs)
{
  core.setup(params[0].mv(), 1.f, 0.f, 0.f, 0.f, (float)fs);
  dirty = false;
}

//=================================================================================================

Amplifier::Amplifier()
{
  type = OpcodeType::Amplifier;
  params.reserve(7);                      // index
  addParameter(Opcode::volumeN);          //   0
  addParameter(Opcode::panN);             //   1
  addParameter(Opcode::widthN);           //   2
  addParameter(Opcode::positionN);        //   3
  addParameter(Opcode::ampN_veltrack);    //   4
  addParameter(Opcode::ampN_keytrack);    //   5
  addParameter(Opcode::ampN_keycenter);   //   6

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
  float volume   = params[0].mv();
  float pan      = params[1].mv();
  float width    = params[2].mv();
  float position = params[3].mv();

  // Extract key/vel modifiers:
  float ampN_veltrack  = params[4].mv();
  float ampN_keytrack  = params[5].mv();
  float ampN_keycenter = params[6].mv();

  // Apply modifiers:
  volume += ampN_veltrack * 0.01f * 40 * log10f(127.f/(float)vel);
  volume += ampN_keytrack * ((float)key - ampN_keycenter);
  // Unit of veltrack is percent and the formula is dB = 20 log (127^2 / Velocity^2). Keytracking 
  // is adjusted in dB per key. See https://sfzformat.com/legacy/

  // Set up the core:
  core.setup(volume, pan, width, position);
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

void EffectPool::allocateEffects()
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

Processor* EffectPool::grabEffect(OpcodeType type)
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

void EffectPool::repositEffect(Processor* p)
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

//=================================================================================================

void ModulatorPool::allocateModulators()
{
  envGens.init(64);
  lowFreqOscs.init(64);
}

Processor* ModulatorPool::grabModulator(OpcodeType type)
{
  using OT = OpcodeType;
  Processor* p = nullptr;
  switch(type)
  {
  case OT::FreeEnv:   p = envGens.grabItem(); break;
  case OT::AmpEnv:    p = envGens.grabItem(); break;
  case OT::FilterEnv: p = envGens.grabItem(); break;
  case OT::PitchEnv:  p = envGens.grabItem(); break;

  case OT::FreeLfo:   p = lowFreqOscs.grabItem(); break;
  case OT::AmpLfo:    p = lowFreqOscs.grabItem(); break;
  case OT::FilterLfo: p = lowFreqOscs.grabItem(); break;
  case OT::PitchLfo:  p = lowFreqOscs.grabItem(); break;
  };
  return p;
  // ToDo: maybe consolidate the cases that return from envGens or lowFreqOscs into some sort of
  // if-statement that checks, if type is within some range
}

void ModulatorPool::repositModulator(Processor* p)
{
  using OT = OpcodeType;
  int i = -1;
  switch(p->getType())
  {
  case OT::FreeEnv:   i = envGens.repositItem(p); break;
  case OT::AmpEnv:    i = envGens.repositItem(p); break;
  case OT::FilterEnv: i = envGens.repositItem(p); break;
  case OT::PitchEnv:  i = envGens.repositItem(p); break;

  case OT::FreeLfo:   i = lowFreqOscs.repositItem(p);  break;
  case OT::AmpLfo:    i = lowFreqOscs.repositItem(p);  break;
  case OT::FilterLfo: i = lowFreqOscs.repositItem(p);  break;
  case OT::PitchLfo:  i = lowFreqOscs.repositItem(p);  break;
  }
  RAPT::rsAssert(i != -1, "Reposited processor was not in pool");
}




// preliminary - just delegate - todo: move actual implementations here

DspResourcePool::DspResourcePool()
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
  effectPool.allocateEffects();
}

Processor* DspResourcePool::grabEffect(OpcodeType type)
{
  return effectPool.grabEffect(type);
}

void DspResourcePool::repositEffect(Processor* p)
{
  effectPool.repositEffect(p);
}

void DspResourcePool::allocateModulators()
{
  modulatorPool.allocateModulators();
}

Processor* DspResourcePool::grabModulator(OpcodeType type)
{
  return modulatorPool.grabModulator(type);
}

void DspResourcePool::repositModulator(Processor* p)
{
  modulatorPool.repositModulator(p);
}





}} // namespaces

//=================================================================================================
/*

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
 -Maybe make also a struct for very basic 1-pole/1-zero filter. It can be realized by the biquad 
 (and by the other structures, too), but maybe it's more efficient to do it like that. I expect 
 that a simple 1st order lowpass is a quite common thing to use.
 -Maybe
 -Use the filter also for the equalizer opcode. No need to define a different class for that. Maybe
 extend sfz to support 4 instead of 3 bands when we later can realize 2 bands per filter...

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
-internal algo parameters should be gLL, gLR, gRL, gRR - the 4 gains for the channel-mix matrix
-user params should be: scale (linear overall scale factor), pan, width and as 3rd...dunno - maybe
 pos? whatever is left to determine the rest.


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
