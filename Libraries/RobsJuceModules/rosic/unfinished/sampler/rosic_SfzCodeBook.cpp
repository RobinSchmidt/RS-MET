namespace rosic { namespace Sampler {


ModulationRouting::ModulationRouting(OpcodeType modSrcType, int modSrcIndex, Opcode modTarget, 
  int modTargetIndex, float modDepth, ModMode modMode)
  : sourceType(modSrcType), sourceIndex(modSrcIndex), target(modTarget)
  , targetIndex(modTargetIndex), depth(modDepth), mode(modMode)
{
  targetType = SfzCodeBook::getInstance()->getOpcodeType(modTarget);
}


SfzCodeBook* SfzCodeBook::instance = nullptr;

SfzCodeBook::SfzCodeBook()
{
  // On construction, we build our database (maybe factor out):
  using OC = Opcode;
  using OF = OpcodeFormat;
  using SP = OpcodeType;
  using OU = OpcodeUnit;
  using OS = OpcodeSpec;
  opcodeEntries.resize((int)Opcode::NumTypes);

  // We need some abbreviations:
  auto add = [this](OC op, OF fmt, const char* name, float minVal, float maxVal, float defVal,
    SP dspType, OU unit, OS spec)
  { addOpcode(op, fmt, name, minVal, maxVal, defVal, dspType, unit, spec); };

  // Opcode value formatting:
  OF Nat = OF::Natural;       // natural numbers (starting at 0) -> indices, controller-numbers
  OF Int = OF::Integer;       // integer numbers 
  OF Flt = OF::Float;         // floating point numbers -> continuous parameters
  OF Txt = OF::String;        // strings -> choice parameters

  // SFZ specification:
  OS Sfz1  = OS::Sfz_1;
  OS Sfz1e = OS::Sfz_1_E;
  OS Sfz2  = OS::Sfz_2;
  OS Sfz2e = OS::Sfz_2_E;
  OS Aria  = OS::Aria;
  OS ARIAe = OS::Aria_E;

  // Player response constraints (aka "Input Control" in the sfz doc):
  SP dsp = OpcodeType::SamplePlayer;  // rename to type or tp or ot bcs not all types relate to DSP objects
  add(OC::Key,   Nat, "key",   0, 127,   0, dsp, OU::MidiKey, Sfz1);
  add(OC::LoKey, Nat, "lokey", 0, 127,   0, dsp, OU::MidiKey, Sfz1);
  add(OC::HiKey, Nat, "hikey", 0, 127, 127, dsp, OU::MidiKey, Sfz1);
  add(OC::LoVel, Nat, "lovel", 0, 127,   0, dsp, OU::RawInt,  Sfz1);
  add(OC::HiVel, Nat, "hivel", 0, 127, 127, dsp, OU::RawInt,  Sfz1);

  // Player playback settings:
  add(OC::Delay,  Flt, "delay",  0,        100, 0, dsp, OU::Seconds, Sfz1);
  add(OC::Offset, Nat, "offset", 0, 4294967296, 0, dsp, OU::Samples, Sfz1);
  //add(OC::Count,  Nat, "count",  0, 4294967296, 0, dsp, OU::RawInt,  Sfz1);

  add(OC::LoopMode, Txt, "loop_mode",  (float)LoopMode::Unknown + 1.f, 
    (float)LoopMode::numModes - 1.f, (float)LoopMode::no_loop, dsp, OU::Text, Sfz1);
  // In sfz, the default value depends on whether or not the sample file itself defines a loop,
  // if yes: continuous, if not: no loop

  add(OC::LoopStart, Nat, "loop_start", 0, 4294967296, 0, dsp, OU::Samples, Sfz1);
  // If not specified, the defined value in the sample will be used, if any, zero otherwise.

  add(OC::LoopEnd,   Nat, "loop_end",   0, 4294967296, 0, dsp, OU::Samples, Sfz1);
  // If not specified, the defined value in the sample will be used, if any.

  // ToDo:
  // -allow floating point numbers for loop-start and end (double precision)
  // -introduce opcode for loop-crossfading (in samples), maybe loop_xfade,
  // maybe another opcode could control the blend-function continuously between linear and 
  // constant power - maybe use a cubic approximation that fades the derivative at 1 from
  // 1 to 0 (for the fade-in function)...see also xf_keycurve...maybe use loop_xf

  // Player pitch:
  add(OC::Transpose,      Int, "transpose",        -127,   127,   0, dsp, OU::Semitones,  Sfz1);
  add(OC::Tune,           Int, "tune",             -100,   100,   0, dsp, OU::Cents,      Sfz1);
  add(OC::PitchKeyCenter, Nat, "pitch_keycenter",  -127,   127,  60, dsp, OU::MidiKey,    Sfz1); 
  add(OC::PitchKeyTrack,  Int, "pitch_keytrack",  -1200, +1200, 100, dsp, OU::CentPerKey, Sfz1); 
  // For pitch_keycenter, the spec actually says, -127..127 for the range but also C-1...G9 and 
  // C-1 would map to 0. So, that's probably a mistake in the sfz documentation and they actually 
  // mean 0..127. I assume this here but should check what other implementation do. In sfz files, 
  // keycenter can be given numerically or textually as e.g. c#2, so the parser should support 
  // midi-note to number conversion. ...or well...it actually may make sense to have keycenters 
  // well in the subaudio range...these are heavily oversampled. Although there is no midi-key 
  // lower than0, a sample's pitch keycenter is not subject to such a restriction...so yeah, we
  // allow the specified range


  // Filter:
  dsp = OpcodeType::Filter;
  add(OC::filN_type, Txt, "filN_type", (float)FilterType::Unknown + 1.f, 
    (float)FilterType::numTypes - 1.f, (float)FilterType::lp_12, dsp, OU::Text, Sfz1e); 
  // sfz default is lpf_2p - maybe rename our enum values to be consistent with sfz

  add(OC::cutoffN, Flt, "cutoffN", 0.f, 22050.f, 22050.f, dsp, OU::Hertz, Sfz1e);
  // Range is 0..fs/2, default is: filter disabled, so perhaps, the default should depend on the 
  // selected type: fs/2 for a lowpass, 0 for a highpass - figure out what sfz+ does
  // ...maybe switch the filter into bypass mode, if cutoff is set to zero in the PlaybackSetting

  add(OC::resonanceN, Flt, "resonanceN", 0.0f, +40.f, 0.f, dsp, OU::Decibels, Sfz1e);
  // In sfz, the resonance is adjusted in terms of the resonance gain in dB. If zero dB is the 
  // minimum, does that mean there is always some resonance? Because without resonance, the gain
  // at cutoff would be -3.01 dB, I think. Does the resonance parameter give the gain at the cutoff
  // freq or at the peak? -> figure out experimentally or by looking at other implementations. 
  // Maybe I should derive a formula for Q in terms of the resonance gain. The formula for the peak
  // gain would probably be more complicated.

  add(OC::filN_keytrack,  Int, "filN_keytrack",     0.f, 1200.f,  0.f, dsp, OU::CentPerKey, Sfz1);
  add(OC::filN_keycenter, Int, "filN_keycenter",    0.f,  127.f, 60.f, dsp, OU::MidiKey,     Sfz1);
  add(OC::filN_veltrack,  Int, "filN_veltrack", -9600.f, 9600.f,  0.f, dsp, OU::Cents,       Sfz1);
  // ToDo: allow floating point values and use -1200 as min for keytrack...actually, we do not need
  // to do anything for this - there is no enforcement of limits anywhere in the code anyway...but
  // maybe there should be? but maybe that can be switched off by another opcode like
  // param_limits=sfz1, param_limits=none, param_limits



  // Player amplifier:
  dsp = OpcodeType::Amplifier; 
  add(OC::amplitudeN, Flt, "amplitudeN",    0.f,  100.f, 100.f, dsp, OU::Percent,  ARIAe);
  add(OC::volumeN,    Flt, "volumeN",    -144.f,   +6.f,   0.f, dsp, OU::Decibels, Sfz1e);
  add(OC::panN,       Flt, "panN",       -100.f, +100.f,   0.f, dsp, OU::RawFloat, Sfz1e);
  add(OC::widthN,     Flt, "widthN",     -100.f, +100.f, 100.f, dsp, OU::Percent,  Sfz1e);
  add(OC::positionN,  Flt, "positionN",  -100.f, +100.f,   0.f, dsp, OU::Percent,  Sfz1e);

  add(OC::ampN_keytrack,  Flt, "ampN_keytrack",  -96.0f, 12.f, 0.f, dsp, OU::DecibelPerKey, Sfz1e);
  add(OC::ampN_keycenter, Flt, "ampN_keycenter",   0.0f,127.f,60.f, dsp, OU::MidiKey,       Sfz1e);
  add(OC::ampN_veltrack,  Flt, "ampN_veltrack", -100.0f,100.f, 0.f, dsp, OU::Percent,       Sfz1e);
  // Wait! The spec says that the default value for width is 0%. 
  // https://sfzformat.com/legacy/
  // But that's not a neutral value! It will monoize the samples by default? is that right or is 
  // that a typo? I have set it to 100% here because I think that's a much saner default.
  // Maybe, if sfz really turns out to have messed up the defaults, we should have an option to
  // switch to neutral defaults. Maybe via an opcode: default_values=neutral or something. Maybe
  // the codebook records need an additional field neutVal. Check behavior of sfz+. Load a stereo
  // sample maybe using sine-waves with different frequencies for left and right channel. It is 
  // also a bit inconvenient, that sfz prescribes a distinction between mono and stereo samples.
  // Further down the signal chain, a mono sample may be stereoized by some DSP along the way.
  // Maybe we should not treat the amplifier settings in the sfz way and instead use 
  // ampN_scale, ampN_pan, ampN_width, ampN_position...or some other parametrization


  // https://sfzformat.com/opcodes/amplitude



  // Equalizer:
  dsp = OpcodeType::Equalizer;
  add(OC::eqN_freq, Flt, "eqN_freq",   0.0f, 30000.f, 1000.f, dsp, OU::Hertz,    Sfz1e);
  add(OC::eqN_gain, Flt, "eqN_gain", -96.f,    +24.f,    0.f, dsp, OU::Decibels, Sfz1e);
  add(OC::eqN_bw,   Flt, "eqN_bw",     0.001f,   4.f,    1.f, dsp, OU::Octaves,  Sfz1e);
  // ToDo:
  // -The upper limit for equalizer center frequencies in sfz is defined to be 30 kHz. That would 
  //  fall above the Nyquist limit for 44.1 kHz sample-rate. I guess, we will need design formulas
  //  that allow this and act as a sort of high-shelver when the cutoff is above fs/2? Agai, we 
  //  need to figure out what other implementations do - by looking at the code where possible and
  //  by doing measurements of ths sfz+ (which i regard as reference implementation). 
  // -The lower limit of 0.001 for the bandwidth is rather low indeed so we need to check, if such
  //  narrow (i.e. high-Q) filters can actually be implemented in single precision. If not, we need
  //  to use double for the equalizer.

  // This is very very preliminary - don't use it yet to define actual instruments - its behavior
  // may be going to change:
  dsp = OpcodeType::WaveShaper;
  OS RsMet = OS::RsMet;
  add(OC::distortN_shape, Nat, "distortN_shape",  0.f,  0.f, 0.f, dsp, OU::RawInt,   RsMet); // not yet used
  add(OC::distortN_drive, Flt, "distortN_drive",  0.f,  8.f, 1.f, dsp, OU::RawFloat, RsMet);
  add(OC::distortN_dc,    Flt, "distortN_dc",    -1.f, +1.f, 0.f, dsp, OU::RawFloat, RsMet);
  // maybe allow negative values for drive like -8...+8...maybe that can be an additional scale
  // parameter - rename opcodes to distortN_shape, etc. cakewalk has some disto_... opcodes defined
  // -> figure out how they work - there's depth, tone, etc. maybe rename opdcodes to something
  // with waveshaper or wvshpr...hmm...naaah


  // ToDo: 
  // -PanLaw, introduce fil_bw for bandwidth parameter for bandpass...actually, we need to figure
  //  out, how the resonance parameter behaves for a bandpass
  // -maybe rename our enum values to map 1:1 to the sfz opcode names
  // -Our single precision floating point representation of integers will run into issues for very 
  //  large integers like in the "offset" opcode where the range goes all the way up to 2^32. Maybe 
  //  switch to double. Or maybe use a union of float32, int32 and uint32. Maybe it's pointless 
  //  trying to minimize the memory footprint of this - the opcodes are not accessed at sample-rate
  //  anyway, so using double should probably be fine
  // -verify everything (min,max,defaults,etc.) by comparing against the sfz spec
  // -try to make the calls shorter by using shorter name in the Opcode enum and maybe define 
  //  abbreviations for the units such as st for OT::Semitones
  // -maybe split the SamplePlayer opcodes into input-controls, amplitude, pitch, etc. as
  //  i done the spec...maybe name them PlayerPitch, PlayerAmp, PlayerResponseCtrl
  // -maybe keep a separate list of (yet) unsupported opcodes
  // -create similar lists for the filter types and other text parameters
  // -maybe factor out the different table creations into separate functions


  // Fill the lookup table with the filter types:
  using FT = FilterType;
  filterTypeEntries.resize((int)FT::numTypes);
  addFilterType(FT::lp_6,   "lpf_1p");
  addFilterType(FT::hp_6,   "hpf_1p");
  addFilterType(FT::lp_12,  "lpf_2p");
  addFilterType(FT::hp_12,  "hpf_2p");
  addFilterType(FT::bp_6_6, "bpf_2p");

  // Filter types:
  // SFZ 1: lpf_1p, hpf_1p, lpf_2p, hpf_2p, bpf_2p, brf_2p
  // SFZ 2: lpf_4p, hpr_4p, lpf_6p, hpf_6p, bpf_1p, brf_1p, apf_1p, pkf_2p, lpf_2p_sv, hpf_2p_sv, 
  //        bpf_2p_sv, brf_2p_sv, comb, pink.
  // https://sfzformat.com/legacy/
  // https://sfzformat.com/opcodes/
  // https://www.linuxsampler.org/sfz/
  // http://ariaengine.com/forums/index.php?p=/discussion/4389/arias-custom-opcodes/p1
  // ...it has 6-pole filters! :-O can we realize that with the current filter implementation 
  // without increasing its memory footprint? maybe using 3 equal biquads in DF2 or TDF1?
  // ...but maybe it would be a better idea to not use one class that does all filter types but
  // instead have for each filter topology an extra class. This will reduce the memory footprint
  // when only simple filters are used in a patch and opens the possibility to later include really
  // fancy filters without blowing up the memory footprint of patches which use only simple 
  // filters


  // Fixed modulators:
  dsp = OpcodeType::AmpLfo;
  add(OC::amplfo_freq,  Flt, "amplfo_freq",    0.0f,  20.f, 0.f, dsp, OU::Hertz,    Sfz1);
  add(OC::amplfo_depth, Flt, "amplfo_depth", -10.0f, +10.f, 0.f, dsp, OU::Decibels, Sfz1);
  // Or should the amplfo_depth be of type ModulationRouting, ModConnection? ..we'll see

  dsp = OpcodeType::FilterLfo;
  add(OC::fillfo_freq,  Flt, "fillfo_freq",       0.0f,   20.f, 0.f, dsp, OU::Hertz, Sfz1);
  add(OC::fillfo_depth, Flt, "fillfo_depth", -12000.f, 12000.f, 0.f, dsp, OU::Cents, Sfz1);

  dsp = OpcodeType::AmpEnv;
  add(OC::ampeg_delay,         Flt, "ampeg_delay",           0.f, 100.f,   0.f,     dsp, OU::Seconds,  Sfz1);
  add(OC::ampeg_start,         Flt, "ampeg_start",           0.f, 100.f,   0.f,     dsp, OU::Percent,  Sfz1);
  add(OC::ampeg_attack,        Flt, "ampeg_attack",          0.f, 100.f,   0.f,     dsp, OU::Seconds,  Sfz1);
  add(OC::ampeg_peak,          Flt, "ampeg_peak",            0.f, 100.f, 100.f,     dsp, OU::Percent,  RsMet); 
  add(OC::ampeg_hold,          Flt, "ampeg_hold",            0.f, 100.f,   0.f,     dsp, OU::Seconds,  Sfz1);
  add(OC::ampeg_decay,         Flt, "ampeg_decay",           0.f, 100.f,   0.f,     dsp, OU::Seconds,  Sfz1);
  add(OC::ampeg_sustain,       Flt, "ampeg_sustain",         0.f, 100.f, 100.f,     dsp, OU::Percent,  Sfz1);
  add(OC::ampeg_release,       Flt, "ampeg_release",         0.f, 100.f,   0.f,     dsp, OU::Seconds,  Sfz1);
  add(OC::ampeg_end,           Flt, "ampeg_end",             0.f, 100.f,   0.f,     dsp, OU::Percent,  RsMet);
  add(OC::ampeg_depth,         Flt, "ampeg_depth",           0.f, 100.f, 100.f,     dsp, OU::Percent,  RsMet);
  add(OC::ampeg_attack_shape,  Flt, "ampeg_attack_shape",  -20.f,  20.f,   0.f,     dsp, OU::RawFloat, Aria);
  add(OC::ampeg_decay_shape,   Flt, "ampeg_decay_shape",   -20.f,  20.f, -10.3616f, dsp, OU::RawFloat, Aria);
  add(OC::ampeg_release_shape, Flt, "ampeg_release_shape", -20.f,  20.f, -10.3616f, dsp, OU::RawFloat, Aria);
  // ToDo: Check the ranges for the shape parameters. It's not documented on sfzformat.com. Also, I think, 
  // they may be defined also in sfz2 not only in Aria, see here:
  //   https://sfzformat.com/opcodes/egN_shapeX
  // In such a case, in may be better to file them under the Sfz2 spec
  // not sure, if depth should be among the AmpEnv opcodes or among the routing opcodes

  dsp = OpcodeType::FilterEnv;
  add(OC::fileg_delay,         Flt, "fileg_delay",             0.f,   100.f,   0.f, dsp, OU::Seconds,  Sfz1);
  add(OC::fileg_start,         Flt, "fileg_start",             0.f,   100.f,   0.f, dsp, OU::Percent,  Sfz1);
  add(OC::fileg_attack,        Flt, "fileg_attack",            0.f,   100.f,   0.f, dsp, OU::Seconds,  Sfz1);
  add(OC::fileg_peak,          Flt, "fileg_peak",              0.f,   100.f, 100.f, dsp, OU::Percent,  RsMet); 
  add(OC::fileg_hold,          Flt, "fileg_hold",              0.f,   100.f,   0.f, dsp, OU::Seconds,  Sfz1);
  add(OC::fileg_decay,         Flt, "fileg_decay",             0.f,   100.f,   0.f, dsp, OU::Seconds,  Sfz1);
  add(OC::fileg_sustain,       Flt, "fileg_sustain",           0.f,   100.f, 100.f, dsp, OU::Percent,  Sfz1);
  add(OC::fileg_release,       Flt, "fileg_release",           0.f,   100.f,   0.f, dsp, OU::Seconds,  Sfz1);
  add(OC::fileg_end,           Flt, "fileg_end",               0.f,   100.f,   0.f, dsp, OU::Percent,  RsMet);
  add(OC::fileg_depth,         Flt, "fileg_depth",        -12000.f, 12000.f,   0.f, dsp, OU::Cents,    Sfz1);
  add(OC::fileg_attack_shape,  Flt, "fileg_attack_shape",    -20.f,    20.f,   0.f, dsp, OU::RawFloat, Aria);
  add(OC::fileg_decay_shape,   Flt, "fileg_decay_shape",     -20.f,    20.f,   0.f, dsp, OU::RawFloat, Aria);
  add(OC::fileg_release_shape, Flt, "fileg_release_shape",   -20.f,    20.f,   0.f, dsp, OU::RawFloat, Aria);
  // yes, the fileg decay and release shapes have 0 as default, different from th -10.3616 of the
  // ampeg. The same is true for the pitcheg decay/release shapes. Only the ameg has these "weird"
  // default values. I'm still not sure about the min/max values.
  // not sure, if depth should be among the FilterEnv opcodes or among the routing opcodes

  // Routable modulators:
  dsp = OpcodeType::FreeEnv;
  add(OC::adsrN_delay,         Flt, "adsrN_delay",           0.f, 100.f,   0.f, dsp, OU::Seconds,  RsMet);
  add(OC::adsrN_start,         Flt, "adsrN_start",           0.f, 100.f,   0.f, dsp, OU::Percent,  RsMet);
  add(OC::adsrN_attack,        Flt, "adsrN_attack",          0.f, 100.f,   0.f, dsp, OU::Seconds,  RsMet);
  add(OC::adsrN_peak,          Flt, "adsrN_peak",            0.f, 100.f, 100.f, dsp, OU::Percent,  RsMet); 
  add(OC::adsrN_hold,          Flt, "adsrN_hold",            0.f, 100.f,   0.f, dsp, OU::Seconds,  RsMet);
  add(OC::adsrN_decay,         Flt, "adsrN_decay",           0.f, 100.f,   0.f, dsp, OU::Seconds,  RsMet);
  add(OC::adsrN_sustain,       Flt, "adsrN_sustain",         0.f, 100.f, 100.f, dsp, OU::Percent,  RsMet);
  add(OC::adsrN_release,       Flt, "adsrN_release",         0.f, 100.f,   0.f, dsp, OU::Seconds,  RsMet);
  add(OC::adsrN_end,           Flt, "adsrN_end",             0.f, 100.f,   0.f, dsp, OU::Percent,  RsMet);
  add(OC::adsrN_attack_shape,  Flt, "adsrN_attack_shape",  -20.f,  20.f,   0.f, dsp, OU::RawFloat, RsMet);
  add(OC::adsrN_decay_shape,   Flt, "adsrN_decay_shape",   -20.f,  20.f,   0.f, dsp, OU::RawFloat, RsMet);
  add(OC::adsrN_release_shape, Flt, "adsrN_release_shape", -20.f,  20.f,   0.f, dsp, OU::RawFloat, RsMet);




  dsp = OpcodeType::FreeLfo;
  add(OC::lfoN_freq,  Flt, "lfoN_freq",    0.f,   20.f,   0.f, dsp, OU::Hertz,   Sfz2);
  add(OC::lfoN_amp,   Flt, "lfoN_amp",  -100.f, +100.f, 100.f, dsp, OU::Percent, RsMet);
  //add(OC::lfoN_amp,   Flt, "lfoN_amp",    -1.f,  +1.f, 1.f, dsp, OU::RawFloat, RsMet);
  // ToDo: figure out what lfoN_amplitude in sfz2 is supposed to do - maybe it does what we want to
  // do with amp here? Maybe we should use percent instead of RawFloat? I suppose, lfoN_amplitude 
  // is just the routing of lfoN to the amplitude. Also, we may want to set it up in percent for
  // consistency with the amplifier's amplitude parameter
  // maybe use "amount" instead of "amp"
  // implement phase, delay, fade, etc.

  // Free modulation routings:
  dsp = OpcodeType::LfoN_ParamX;
  add(OC::lfoN_volumeX, Flt, "lfoN_volumeX",    -10.0f,  10.f, 0.f, dsp, OU::Decibels, Sfz2e);
  //add(OC::lfoN_volumeX, Flt, "lfoN_amplitudeX",   0.0f,  20.f, 0.f, dsp, OU::RawFloat, Sfz2e);
  //add(OC::lfoN_volumeX, Flt, "lfoN_panX",         0.0f,  20.f, 0.f, dsp, OU::RawFloat, Sfz2e);
  // todo: veriyf ranges, units, defaults, etc.

  // SFZ2 has: lfoN_freq, lfoN_amplitude, lfoN_volume, lfoN_pan, lfoN_cutoff, lfoN_cutoff2
  // we want:  lfoN_amplitudeX, lfoN_volumeX, lfoN_cutoffX, ...

  dsp = OpcodeType::MidiCtrl;
  add(OC::controlN_index,   Int, "controlN_index",   0.f, 127.f, 0.f, dsp, OU::Index,    RsMet);
  add(OC::controlN_neutral, Int, "controlN_neutral", 0.f, 127.f, 0.f, dsp, OU::Index,    RsMet);


  dsp = OpcodeType::Setup;
  add(OC::set_ccN,   Flt, "set_ccN",   0.f, 127.f, 0.f, dsp, OU::RawFloat, Sfz2);  // default, initial value
  add(OC::label_ccN, Txt, "label_ccN", 0.f,   0.f, 0.f, dsp, OU::Text,     Aria);  



  // Hardwired modulation routings:
  dsp = OpcodeType::HardwiredModRouting;




  // We also want LFO phases and we want to be able able to specify the LFO frequency in 
  // temp-synced units - maybe an opcode lfoN_unit=hertz (cycles per second), cycles per beat or 
  // maybe lfoN_sync = true/false. Or maybe allow the user to optionally give a unit like:
  // like lfoN_freq=20_Hz ..here, _Hz would be optional because Hz is the default unit for LFO 
  // freqs but it could also be a different unit like 20_perBeat or 20_perMeasure. This feature 
  // could also be used to allow the user to specify other opcodes in different units like 
  // cutoff=69_key instead of cutoff=440_Hz, egN_attack=200_ms instead of egN_attack=0.02 (which 
  // would default to egN_attack=0.02_s), etc.
  // Maybe go to KVR to talk about such proposals. Make a thread:
  // "Extending SFZ: existing specifications, new proposals, general discussion"
  // Maybe let's call the extended format SFZ++ :-)


  int dummy = 0;
}

template<class T>
inline void rsEnsureSize(std::vector<T>& v, size_t s)
{
  if(v.size() < s)
    v.resize(s);
} // maybe move to rapt
void SfzCodeBook::addOpcode(Opcode op, OpcodeFormat type, const std::string& sfzStr,
  float minVal, float maxVal, float defVal, OpcodeType dspType, OpcodeUnit unit, OpcodeSpec spec)
{
  int i = (int)op;
  rsEnsureSize(opcodeEntries, size_t(i+1));
  opcodeEntries[i] = 
    OpcodeEntry({ op, type, sfzStr, minVal, maxVal, defVal, dspType, unit, spec });
  // Actually, storing the "op" is redundant because it's implicitly given by the array index, so
  // maybe remove that field...but maybe it's useful in other contexts
}

void SfzCodeBook::addFilterType(FilterType type, const std::string& sfzStr)
{
  int i = (int)type;
  rsEnsureSize(filterTypeEntries, size_t(i+1));
  filterTypeEntries[i] = FilterTypeEntry({ type, sfzStr });
}

bool SfzCodeBook::isFilterRelated(Opcode op) const
{
  RAPT::rsError("Not yet correctly implemented");
  return op == Opcode::cutoffN || op == Opcode::resonanceN || op == Opcode::filN_type;
  // ...these are not all - there are actually many more! look up, which of these need special
  // treatment with regard to interpreting absence of a number as 1. Make sure the filter-related
  // opcodes have contiguous indices and use >= and <= comparison here.

  // Why do we actually need that function? May it be obsolete, once we handle indexed opcodes 
  // generally?
}

OpcodeType SfzCodeBook::getOpcodeType(Opcode op)
{
  if((int)op < 0 || (int)op >= (int)opcodeEntries.size()) {  // todo: use if(!isValidOpcode(op))
    RAPT::rsError("Unknown opcode in SfzCodeBook::getOpcodeType");
    return OpcodeType::Unknown; 
  }
  return opcodeEntries[(int)op].dsp;  // ToDo: rename the "dsp" field to "type"
}

OpcodeFormat SfzCodeBook::getOpcodeFormat(Opcode op)
{
  if((int)op < 0 || (int)op >= (int)opcodeEntries.size()) {  // todo: use if(!isValidOpcode(op))
    RAPT::rsError("Unknown opcode in SfzCodeBook::getOpcodeFormat");
    return OpcodeFormat::Unknown; 
  }
  return opcodeEntries[(int)op].format;
}

float SfzCodeBook::opcodeMinValue(Opcode op)
{
  if(!isValidOpcode(op)) {
    RAPT::rsError("Unknown opcode in SfzCodeBook::opcodeMinValue");
    return 0.f;  }
  return opcodeEntries[(int)op].minVal;
}

float SfzCodeBook::opcodeMaxValue(Opcode op)
{
  if(!isValidOpcode(op)) {
    RAPT::rsError("Unknown opcode in SfzCodeBook::opcodeMaxValue");
    return 0.f;  }
  return opcodeEntries[(int)op].maxVal;
}

float SfzCodeBook::opcodeDefaultValue(Opcode op, int index)
{
  if((int)op < 0 || (int)op >= (int)opcodeEntries.size()) {  // todo: use if(!isValidOpcode(op))
    RAPT::rsError("Unknown opcode in SfzCodeBook::opcodeDefaultValue");
    return 0.f;
  }

  // For certain opcodes, the default-value depends on the index:
  using OC = Opcode;
  switch(op)
  {
  case OC::eqN_freq:
  {
    switch(index)
    {
    case 1: return   50.f;
    case 2: return  500.f;
    case 3: return 5000.f;
    }
  }
  }

  // ToDo: for loop_mode, the default value actually depends on whether or not the sample-file 
  // defines a loop...that should probably be handled by the caller as special case

  // For all others, we return the stored default value:
  return opcodeEntries[(int)op].defVal;
}

ModMode SfzCodeBook::opcodeDefaultModMode(Opcode op)
{
  using MM = ModMode;
  using OC = Opcode;

  switch(op)
  {
    // Filter:
  case OC::cutoffN:    return MM::cents;
  case OC::resonanceN: return MM::absolute;  // check, if this is right

    // Amplifier:
  case OC::volumeN:    return MM::absolute;
  case OC::panN:       return MM::absolute;  // check, if this is right
  case OC::amplitudeN: return MM::absolute;
  }

  RAPT::rsError("Unknown opcode in opcodeDefaultModMode");
  return MM::unknown;

  // Maybe all the MM::absolute should be subsumed in the default branch. But maybe not because it
  // would thwart our assertion which is very useful during development.
}


std::string SfzCodeBook::opcodeToString(Opcode op, int index, bool withIndex1) const
{
  if((int)op < 0 || (int)op >= (int)opcodeEntries.size()) {  // todo: use if(!isValidOpcode(op))
    RAPT::rsError("Unknown opcode in SfzCodeBook::opcodeToString");
    return dummyString; }

  //RAPT::rsAssert(index > 0 || index == -1, "Invalid index in SfzCodeBook::opcodeToString");
  // Either index is actually a valid index in which case must be a positive natural number or 
  // indexing doesn't apply to the given opcode in which case we use -1 as code to indicate this
  // situation.

  // Update 2023/05/22:
  RAPT::rsAssert(index >= -1, "Invalid index in SfzCodeBook::opcodeToString");
  // We may actually need to allow N to be 0 as well because the set_ccN opcode may be set_cc0
  // because zero is actually a valid midi controller index

  if(index != -1) {                              // if we are dealing with an indexed opcode...
    std::string s = opcodeEntries[(int)op].text; //   retrieve opcode template (e.g. eqN_freq)
    rsReplace(s, "N", std::to_string(index));    //   replace placeholder "N" with actual index
    if(!withIndex1)
      makeExplicitIndexImplicit(s);              //   turn e.g. fil1_type into fil_type
    return s; }
  else
    return opcodeEntries[(int)op].text;
}

//-------------------------------------------------------------------------------------------------

Opcode SfzCodeBook::stringToOpcode(const std::string& strIn, int* index)
{
  // Pre-process the opcode string to handle indexed opcodes. In such a case, we need to translate
  // the concrete opcode string containing the index into its corresponding opcode template, i.e. 
  // replace the index with the placeholder "N" and we also need to figure out, what that index is.
  // For example, the string eq13_freq should be turned into eqN_freq and we want to extract the 13
  // as integer:
  std::string str = strIn;              // maybe avoid this by taking parameter by value
  makeImplicitIndexExplicit(str);       // changes e.g. fil_type into fil1_type
  *index = getIndexAndReplaceByN(str);

  // Figure out, which opcode (or opcode-template) it is:
  for(int i = 0; i < opcodeEntries.size(); i++)
    if(opcodeEntries[i].text == str)
      return opcodeEntries[i].op;    // op should be equal to i

  //RAPT::rsError("Unknown opcode in SfzCodeBook::stringToOpcode");
  // nope - we need to allow that to support the mod-routing opcodes on a higher level

  return Opcode::Unknown;

  // This lookup has currently linear complexity in the number of opcodes. Maybe bring this down 
  // to at most O(log(N)) by maintaining a map of indices into the opcodeEntries array that is 
  // sorted lexicographically according to the opcode string. I'm not sure, if it's worth it 
  // though. This is called only on patch loading and maybe it's fast enough as is. We'll see.
  // Another option could be to use a hash-table as in std::map.

  // ToDo:
  // -We need a way to return the index, too! ...done?
  // -There are opcodes like pitcheg_vel2hold where the 2 is actually not an index!. In this case, 
  //  the 2 should not be removed. Maybe whenever a 2 is prefixed by vel, we should not remove it.
  //  This is a messy business with all sorts of special cases and needs thorough unit tests!
  // -For the opcodes that contain an index (or two), we perhaps need to preprocess the string to
  //  remove it or to replace it by "N". In sfz 2, there are also opcodes with two indices, like
  //  egN_timeX. Maybe first scan through the string from the beginning and the first number that
  //  is found is replaced by N, then do it again and when another number is found, replace it by 
  //  X. Maybe the 2nd scan could also go backward. 
}

/*
int SfzCodeBook::stringToIndex(const std::string& str)
{
  return -1;  // preliminary

  // This function should, for example, extract the 74 in cutoff_cc74. For eg3_time5, it should 
  // extract the 3. Maybe rename it to stringToIndexN and have another function stringToIndexX that
  // would extract the 5 in the 2nd case. N and X are used as placeholders here:
  // https://www.linuxsampler.org/sfz/
  // It seems like the caller (SfzInstrument::getSettingFromString) knows the opcode represented by
  // str, so if that's helpful, it could be included as parameter for the function. Perhaps that's
  // useful to implement the special rules for the filter-related params (to interpret cutoff as 
  // cutoff1, fil_type as fil1_type etc.). Maybe, we should have a function isFilterRelated
}
// may be obsolete because stringToOpcode now extracts both, the opcode-template and the index, so
// this may not be needed anymore
*/

const std::string& SfzCodeBook::filterTypeToString(FilterType ft)
{
  if((int)ft < 0 || (int)ft >= (int)filterTypeEntries.size()) {
    RAPT::rsError("Unknown type in SfzCodeBook::filterTypeToString");
    return dummyString;  }
  return filterTypeEntries[(int)ft].sfzStr;
}

FilterType SfzCodeBook::stringToFilterType(const std::string& str)
{
  for(int i = 0; i < filterTypeEntries.size(); i++)
    if(filterTypeEntries[i].sfzStr == str)
      return filterTypeEntries[i].typeId;    // should be equal to i
  RAPT::rsError("Unknown type in SfzCodeBook::stringToFilterType");
  return FilterType::Unknown;
}
// this code is repetitive! try to refactor! ...or use an implementation similar to the one for 
// LoopMode

LoopMode SfzCodeBook::stringToLoopMode(const std::string& str)
{
  if(str == "no_loop")         return LoopMode::no_loop;
  if(str == "loop_continuous") return LoopMode::loop_continuous;
  if(str == "one_shot")        return LoopMode::one_shot;
  return LoopMode::Unknown;
  // ToDo: use a switch statement
}

std::string SfzCodeBook::loopModeToString(LoopMode lm)
{
  if(lm == LoopMode::no_loop)         return "no_loop";
  if(lm == LoopMode::loop_continuous) return "loop_continuous";
  if(lm == LoopMode::one_shot)        return "one_shot";
  return "unknown";
  // ToDo: use a switch statement
}

std::string SfzCodeBook::modSourceToString(OpcodeType sourceType, int index)
{
  using OT = OpcodeType;
  std::string tmp;

  switch(sourceType)
  {
  case OT::FreeEnv:   tmp = "adsr";    break;
  case OT::AmpEnv:    tmp = "ampeg";   break;
  case OT::FilterEnv: tmp = "fileg";   break;

  case OT::FreeLfo:   tmp = "lfo";     break;
  case OT::AmpLfo:    tmp = "amplfo";  break;
  case OT::FilterLfo: tmp = "fillfo";  break;

  case OT::MidiCtrl:  tmp = "control"; break;   // new, needs tests

  default:
  {
    RAPT::rsError("Unknown type of modulation source.");
    return "";
  }
  }

  tmp += to_string(index);
  return tmp;
}

std::string SfzCodeBook::modTargetToString(OpcodeType type, int index, Opcode opcode)
{
  // Maybe the type parameter is superfluous. If so, remove it. The whole function my be 
  // superfluous and client code can just use opcodeToString instead. At the moment, we simply 
  // delegate to that anyway. But maybe later, we may need to do something more complex?

  return opcodeToString(opcode, index);
}

std::string SfzCodeBook::modDepthToString(float depth, ModMode mode, Opcode opcode)
{
  std::string s = to_string(depth);

  // ToDo: 
  // If the mode is not the default modulation mode for the given target opcode, append a 
  // unit suffix to the number string s like ct, Hz, %, st, oct. For example, for the cutoffN
  // opcode, the default mode is cents, but if the user wants to modulate cutoff in absolute Hz or
  // in percent, the value needs a unit suffix. We then also need to implement parsing of these
  // suffixes. Such suffixes may be applied to other opcodes, too (not only to modulation 
  // routings).

  return s;
}

OpcodeType SfzCodeBook::stringToModSource(const std::string& str, int* index)
{
  int i = lastNonDigit(str);

  std::string srcStr = str.substr(0,   i+1);
  std::string idxStr = str.substr(i+1, str.length()-i-1);
  // Test if this works also if there is no number at all. In this case, the idxStr should be empty

  if(idxStr.empty())
    *index = 1;
  else
    *index = parseNaturalNumber(idxStr, 0, (int)idxStr.length()-1);

  using OT = OpcodeType;
  if(srcStr == "adsr")     return OT::FreeEnv;
  if(srcStr == "ampeg")    return OT::AmpEnv;
  if(srcStr == "fileg")    return OT::FilterEnv;
  if(srcStr == "pitcheg")  return OT::PitchEnv;
  if(srcStr == "lfo")      return OT::FreeLfo;
  if(srcStr == "amplfo")   return OT::AmpLfo;
  if(srcStr == "fillfo")   return OT::FilterLfo;
  if(srcStr == "pitchlfo") return OT::PitchLfo;
  if(srcStr == "control")  return OT::MidiCtrl;   // new, needs tests

  RAPT::rsError("Unknown string for modulation source");
  return OpcodeType::Unknown;
}

float SfzCodeBook::stringToModDepth(const std::string& str, ModMode* modMode, Opcode target)
{
  int i = suffixStart(str);

  std::string numStr = str.substr(0, i);
  std::string sfxStr = str.substr(i, str.length()-i);
  // Test if this works also if there is no suffix at all. In this case, the sfxStr should be empty

  if(sfxStr.empty())
    *modMode = opcodeDefaultModMode(target);
  else
    *modMode = stringToModMode(sfxStr);

  return rsStringToFloat(numStr);
}

ModMode SfzCodeBook::stringToModMode(const std::string& str)
{
  using MM = ModMode;

  if(str == "ct") return MM::cents;
  // ...more to come


  RAPT::rsError("Unknown string for modulation mode");
  return MM::unknown;
}


std::string SfzCodeBook::valueToString(Opcode op, float val)
{
  switch(op)
  {
  case Opcode::filN_type: return filterTypeToString((FilterType)(int)val);
  case Opcode::LoopMode:  return loopModeToString(  (LoopMode)  (int)val);
  default: return to_string(val);
  }

  /*
  if(op == Opcode::filN_type)
  {
    FilterType i = (FilterType)(int)val;
    return filterTypeToString(i);
  }
  if(op == Opcode::LoopMode)
  {
    LoopMode i = (LoopMode)(int)val;
    return loopModeToString(i);
  }
  */


  //return to_string(val);

  // Maybe use custom string conversion functions because the std::to_string just uses a 
  // fixed number of 6 decimal digits after the point. Maybe that's suitable, but maybe not:
  // https://www.cplusplus.com/reference/string/to_string/
  // ...well, i think, it's not suitable for int params, but we may convert to int. I think, a 
  // fixed number (maybe 8 or 9..whatever number ensures lossless roundtrips) of total decimal 
  // digits is better
  // Maybe use a switch statement later when we have more such special cases
}



float SfzCodeBook::stringToValue(Opcode op, const std::string& str)
{
  switch(op)
  {
  //case Opcode::Unknown:   return 0.f;
    // Workaround to fix crashes with malformed sfz files. Actually, we should fall through to 
    // default and rsStringToFloat should return 0 when the string contains nonsense but that 
    // doesn't work yet because at the moment rsStringToFloat just calls stof which can't handle
    // nonsense gracefully.

  case Opcode::filN_type: return (float) stringToFilterType(str);
  case Opcode::LoopMode:  return (float) stringToLoopMode(str);
  default: return rsStringToFloat(str);
  //default: return std::stof(str);
  }
  // in the default path, we must first check, if the string represents a number - implement an 
  // rsStringToDouble function (may already exist), that is safe: it should return a sane default 
  // value when the string is not parsable as a number (the default value should be passed in and
  // default to zero)

  /*
  if(op == Opcode::filN_type)
  {
    FilterType ft = stringToFilterType(str);
    return (float) ft;
  }
  return std::stof(str);
  */
}

std::string SfzCodeBook::settingToString(const PlaybackSetting& s)
{
  //return std::string();

  using PST = Opcode;
  PST    op = s.getOpcode();
  float val = s.getValue();
  int   idx = s.getIndex();
  return opcodeToString(op, idx) + std::string("=") + valueToString(op, val);
}

std::string SfzCodeBook::modRoutingToString(const ModulationRouting& r)
{
  //return std::string();

  // Retrieve data from routing object:
  using OT = OpcodeType;
  using OC = Opcode;
  using MM = ModMode;
  OT  srcTp  = r.getSourceType();
  OT  tgtTp  = r.getTargetType();
  OC  tgtOp  = r.getTargetOpcode();
  int srcIdx = r.getSourceIndex();
  int tgtIdx = r.getTargetIndex();
  MM  mode   = r.getMode();

  // Handle special cases:
  if(srcTp == OT::AmpLfo && tgtOp == OC::volumeN && mode == MM::absolute)
  {
    RAPT::rsAssert(srcIdx == 1, "There should be at most one AmpLfo");
    return "amplfo_depth=" + rsFloatToString(r.getDepth());
  }
  if(srcTp == OT::AmpEnv && tgtOp == OC::amplitudeN && mode == MM::absolute)
  {
    // The ampeg is always implicitly routed to the last Amplifier in the chain.
    // Maybe an additional constraint should be that tgtIdx is the index of the last Amplifier such
    // that routings of the ampeg_ to other amplifiers will be stored. Currently we may throw away
    // too many ampeg routings
  }
  if(srcTp == OT::FilterEnv && tgtOp == OC::cutoffN && mode == MM::cents)
  {
    RAPT::rsAssert(srcIdx == 1, "There should be at most one FilterEnv");
    // ...because in the sfz spec, the fileg_ opcodes are not indexed, i.e. there's no such thing
    // as fileg2_cutoff3 or anything like that. There is just fileg_depth which we interpret as
    // fileg1_cutoff1 at the moment. Maybe later, we should interpret it as fileg1_cutoffAll but
    // for that, we will first generally have to support such an "All" syntax.

    RAPT::rsAssert(tgtIdx == 1);
    // Hmm...we currently interpret the fileg_depth opcode as applying to the first cutoff, i.e. 
    // the cutoff of the first filter in the chain. But maybe it should apply to all filters? But 
    // then we should probably only store it once in the sfz string...hmm...

    return "fileg_depth=" + rsFloatToString(r.getDepth()); 
  }
  // ToDo: Add other special cases here one by one. For example, fillfo, pitcheg, amplfo, etc.
  // Maybe factor out a function for handling the special cases, returning true, if the case was 
  // handled, false if not and call it here...or maybe not...

  // Handle the general case:
  std::string tmp;
  tmp  = modSourceToString(srcTp, srcIdx) + "_";
  tmp += modTargetToString(tgtTp, tgtIdx, tgtOp) + "=";
  tmp += modDepthToString(r.getDepth(), mode, tgtOp);
  return tmp;
}

SfzCodeBook* SfzCodeBook::getInstance()
{
  RAPT::rsAssert(instance != nullptr);
  // Client code is supposed to explicitly create the singleton instance using createInstance() 
  // before using it. It should also clean up by calling deleteInstance(), when the object is not 
  // needed anymore. We need this explicit lifetime management (in particular, the clean up) of 
  // the singleton to prevent false positives from the memory leak checker. Well, it's actually
  // a valid positive - the (GoF) textbook version of the pattern doesn't do any clean up. The book
  // "Pattern Hatching" addresses this omission.

  if(instance == nullptr)  // ...yeah, ok - just in case - defensive programming. But it's really 
    createInstance();      // cleaner to do an explicit creation somewhere before usage.
  return instance;
}

void SfzCodeBook::createInstance()
{
  RAPT::rsAssert(instance == nullptr);
  // Don't create a new instance before deleting the old one. That's a memory leak because the 
  // instance pointer will be overwritten and the old object to which it previously pointed will
  // never be deleted.

  instance = new SfzCodeBook;
}

void SfzCodeBook::deleteInstance()
{
  delete instance;
  instance = nullptr;
}

int SfzCodeBook::getIndexAndReplaceByN(std::string& str) const
{
  int is = firstDigit(str);                    // start position of first found number
  if(is == -1) return -1;                      // str does not contain an index
  int ie = lastDigitInSeq(str, is);            // end position of the number
  int num = parseNaturalNumber(str, is, ie);   // figure out the number
  str.replace(is, ie-is+1, "N");               // replace number by placeholder "N"
  return num;                                  // return the number
}
// move to SamplerTools

void SfzCodeBook::makeImplicitIndexExplicit(std::string& str) const
{
  if(     str == "fil_type")      str = "fil1_type";
  else if(str == "cutoff")        str = "cutoff1";
  else if(str == "resonance")     str = "resonance1";
  else if(str == "fil_keytrack")  str = "fil1_keytrack";
  else if(str == "fil_keycenter") str = "fil1_keycenter";
  else if(str == "fil_veltrack")  str = "fil1_veltrack";

  else if(str == "volume")        str = "volume1";
  else if(str == "pan")           str = "pan1";
  else if(str == "width")         str = "width1";
  else if(str == "position")      str = "position1";
}
void SfzCodeBook::makeExplicitIndexImplicit(std::string& str) const
{
  if(     str == "fil1_type")      str = "fil_type";
  else if(str == "cutoff1")        str = "cutoff";
  else if(str == "resonance1")     str = "resonance";
  else if(str == "fil1_keytrack")  str = "fil_keytrack";
  else if(str == "fil1_keycenter") str = "fil_keycenter";
  else if(str == "fil1_veltrack")  str = "fil_veltrack";

  else if(str == "volume1")        str = "volume";
  else if(str == "pan1")           str = "pan";
  else if(str == "width1")         str = "width";
  else if(str == "position1")      str = "position";
}
// This is stupid! Maybe we can do something more clever later. Maybe in the constructor, fill 
// two parallel arrays of strings, one with the parameter name with the index included and one 
// without (at corresponding positions). The addition of a parameter to the array could be done
// by a helper function that could be called like addImplicitlyIndexedParam("volume1"). The 
// function would add the string as is into one of the arrays and a version with the 1 removed
// into the other. Here in this function, we could then just loop through one array and if we 
// hit a match, return the corresponding string from the other. Still not pretty but much less 
// boilerplate to write. The function name should be abbreviated. We could even get away with 
// storing just one array of the strings without the index along with an integer for the 
// position where it should be inserted or removed...at the cost of more complex code for 
// finding a match.







bool SfzCodeBook::hasImplicitFirstGroup(const std::string& code)
{
  // An implicit first group definition (starting at the begin of the string) occurs, if the first
  // <region> header appears before the first <group> header.

  size_t firstGroupStart  = code.find("<group>");
  if(firstGroupStart == string::npos)
    return true;
  // No <group> header was found at all. In this case, there is just one group and that group is
  // indeed implicit.

  size_t firstRegionStart = code.find("<region>");
  if(firstRegionStart == string::npos)
    return false;
  // There is a <group> header but no <region> header. In this case, we have an explicit but 
  // empty (first) group.

  return firstGroupStart > firstRegionStart;
  // When both <group> and <region> headers actually do exist, we need to compare their positions. 
  // If the first region header comes befeore the first group header, the first group is implicit.
}

void SfzCodeBook::findGroup(const std::string& code, int groupIndex, 
  int* startIndex, int* endIndex)
{
  *startIndex = -1;
  *endIndex   = -1;
  if(code.empty())
    return;

  std::string pattern = "<group>";       // Search pattern
  size_t L = pattern.length();
  int foundIndex = -1;
  size_t start = 0;

  // Handle special case when the first group header is missing and therefore implicit:
  if(hasImplicitFirstGroup(code)) {
    *startIndex = 0;
    foundIndex  = 0; }

  // Find the start:
  while(foundIndex < groupIndex) {
    start = code.find(pattern, start);
    if(start != string::npos) {          // npos is returned when no match was found
      foundIndex++;
      *startIndex = (int) start;
      start += L;
      if(foundIndex == groupIndex)
        break; }
    else {
      *startIndex = -1;
      return;  }}           // No group with given index could be found in the given sfz code

  // Find the end. The end is defined to be either the last character in the string or the 
  // character immediately before the subsequent <group> header:
  start = code.find(pattern, start);
  if(start != string::npos)
    *endIndex = (int) start - 1;
  else
    *endIndex = (int) code.length() - 1;
}

void SfzCodeBook::findRegion(const std::string& code, int regionIndex, 
  int searchStart, int searchEnd, int* startIndex, int* endIndex)
{
  *startIndex = -1;
  *endIndex   = -1;
  if(code.size() < searchStart)  // maybe <=?
    return;

  std::string pattern = "<region>";       // Search pattern
  size_t L = pattern.length();
  int foundIndex = -1;
  size_t start = (size_t) searchStart;

  // Find the start:
  while(foundIndex < regionIndex && (int) start <= searchEnd) {
    start = code.find(pattern, start);
    if(start != string::npos) {          // npos is returned when no match was found
      foundIndex++;
      *startIndex = (int) start;
      start += L;
      if(foundIndex == regionIndex)
        break; }
    else {
      *startIndex = -1;
      return;  }}  

  // Handle case when we left the loop due to the 2nd loop condition. In this case the caller is
  // searching for a region with an index higher than present in the code:
  if(foundIndex < regionIndex) {
    *startIndex = -1;
    return; }
  // ToDo: document, why we don't have/need a similar handling in findSfzGroup. I think, it's 
  // because there, the 2nd loop condition doesn't exist

  // Find the end. The end is defined to be either the searchEnd or the character immediately 
  // before the subsequent <region> header:
  start = code.find(pattern, start);
  if(start != string::npos)
    *endIndex = (int) start - 1;
  else
    *endIndex = searchEnd;

  // ToDo:
  // -Maybe add some assertions to the beginning, namely that the searchEnd is not beyond the end
  //  of the code string and that no <group> header appears between searchStart and searchEnd.
}
// needs tests
// Can we get rid of the code duplication? Maybe like this:
// -Pass the search pattern as argument
// -Let searchstart/End be optional parameters that defualt to 0 and code.length()-1
// -I think the conditionals:
//    if(code.size() < searchStart) 
//    while(foundIndex < regionIndex && (int) start <= searchEnd)
//  in findSfzRegion are compatible with those in findSfzGroup. 
// -The if(hasImplicitFirstGroup(code)) in findSfzGroup could be problematic, though



void SfzCodeBook::findOpcode(const std::string& code, Opcode opcode, int opcodeIndex,
  int searchStart, int searchEnd, int* startIndex, int* endIndex)
{
  *startIndex = -1;
  *endIndex   = -1;

  //return;  // preliminary


  // If the index is 1, we have in some cases to search for two alternative search patterns, e.g.
  // "cutoff" or "cutoff1". It's a bit messy but needs to be done to support both syntaxes.
  std::string ptn1 = opcodeToString(opcode, opcodeIndex, true);  // This has the optional index 1
  std::string ptn2 = opcodeToString(opcode, opcodeIndex, false); // This hasn't
  size_t L1 = ptn1.length();
  size_t L2 = ptn2.length();
  size_t end = (size_t)searchEnd;



  // Returns true, iff the given position in the code is located within a comment:
  auto isInComment = [&](const std::string& code, size_t pos)
  {
    while(pos > searchStart)
    {
      char c = code[pos];
      if(c == '/')
        return true;
      if(c == '\n')
        return false;
      --pos;
    }
    return false;
    // In sfz, everything after a '/' on a given line is a comment, so we may detect a position
    // within a comment by going leftward. If we encounter a '/' before we encounter a line-break,
    // pos is indeed in a comment. If, on the other hand, we encounter a line-break before a '/', 
    // we may conclude that ps is not within a comment because the above way is the only way to
    // write comments in sfz.
  };
  // Maybe factor this out. This may be needed in findGroup/findRegion, too. SFZ authors may want
  // to comment out groups or regions. Our code section finder is not yet prepared to weed out
  // commented groups and regions.

  // Returns true, iff the given position in the code is located within a right-hand-side of an 
  // assignment:
  auto isInAssignment = [&](const std::string& code, size_t pos)
  {
    while(pos > searchStart)
    {
      char c = code[pos];
      if(c == '=')
        return true;
      if(c == '\n' || c == ' ' || c == '>' || c == '\t')
      //if(c == '\n' || c == '>' || c == '\t') // See comment below why we don't test against ' ' 
        return false;
      --pos;
    }
    return false;
  };
  // I think, testing against ' ' may break it when blank spaces occur within a filename, which is 
  // allowed according to the sfz spec, see: https://sfzformat.com/opcodes/sample which says:
  //   "...names with blank spaces and other special characters (excepting the = character) are 
  //    allowed..."
  // OK - so we may use a '=' character that appears somewhere to the left to identify a RHS
  // of an assignment, i.e. the "if(c == '=') return true;" segment should be ok. But the allowance
  // of blank spaces makes the other if-statement a bit more difficult. We should not actually 
  // conclude that we are not within a right-hand-side, as soon as we see a ' '. ...so what can we 
  // do instead? Hmm...I think, we can still immediately return false as soon as we see any of the 
  // other characters ('\n', '>', 't') but when we encounter a ' ', we should not immediately jump
  // to the conclusion that this is not a RHS but instead need some more complex test to figure 
  // out, if the space is part of a RHS freeform-string. Maybe introduce a 3rd if-statement like
  //   if(c == ' ' !isSpaceInRhsString(code, pos) ) return false
  // or something like that. Allowing spaces in RHS strings is bit of a burden here, but I actually
  // think, it's a good thing to allow this, especially when we later want to introduce a formula
  // opcode (in which spaces really should be allowed for nicer formatting). Oh - actually, we want
  // to allow the '=' in formulas as well...hmmm....maybe we should require strings to be entered
  // in quotes like "This is a string", except for sample-names where the quotes are optional for
  // historic reasons...but what if the string itself contains quotes?
  // Maybe factor out both functions isInComment and isInAssignment and test them on their own.
  // Actually, just removing the "|| c == ' '" might actually be enough in this context because 
  // before we even call "isInAssignment", we check, if the next character to the right is a '=' 
  // and if it isn't, we are already done and don't even get here. I think, within the current 
  // feature set, the only way that freeform strings can possibly occur in an sfz is in comments 
  // and in filenames where in filenames, the '=' is explicitly forbidden, so, I think, we shold be
  // on the safe side. If, however, we want to extend the feature set with more opcodes that allow 
  // freeform strings (such as a "formula" opcode within which we definitely want to allow '=', we
  // will have to rethink/revise this code here. Maybe we will indeed need some way to create 
  // tree-like document-structure, i.e. a class like SfzDocument that contains nodes which can be 
  // group/region/opcode/comment/etc. and maybe a proper class SfzParser that creates such a 
  // structured document from a string. Maybe the class should contain the orginal string together
  // with metadata about the structure. The "Nodes" could contain a start (and maybe end) position
  // within the string, etc. ...we'll see....
  // DAMN: leaving out the test aginst ' ' makes it fail when loading the BassdrumPsy1.sfz patch
  // in ToolChain, selecting the volume opcode of the 1st region and tweaking it (we hit an 
  // assert). The unit test with " pan2=0 pan1=0" even hangs. So for the time being, I leave the 
  // check in which may imply that locating the correct code segment may fail when the patch 
  // contains space in sfz filenames, so for the time being, we disallow spaces in filenames. This
  // needs to be fixed some other day....We should file this in a "known bugs" list.

  // Returns true, if the string starting at the given position is not the suffix of some longer 
  // string:
  auto isNoSuffix = [&](const std::string& code, size_t pos)
  {
    if(pos == 0) 
      return true;
    else
      return code[pos-1] == ' ' || code[pos-1] == '\n' || code[pos-1] == '\t' 
          || code[pos-1] == '>';
    // Any substring of the code that doesn't have a whitespace, newline, tabulator or '>' 
    // immediately to the left, is considered to be a suffix of a longer string. The case for which
    // no character to the left exists is treated separately. With this test, we avoid false 
    // positives in cases  such as when "tune" appears as suffix of some "...detune" parameter or 
    // cutoff appears as suffix of lfo1_cutoff. The '>' character is allowed as left neighbor 
    // because opcodes can appear immediately after the closing angle bracket of a <group> or 
    // <region> tag.
  };

  // Helper function to determine whether a found susbtring s that *looks like* an instance of the 
  // desired opcode definition really *is* one. This function is used to weed out the false 
  // positives that may occur due to the following conditions:
  // -s may occur as part of another opcode string. For example, the string "cutoff" appears in 
  //  "cutoff_ccN". These are weeded out by checking, that the next character is a '='.
  // -s may occur as part of a comment -> weeded out by isInComment()
  // -s may occur as part of a filename or other kind of free-form string that occurs on a 
  //  right-hand-side of an assignment -> weeded out by isInAssignment()
  auto meetsCriteria = [&](const std::string& code, size_t startPos, size_t endPos)
  {
    if(startPos == string::npos) return false; // This check may be redundant
    if(endPos > searchEnd-1)     return false; // Must be in search range
    if(code[endPos+1] != '=')    return false; // Must be followed by '='
    return !isInComment(code, startPos) && !isInAssignment(code, startPos)
      && isNoSuffix(code, startPos);
    // Are these constraints really enough to catch all false positives or do we need to impose 
    // further constraints? -> Implement more unit test! We have some but the coverage of all the
    // possible scenarios by our tests is still far from being exhaustive....
  };
  // Maybe move this function also out of SfzCodeBook::findOpcode and give it a more proper name
  // like "isOpcodeDefinition" or "isActualOpcodeDefinition"



  while(true) // The stopping conditions are complex and handled by breaks
  {
    size_t pos1 = code.rfind(ptn1, end-L1+1);  // or should it be end-L1?

    size_t pos2 = string::npos;
    if(opcodeIndex == 1)                    // Special case for optional index 1
      pos2 = code.rfind(ptn2, end-L2+1);    // or should it be end-L2?

    // If neither ptn1 nor ptn2 was found within the current search range, we return without having
    // found anything:
    if(  (pos1 == string::npos || (int) pos1 < searchStart)
      && (pos2 == string::npos || (int) pos2 < searchStart) )
      break;

    // If only one of the two found positions meets our additional constraints/criteria, then this
    // is our result. If both meet the criteria, the one that appears later is our result. If none 
    // of them meets the criteria, we set "end" to the lower of the two positions and keep 
    // searching:
    bool pos1ok = meetsCriteria(code, pos1, pos1 + L1 - 1);
    bool pos2ok = meetsCriteria(code, pos2, pos2 + L2 - 1);
    if(pos1ok && !pos2ok)
    {
      *startIndex = (int)  pos1;
      *endIndex   = (int) (pos1 + L1) - 1;
      break;
    }
    else if(pos2ok && !pos1ok)
    {
      *startIndex = (int)  pos2;
      *endIndex   = (int) (pos2 + L2) - 1;
      break;
    }
    else if(pos1ok && pos2ok)
    {
      if(pos1 > pos2)
      {
        *startIndex = (int)  pos1;
        *endIndex   = (int) (pos1 + L1) - 1;
        break;
      }
      else if(pos2 > pos1)
      {
        *startIndex = (int)  pos2;
        *endIndex   = (int) (pos2 + L2) - 1;
        break;
      }
      else
      {
        // pos1 == pos2. This happens when the opcode without 1 is a prefix of the one with 1. In 
        // this case, both variants were found and the one *with* the 1 is the one we are 
        // interested in because the one without is just its substring ...i think.
        *startIndex = (int)  pos1;
        *endIndex   = (int) (pos1 + L1) - 1;
        break;
        // This branch can actually be absorbed into 1st branch by using >= instead of >
      }

    }
    else if(!pos1ok && !pos2ok)
    {
      end = std::min(pos1, pos2);
      if(end <= searchStart)
        break;
    }
  }

  // Notes:
  // -The important difference to findGroup/Region is how we interpret the index. In the methods
  //  above, the index was just a count of how many times the <group> or <region> opcode already 
  //  had appeared within the search region. Here, the opcodeIndex is a part of the search-pattern
  //  (and it is optional when it's equal to 1) and we always try to find the *last* occurence of a
  //  matching string as opposed to the index-th occurence as we did in the other methods.
}

void SfzCodeBook::findOpcodeValueString(const std::string& code, int groupIndex, int regionIndex, 
  Opcode op, int opIndex, int* startPos, int* endPos)
{
  *startPos = -1;
  *endPos   = -1;

  // Find locations in the code, where the group definition starts and ends:
  int groupStart, groupEnd;
  findGroup(code, groupIndex, &groupStart, &groupEnd);
  if(groupStart == -1)
    return;
  //RAPT::rsAssert(groupStart != -1, "Group not found in code");

  // Find locations in the code, where the region definition starts and ends:
  int regionStart, regionEnd;
  findRegion(code, regionIndex, groupStart, groupEnd, &regionStart, &regionEnd);
  if(regionStart == -1)
    return;
  //RAPT::rsAssert(regionStart != -1, "Region not found in code");

  // ToDo - just call:
  //findOpcodeValueString(code, op, opIndex, regionStart, regionEnd, startPos, endPos);
  // here and return - avoid a lot of duplication


  // Find locations in the code, where the opcode definition starts and ends (this does not include 
  // value):
  int opcodeStart, opcodeEnd;
  findOpcode(code, op, opIndex, regionStart, regionEnd, &opcodeStart, &opcodeEnd);
  if(opcodeStart == -1)
    return;
  //RAPT::rsAssert(opcodeStart != -1, "Opcode not found in code");

  // Outdated:
  // This fails for some patches because we use the forward slash as seperator in the file-names, 
  // which is wrong anyway (rgc:sfz doesn't accept it, for example). ToDo: change the 
  // forward-slashes to backslashes in the patches and see, if this fixes the problem. OK - yes
  // it does (done only for the 1st patch - ToDo: fix this for all patches). However, the returned 
  // range currently only inludes the string for the opcode name itself. The part of the code that 
  // we need to modify is actually the value that comes after it. For example for something like 
  // volume=-6.02, opcodeStart returns the position of the 'v' and opcodeEnd the position opf the 
  // 'e'. We now need to figure out the range of the "-6.02" string. ...that should be easy, 
  // though. Maybe all of that can be wrapped into a convenience function 

  // To figure out the value's start/end position, scan rightward until we encounter a 
  // ' ', '\n', '\t', '/'
  int codeLength = (int) code.length();
  RAPT::rsAssert(opcodeEnd <= regionEnd);
  RAPT::rsAssert(opcodeEnd <= codeLength-3); // ex.: "pan=2" has length 5, the 'n' is at 2
  int i = opcodeEnd+1;
  RAPT::rsAssert(code[i] == '=');
  while(i <= regionEnd)
  {
    char c = code[i];
    if(c == ' ' || c == '\n' || c == '\t' || c == '/')
      break;
    ++i;
  }
  *startPos = opcodeEnd+2;
  *endPos   = i-1;
  // factor this out! It's used also in the func below


  // ToDo:
  // -Maybe remove some of the assertions and just return -1,-1 in such a cases. I think, we should 
  //  be able to handle this error condition gracefully.
  int dummy = 0;
}

void SfzCodeBook::findOpcodeValueString(const std::string& code, Opcode op, int opIndex,
  int searchStart, int searchEnd, int* startPos, int* endPos)
{
  *startPos = -1;
  *endPos   = -1;
  int opcodeStart, opcodeEnd;
  findOpcode(code, op, opIndex, searchStart, searchEnd, &opcodeStart, &opcodeEnd);
  if(opcodeStart == -1)
    return;
  //RAPT::rsAssert(opcodeStart != -1, "Opcode not found in code");

  // Factor out! Has been copied and pasted from the function above - but we have added the '\r'
  // here because in ToolChain in the debugger, I've seen such an '\r'
  // To figure out the value's start/end position, scan rightward until we encounter a 
  // ' ', '\n', '\t', '/'
  int codeLength = (int) code.length();
  RAPT::rsAssert(opcodeEnd <= searchEnd);
  RAPT::rsAssert(opcodeEnd <= codeLength-3); // ex.: "pan=2" has length 5, the 'n' is at 2
  int i = opcodeEnd+1;
  RAPT::rsAssert(code[i] == '=');
  while(i <= searchEnd)
  {
    char c = code[i];
    if(c == ' ' || c == '\n' || c == '\t' || c == '/' || c == '\r')
      break;
    ++i;
  }
  *startPos = opcodeEnd+2;
  *endPos   = i-1;


  // ToDo:
  // -Maybe remove some of the assertions and just return -1,-1 in such a cases. I think, we should 
  //  be able to handle this error condition gracefully.
  int dummy = 0;

}

// Get rid of this function:
void SfzCodeBook::findOpcodeAssignment(const std::string& code, Opcode opcode, int opcodeIndex,
  int searchStart, int searchEnd, int* startIndex, int* endIndex)
{
  // Old;
  findOpcode(code, opcode, opcodeIndex, searchStart, searchEnd, startIndex, endIndex);
  // ToDo:
  // -Adjust endIndex to point to the last digit of the numeric value in case of numeric parameters
  //  or the last letter of a choice or string parameter. On return of findOpcode, it should point
  //  to the character immediately before the '='.

  /*
  // New 2023/05/22 - does not yet work:
  int opcodeStart, opcodeEnd;
  findOpcode(code, opcode, opcodeIndex, searchStart, searchEnd, &opcodeStart, &opcodeEnd);
  RAPT::rsAssert(opcodeStart != -1, "Opcode not found in code");
  int i = opcodeEnd+1;
  if(i > searchEnd) {
    *startIndex = *endIndex = -1; }
  RAPT::rsAssert(code[i] == '=');
  while(i <= searchEnd)
  {
    char c = code[i];
    if(c == ' ' || c == '\n' || c == '\t' || c == '/')
      break;
    ++i;
  }
  //*startIndex = opcodeEnd+2;
  *endIndex   = i-1;
  */

  int dummy = 0;
}


}}

/*

ToDo:
-replace all calls to to_string with rsFloatToString

-Turn it into a Singleton: 
 -mostly done, still to do: make constructor and assignment operator protected
-Add unit and spec fields to the opcode records. Maybe it should be possible to configure 
 the sampler such that ignores certain kinds of opcodes, e.g. recognizes only sfz1 opcodes. That 
 way, we could audit how an sfz patch would sound on a sampler that doesn't support one of the
 extended specifications.
-Write a unit test that loops through all opcodes and translates them to a string and back. It 
 should also ensure that the op member of the record is equal to the array index. Wrtie similar
 tests for back-and-forth conversion of the string-type parameters.
-Maybe write a little program that generates a string which lists all opcodes with their 
 properties, perhaps sorted by name, dsp-type, etc. This could be used to generate a little 
 textfile as reference manual which would be a convenient thing to have. Maybe it could even be 
 member function here: generateReferenceManual() or something. It could perhaps take a couple of
 parameters to select formatting and sorting options. Maybe it could also be possible to show
 all opcodes that belong to a particular DSP or a particular spec, ie. filtering options.
-To support the opcodes with an attached "N", i.e. those to which the author can append a number,
 we should strip off the number in stringToOpcode and maybe return a std::pair<Opcode, int> 
 instead of just an opcode. In the toString conversion, maybe have an optional index parameter
 defaulting to -1 and append it to the string when != -1. Test this by trying two filters in 
 series (maybe highpass and lowpass). Maybe test the following dsp chain: 
   lowpass -> waveshaper -> lowpass -> waveshaper -> equalizer
 with 3 filters and 2 waveshapers

To add support for a new opcode, the following things need to be done:
-If it's not already there, a new entry must be added to the Opcode enum
-In our constructor, we'll need to add a corresponding add(...) call
-PlayStatus may need to get a new field for it
-RegionPlayer::setupPlayerSetting and SampleBusPlayer::setupPlayerSetting may need to add a new
 case in their switch statements
-If the opcode affects one of our existing DSPs, we may need to add a parameter to this DSPs list,
 check rsSamplerProcessors::Amplifier::Amplifier() for this. If the new opcode is for some new type 
 of DSP which does not yet exist, a new DSP class may need to be introduced - see below
-If the opcode affects the RegionPlayer, we must account for it in RegionPlayer::setupDspSettingsFor
 (or rather one of the function it calls...-> explain better)
-quite probably, prepareToPlay of the DSP class also needs to be updated
-...what about controller handling?
-If the opcode's values are of textual type, new enums for the possible choices have to be defined 
 and appropriate conversion functions need to be added to SfzCodeBook, see enum LoopMode, 
 SfzCodeBook::stringToLoopMode, SfzCodeBook::loopModeToString
-Probably a unit test should be written and added to the suite

To add a new DSP class, the following things need to be done:
...tbc...


SFZ - Resources:
https://sfzformat.com/legacy/   opcode reference for sfz 1.0
https://www.linuxsampler.org/sfz/    has convenient list of opcodes, also for sfz v2
https://en.wikipedia.org/wiki/SFZ_(file_format)
https://github.com/sfz/tests/   test sfz files demonstrating various features
https://sfzformat.com/headers/  reference for organization levels (group/region/...) sfz files
http://www.drealm.info/sfz/plj-sfz.xhtml  description of the sfz format
https://www.kvraudio.com/forum/viewtopic.php?f=42&t=508861  kvr forum thread with documentation
https://sfzinstruments.github.io/  collection of sfz instruments
http://ariaengine.com/overview/sfz-format/
http://doc.linuxsampler.org/sfz/
https://noisesculpture.com/cakewalk-synthesizers/
https://noisesculpture.com/cakewalk-synthesizers-downloads/

SFZ Players:
https://sfzformat.com/software/players/  players (also open source)

https://plugins4free.com/plugin/217/   sfz by rgcaudio, the original legacy vst2 plugin (+ standalone)
installation path: C:\Program Files (x86)\Vstplugins

https://sfz.tools/sfizz/downloads   sfizz, open source, vst3 + lv2, supports many sfz2 opcodes
installation paths: C:\Program Files\Common Files\VST3\sfizz.vst3
                    C:\Program Files\Common Files\LV2\sfizz.lv2


open source sfz players:
https://github.com/swesterfeld/liquidsfz/
https://sfz.tools/sfizz/downloads
https://github.com/altalogix/SFZero/
https://github.com/s-oram/Grace/

sfz compatible samplers
https://github.com/christophhart/HISE/
https://github.com/surge-synthesizer/shortcircuit-xt (not sure if it supports sfz)

deeper into the codebases:
https://github.com/swesterfeld/liquidsfz/tree/master/lib
https://github.com/swesterfeld/liquidsfz/blob/master/lib/synth.hh
This seems to do it the simple way: it has a fixed number of voices and if they are used up, no
more can be added - if i understand it correctly (see alloc_voice, line 230)

SFZ instruments (i.e. sample packs):
https://www.pianobook.co.uk/  see: https://www.youtube.com/watch?v=os6fg2t8uK4   ...website seems broken

*/