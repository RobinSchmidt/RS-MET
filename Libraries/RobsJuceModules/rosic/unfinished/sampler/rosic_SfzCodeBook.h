#ifndef rosic_SfzCodeBook_h
#define rosic_SfzCodeBook_h
namespace rosic { namespace Sampler {

//-------------------------------------------------------------------------------------------------
/** Enumeration of possible types of settings. These types correspond to the opcodes defined
in the sfz specification. At the moment, only a small fraction of the opcodes are supported by the
engine but I have already integrated the full list into the enum.
References:
  https://sfzformat.com/legacy/     */

enum class Opcode
{
  Unknown = 0, Unsupported,

  // Sample Definition:
  Sample,  // includes file extension (e.g. Piano.wav). path is realtive to the .sfz file.

  // Response Constraints (aka Input Controls):
  LoChan, HiChan, LoKey, HiKey, Key, LoVel, HiVel, LoCtrlN, HiCtrlN, LoBend, HiBend, 
  LoChanAft, HiChanAft, LoPolyAft, HiPolyAft, LoRand, HiRand, LoBpm, HiBpm,
  SeqLength, SeqPosition, SwLoKey, SwHiKey, SwLast, SwDown, SwUp, SwPrevious, SwVel,
  Trigger, Group, OffBy, OffMode, OnLoCtrlN, OnHiCtrlN,

  // Sample Player:
  Delay, DelayRandom, DelayCtrlN, Offset, OffsetRandom, OffsetCtrlN,
  End, Count, LoopMode, LoopStart, LoopEnd, SyncBeats, SyncOffset,

  // Pitch:
  Transpose, Tune, PitchKeyCenter, PitchKeyTrack, PitchVelTrack, PitchRandom,
  BendUp, BendDown, BendStep,

  // Pitch Envelope:
  PitchEnvDelay, PitchEnvStart, PitchEnvAttack, PitchEnvHold, PitchEnvDecay, PitchEnvSustain,
  PitchEnvRelease, PitchEnvDepth, PitchEnvVel2Delay, PitchEnvVel2Attack, PitchEnvVel2Hold, 
  PitchEnvVel2Decay, PitchEnvVel2Sustain, PitchEnvVel2Release, PitchEnvVel2Depth,

  // Pitch LFO:
  PitchLfoDelay, PitchLfoFade, PitchLfoFreq, PitchLfoDepth, PitchLfoDepthCtrlN,  
  PitchLfoDepthChanAft, PitchLfoDepthPolyAft, PitchLfoFreqCtrlN, PitchLfoFreqChanAft, 
  PitchLfoFreqPolyAft,

  // Filter:
  /*FilType, Cutoff,*/ CutoffCtrlN, CutoffChanAft, CutoffPolyAft, /*Resonance,*/
  /*FilKeyTrack, FilKeyCenter, FilVelTrack,*/ FilRandom,

  // new filter opcodes, adorned with index:
  filN_type, cutoffN, resonanceN,
  filN_keytrack, filN_keycenter, filN_veltrack, 

  // todo: complete the list, then delete old, unindexed opcodes

  // Filter Envelope:
  fileg_delay, fileg_start, fileg_attack, fileg_peak, fileg_hold, fileg_decay, fileg_sustain, 
  fileg_release, fileg_end, fileg_attack_shape, fileg_decay_shape, fileg_release_shape,
  fileg_depth,
  // Should we use filegN? SFZ2 does actually not define fileg2 opcodes! It defines stuff like 
  // cutoff2, resonance2, fil2_type, fil2_gain. What is the supposed behavior? Will the fileg_
  // opcode affect both filters or only filter 1? If it affects both - should we generalize it to
  // affect all filters? ...that's what we currently do

  // Need to be converted to the new style:
  // FilEnvVel2Delay, FilEnvVel2Attack, FilEnvVel2Hold, 
  // FilEnvVel2Decay, FilEnvVel2Sustain, FilEnvVel2Release, FilEnvVel2Depth,

  // Filter LFO:
  FilLfoDelay, FilLfoFade, FilLfoFreq, FilLfoDepth, FilLfoDepthCtrlN, 
  FilLfoDepthChanAft, FilLfoDepthPolyAft, FilLfoFreqCtrlN, FilLfoFreqChanAft, FilLfoFreqPolyAft,
  fillfo_freq, fillfo_depth,


  // Amplifier:
  /*Volume, Pan, Width, Position, */  // replaced by volumeN, etc.
  /*AmpKeyTrack, AmpKeyCenter, AmpVelTrack,*/ AmpVelCurveN, AmpRandom,
  RelTrigDecay, Output, GainCtrlN,
  FadeInLoKey, FadeInHiKey, FadeOutLoKey, FadeOutHiKey, FadeCurveKey,
  FadeInLoVel, FadeInHiVel, FadeOutLoVel, FadeOutHiVel, FadeCurveVel,
  FadeInLoCtrlN, FadeInHiCtrlN, FadeOutLoCtrlN, FadeOutHiCtrlN, FadeCurveCtrl,

  // new amp-opcodes - will make the above ones obsolete:
  amplitudeN,
  volumeN, panN, widthN, positionN, 
  ampN_keytrack, ampN_keycenter, ampN_veltrack,
  // ampN_velcurve_X /* ? */, ampN_random,

  // Amplifier Envelope:
  ampeg_delay, ampeg_start, ampeg_attack, ampeg_peak, ampeg_hold, ampeg_decay, ampeg_sustain, 
  ampeg_release, ampeg_end, ampeg_attack_shape, ampeg_decay_shape, ampeg_release_shape,
  ampeg_depth,
  // We are not using ampegN_delay etc. because there is only one hardwired ampeg which is always
  // wired to the last Amplifier in the effect chain.

  // Need to be converted to the new style:
  // /*AmpEnvDepth,*/ AmpEnvVel2Delay, AmpEnvVel2Attack, AmpEnvVel2Hold, 
  //AmpEnvVel2Decay, AmpEnvVel2Sustain, AmpEnvVel2Release, /*AmpEnvVel2Depth,*/
  //AmpEnvDelayCtrlN, AmpEnvStartCtrlN, AmpEnvAttackCtrlN, AmpEnvHoldCtrlN, AmpEnvDecayCtrlN, 
  //AmpEnvSustainCtrlN, AmpEnvReleaseCtrlN, 
  // depth-parameters do not exist in sfz for amp-env. wouldn't make much sense, i guess

    // Amplifier LFO:
  AmpLfoDelay, AmpLfoFade, /*AmpLfoFreq, AmpLfoDepth,*/ AmpLfoDepthCtrlN, 
  AmpLfoDepthChanAft, AmpLfoDepthPolyAft, AmpLfoFreqCtrlN, AmpLfoFreqChanAft, AmpLfoFreqPolyAft,
  amplfo_freq, amplfo_depth,

  //amplfo_amp, 
  // For consistency with lfoN_amp opcode. It doesn't really serve a meaningful purpose because we
  // already have amplfo_depth. It's redundant and not in the original sfz spec.


  // Equalízer:
  eqN_freq, eqN_freqccX, eqN_vel2freq, eqN_bw, eqN_bwccX, eqN_gain, eqN_gainccX, eqN_vel2gain,

  // Effects:
  effN,         // effect send levels in percent (eff1: reverb, eff2: chorus)

  // SFZ 2.0:
  // fil2_type, cutoff2, etc. ...generalize to filN_type, cutoffN, etc. - maybe for translations
  // and lookup, we'll need a special rule for filter-related parameters that leaves out the 1,
  // i.e. fil1_type, cutoff1, etc. should not be produced or accepted - instead, we'll use 
  // fil_type, cutoff, etc. because that's what sfz 1 has. it doesn't have a concept of multiple
  // filters

  // ADSR+ envelope generators (opcodes not defined in SFZ except sustain):
  //egN_start, egN_delay, egN_attack, egN_peak, egN_hold, egN_decay, egN_sustain, egN_release, 
  //egN_end,

  adsrN_start, adsrN_delay, adsrN_attack, adsrN_peak, adsrN_hold, adsrN_decay, adsrN_sustain, 
  adsrN_release, adsrN_end, adsrN_attack_shape, adsrN_decay_shape, adsrN_release_shape,



  // egN_timeX, egN_levelX
  lfoN_freq,      // lfoN_delay, lfoN_fade, lfoN_phase, lfoN_wave,
  lfoN_volumeX,
  lfoN_amplitudeX, 
  // is this supposed to control the amplitude of the LFO..or the routing of the LFO output to an
  // amplitude parameter? in the former case, we need nee no second index X, in the latter, we do



  // ARIA:
  PanLaw,


  //  RS-MET
  // My own extensions - preliminary, for experimentation. Before defining extensions, check what 
  // is already there in SFZ 2.0 and in other sfz engines with extensions. Try to be compatible 
  // with the  largest possible range of other engines. Discuss extensions on KVR before 
  // implementing them in production code:
    
  distortN_shape, distortN_drive, distortN_dc, 
  //distortN_shape, distortN_drive, distortN_dc, // DistGain...may be redundant with Volume

    
  lfoN_amp, 
  // ToDo: Check, if there is an opcode in sfz2 that controls the amplitude of the LFO. Maybe
  // lfoN_amplitude is made for that? But I suppose, that controls the routing of lfoN to the
  // signal amplitude


  // SFZ2 has opcode egN_driveshape but no driveshape as such? same with lfoN_drive. does that 
  // mean the LFO signal is drive into a saturator?
  // check: https://www.plogue.com/products/sforzando.html

  // eff1_type, eff2_type (can be reverb, chorus, echo, convolution)
  // reverb_time_scale, reverb_density, reverb_size, reverb_time_lo, reverb_time_mid, 
  // reverb_time_high, reverb_freq_lo, reverb_freq_hi, reverb_algo (fdn16, schroeder, ...)
  // chorus_voices, chorus_freq, chorus_depth, ...
  // convolve_sample (or maybe just convo_sample)
  // sample_dir: default: directory of sfz-file, may contain ../../Samples to move out of the 
  // sfz-dir...maybe using ../../Samples directly in the sample-path even works already wihtout
  // doing anything?
  // oversample=1..16 per region oversampling factor, always has override behavior
  // extend_ranges=true/false allow extended parameter ranges such a width of 150, volume of +12 dB
  // etc. or maybe use a param_range=sfz/extended that would allow to later define more specific
  // extensions
  // add: loop_start_ccN, loop_end_ccN

  // What about routing?
  // is an opdoce like lfoN_volume meant to route lfo N to volume? and what isegN_amplitude 
  // supposed to do? Will the envelope with index N just get added to the already existing amp 
  // envelope?


  // Muted: convenient to switch regions or groups off wihthout removing them - but maybe, we can
  // just set the volume to -inf ...should be support -inf in the parser? maybe...

  NumTypes  // rename to NumOpcodes
};
// The underlying integers for the opcodes in this enum must start at zero and increment by one.
// Dont do something like FilterCutoff = 1000 etc. They must be usable as array indices and we dont
// want to waste space
// 
// maybe don't capitalize first letter - make it conistent with other (newer) enums in the 
// library. see community stadards recommendations...or maybe use names equal to the sfz opcode
// names. i think, that would be more convenient

//-------------------------------------------------------------------------------------------------

/** Enumeration of the waveforms that are available in the LFOs. */
enum class WaveForm
{
  unknown = 0,

  sine,
  triangle,
  square,
  saw_up,
  saw_down,

  sample,          // let's the user specify a cingle-cycle .wav file

  numWaveForms
};
// ToDo: compare names with sfz-spec - maybe the names are different there - if so, rename


/** Enumeration of the different filter types that are available. */
enum class FilterType // maybe rename to fil_type for consistency with sfz
{
  Unknown = 0,
  off,           // maybe remove - it's not defined in sfz

  lp_6, lp_12, hp_6, hp_12, bp_6_6, br_6_6, // SFZ 1.0 types
  // use lpf_1p, hpf_1p, lpf_2p, hpf_2p, bpf_2p, brf_2p  as in sfz spec


  // My own additional types: ls stands for "low-shelf", hs for "high-shelf", pk for "peak" aka
  // bell or sometimes band-shelf. The equalizer opcodes of sfz are realized using the same DSP 
  // objects as for the filter opcode, just with the pk_2p filter type.
  ls_1p, ls_2p, hs_1p, hs_2p, pk_2p, 
  // Hmm... ARIA has lsh, hsh, peq - maybe use these. Are the shelvers there 1-poles or 2-poles?


  numTypes
};
// see: https://sfzformat.com/opcodes/fil_type

enum class LoopMode   // maybe rename to loop_mode (as in sfz)
{
  Unknown = 0,

  // SFZ1:
  no_loop,         // play from start to end or until note off
  one_shot,        // play from start to end, ignore note off, engaged if count opcode is defined
  loop_continuous, // when player reaches sample loop point, loop will play until note expiration
  loop_sustain,    // play loop while note or sustain pedal is held. rest will play after release

  // RS-MET:
  //single_cycle,
  // Interacts with modulation of loop_start/loop_end to preserve pitch when these are modulated. 
  // Maybe we should also be able to say, how many cycles the loop represents, i.e. support 
  // multi_cycle. Maybe in this mode, pitch_keycenter should be ignored


  numModes
};

enum class PanRule  // maybe rename to pan_law (as in aria engine)
{
  linear, sinCos,

  numPanRules

  // Aria engine defines 2 laws: balance and mma:
  //   https://sfzformat.com/opcodes/pan_law
  // but here mma is capitalized as MMA:
  //   http://ariaengine.com/forums/index.php?p=/discussion/4389/arias-custom-opcodes/p1
  // so we should perhaps be case insensitive in the parser...maybe quite generally so? But what do
  // these laws actually mean? is balance the sinCos rule and mma the linear? Or the other way 
  // around? ...measure it! maybe we could provide a continuously adjustabel pan law which is given
  // as float in dB center attenuation? maybe the opcdoe whould be pan_center_level
  // https://www.soundonsound.com/sound-advice/q-what-pan-law-setting-should-use
};
// Or maybe it should be called PanLaw? Maybe have different variations with respect to total
// gain - for linear: either factor 2 for hard left/right setting or a factor of 0.5 for a 
// center setting. The former would imply that with neutral default settings, stereo samples 
// are played as is. The latter would imply that hard-panned mono samples would be played as is
// on their respective channel. Both behaviors may be useful, although, it would be a bit 
// redundant because we also have an overall gain parameter as well which can always be used to 
// compensate...although a factor of exactly 2 or 0.5 may be hard to achieve because gain is 
// given in dB, so the sfz file would have to specify +-6.0205999132796239....., which is 
// inconvenient


enum class DistortShape
{
  Unknown = 0,

  linear,           // no distortion at all
  clip,             // hard-clipper at +-1
  tanh, 
  soft_fold,        // x / (1 + x^2)

  // todo: atan, soft_clip, sin, asinh, sinh, sinh(a * asinh(x)/a), asinh(a * sinh(x)/a),

  numDistortShapes
};



//-------------------------------------------------------------------------------------------------
/** Enumeration of the different data formats of the values of the opcode. */

enum class OpcodeFormat
{
  Unknown = 0,

  Boolean,
  Natural,  // a.k.a. unsigned int
  Integer,
  Float,
  String,    // rename to Text and use it only for freeform strings, add a "choice" format for
             // strings taken from a specified set (such as lpf_1p, hpf_1p, etc.)

  NumTypes
};

//-------------------------------------------------------------------------------------------------
/** Enumeration of the physical (or mathematical) units that apply to the value of the opcode. */

enum class OpcodeUnit
{
  Unknown = 0,
  Text, Index, RawInt, RawFloat,                                 // Raw
  Hertz, Octaves, Semitones, Cents, MidiKey,                     // Frequency
  Seconds, Samples, //Beats, //Milliseconds, //BeatsOrSeconds,   // Time
  Percent, Decibels, DecibelPerKey, CentPerKey, //Degrees,       // Misc
  NumUnits

  // The MidiKey can be given as midi note number, e.g. 61 or as string e.g. c#4. We need to allow
  // both syntaxes in the parser (and maybe also C#4 with capital C).
};

//-------------------------------------------------------------------------------------------------
/** Enumeration of the different specifications in which the respective opcode was defined. */

enum class OpcodeSpec
{
  Unknown,

  Sfz_1,
  Sfz_1_E,      // extended sfz 1, for example with arbitrary number of eq bands
  Sfz_2,
  Sfz_2_E,      // e.g.: lfoN_volumeX (sfz2 has only lfoN_volume)
  Aria,
  Aria_E,
  RsMet,
  //...etc...

  NumSpecs
};

//-------------------------------------------------------------------------------------------------
/** Enumeration of the different signal processor types that may be used in the definition of
instruments. What kinds of processors are used within a region is implicitly determined by the sfz 
opcodes, e.g. the presence of a FilterCutoff opcode dictates the presence of a filter within the 
respective region. In order to facilitating to build the DSP chain for a region player,
we also need an explicit representation of the DSP processor types. */

enum class OpcodeType   // Maybe rename to OpcodeTarget
{
  Unknown,

  // Sound production:
  SamplePlayer,

  // Effects:
  _TagEffectsStart,
  Amplifier,
  Filter,
  Equalizer,
  WaveShaper,
  _TagEffectsEnd,

  // Fixed modulators:
  _TagFixedModulatorsStart,
  AmpEnv,    AmpLfo,
  FilterEnv, FilterLfo,
  PitchEnv,  PitchLfo,
  _TagFixedModulatorsEnd,

  // Freely routable modulators:
  _TagFreeModulatorsStart,
  FreeEnv,    // rename to EnvGen or FreeAdsr
  FreeLfo,    // rename to LowFreqOsc

  // Allow midi-inputs to be modulation sources, too:
  //MidiCtrl, MidiKey, MidiOnVel, MidiOffVel, MidiAftertouch, ...

  _TagFreeModulatorsEnd,
  // Free modulators need to come immediately after the fixed modulators. Some code relies on that.
  // ToDo: maybe avoid this dependency by using _TagModulatorsStart, _TagModulatorsEnd and wrap 
  // them all between these tags ...but maybe some code needs to distinguish the cases...we'll see
  // wehn the mod-system is finished...we coul actually also have both by nesting those tags

  // Routing of modulators:
  _TagModRoutingStart,
  EnvN_ParamX,
  LfoN_ParamX,
  HardwiredModRouting,
  _TagModRoutingEnd,


  // Opcodes for controlling key- and velocity tracking:
  //Tracking,           // no - we don't need this - the tracking parameters bleong into the same
  //_TagTackingStart, // category as their target parameters
  //_TagTrackingEnd,



  NumDsps
};

//=================================================================================================

/** A class to represent various playback settings of a region, group or instrument. Such settings 
include constraints for the circumstances under which a particular sample should be played. Key- 
and velocity ranges are the obvious primary constraints (and they therefore are directly baked into
the Region class below), but sfz defines many more. But settings don't need to be playback 
constraints. Other types are things like envelope settings, filter frequencies, etc. */

class PlaybackSetting
{

public:

  PlaybackSetting(Opcode type = Opcode::Unknown, float value = 0.f, int index = -1) 
  { 
    this->type  = type; 
    this->value = value; 
    this->index = index;
  }
  // maybe remove the -1 default argument. callers should be explicit...maybe...or: check with an
  // assert that no index applies to the given Opcode when index == -1

  Opcode getOpcode() const { return type; }


  /** Returns the stored value for this setting. Values are always stored as floats and it is
  understood that in cases, where the corresponding parameter in the sfz spec is defined to be an
  integer, we just represent it as that integer with all zeros after the decimal dot and the
  caller is supposed to convert to int by writing e.g.:

  int intValue = (int)setting.getValue();

  For settings whose value is represented by text that indciates a particular choice, the caller
  has to look up the int in an enumeration corresponding to the type of the parameter...tbc... */
  float getValue() const { return value; }
  // todo: (decide and) document, how choice parameters like filter-type are represented. In sfz,
  // they are just represented as text. Maybe for each such choice parameter, we need another 
  // enum to represent its allowed values.

  /** Some settings need to specify an index in addtion to the value. An example is a setting
  involving midi control changes. In the sfz file, they are written as e.g. loccN where the N is
  replaced by the actual controller number, like locc74=20 to indicate, that the sample should
  only play, if the last received value for CC#74 was >= 20. If indexing not applicable to a
  particular setting/opcode, this will return -1. */
  int getIndex() const { return index; }

  /** Returns the default value for a given type of setting. */
  static float getDefaultValue(Opcode op, int index);

  /** For a given type of opcode, this function returns the type of signal processor to which the
  opcode applies, e.g. SignalProcessorType::Filter for Type::FilterCutoff. */
  static OpcodeType getTargetOpcodeType(Opcode op);
  // rename to getTargetOpcodeType


  OpcodeType getTargetOpcodeType() const { return getTargetOpcodeType(type); }


  bool operator==(const PlaybackSetting& rhs) const
  {
    return type == rhs.type && value == rhs.value && index == rhs.index;
  }


private:

  Opcode type  = Opcode::Unknown;  // rename type to opcode
  float  value = 0.f;
  int    index = -1;  //< Used e.g. for conrol-change settings. Is -1, if not applicable.
  // Maybe use 1 as default. If there's only one such setting anyway, that seems appropriate
  // index should always be a positive real number. But maybe that's not such a good idea - see
  // comment in SfzInstrument::writeSettingToString in the cpp file. For certain things, we need
  // a code for "not applicable".
  
  // Maybe store, where the opcode appears in the to facilitate highlighting the relevant 
  // section - but actually that would violate separation of concerns
};

//=================================================================================================

/** A class to represent the settings of a modulation routing. This class is used only to store the
settings, not to actually connect modulator objects to parameters during processing. For this 
purpose, the class ModulationConnection is used. ..tbc.... */

class ModulationRouting 
{

public:

  //---------------------------------------------------------------------------------------------
  // \name Lifetime

  ModulationRouting(OpcodeType modSrcType, int modSrcIndex, Opcode modTarget, int modTargetIndex,
    float modDepth, ModMode modMode);

  /** Standard constructor. Will create an invalid object. All enum members are initialized to 
  "unknown", indices to -1, etc. The object needs to be either set up correctly later via the 
  setters or can serve as a sort of null-object or dummy-object that can be used, when a 
  mod-routing string could not be parsed, for example. */
  ModulationRouting() {}


  //---------------------------------------------------------------------------------------------
  // \name Setup

  void setDepth(float newDepth) { depth = newDepth; }

  void setMode(ModMode newMode) { mode = newMode; }


  //---------------------------------------------------------------------------------------------
  // \name Inquiry


  /** Returns true, iff the type and index of the source of the given routing matches ours. */
  bool hasMatchingSource(const ModulationRouting& r) const
  {
    return r.sourceType == sourceType && r.sourceIndex == sourceIndex;
  }

  bool hasMatchingTarget(const ModulationRouting& r) const
  {
    return r.target == target && r.targetIndex == targetIndex;
  }

  bool hasMatchingEndpoints(const ModulationRouting& r) const
  {
    return hasMatchingSource(r) && hasMatchingTarget(r);
  }
  // try to find better name maybe hasMatchingPorts/Pins...dunno


  OpcodeType getSourceType() const { return sourceType; }

  int getSourceIndex() const { return sourceIndex; }


  OpcodeType getTargetType() const { return targetType; }

  int getTargetIndex() const { return targetIndex; }


  Opcode getTargetOpcode() const { return target; }
  // maybe rename.. the "Opcode" here means what sort of Parameter is being modulated

  ModMode getMode() const { return mode; }

  float getDepth() const { return depth; } 

  bool isInvalid() const
  {
    return sourceType == OpcodeType::Unknown || sourceIndex == -1 || mode   == ModMode::unknown
        || targetType == OpcodeType::Unknown || targetIndex == -1 || target == Opcode::Unknown;
  }

  bool operator==(const ModulationRouting& rhs) const
  {
    return sourceType == rhs.sourceType && sourceIndex == rhs.sourceIndex 
      && targetType == rhs.targetType && targetIndex == rhs.targetIndex && target == rhs.target
      && depth == rhs.depth && mode == rhs.mode;
  }
  // Why do we need to implement this ourselves? Shouldn't the compiler auto-generate it? All 
  // members are trivial data fields. However, we do get strange compiler errors in MSVC saying 
  // something about 
  //   "operator __surrogate_func': no matching overloaded function found"
  // if we don't implement it manually :-/


private:

  OpcodeType sourceType  = OpcodeType::Unknown; // Type of modulator, e.g. EnvGen, LowFreqOsc, etc.
  int        sourceIndex = -1;                  // Index like the 3 in lfo3_freq 
  OpcodeType targetType  = OpcodeType::Unknown; // Type of receiver object, e.g. Filter, Amplifier,
  int        targetIndex = -1;                  // Index like the 2 in lfo3_cutoff2, lfo3_eq2gain
  Opcode     target      = Opcode::Unknown;     // Target parameter like cutoffN, volumeN, etc.   // rename to targetPar
  float      depth       = 0.f;                 // Modulation depth
  ModMode    mode        = ModMode::unknown;    // absolute, relative, exponential, etc.

  // ToDo:
  // -Maybe have a member for modulationMode (absolute, relative, exponential, etc.) but maybe
  //  that should be determined by the target...not sure yet maybe one of the mode options could be
  //  "automatic" or "default"
  // -Maybe use sourceIndex values < 1 for the special fixed modulators, like
  //  0: amp, -1: filter, -2 pitch - maybe define constants for them. 
  //  Maybe -1: amp, -2: filter, -3: pitch is more convenient
};

//=================================================================================================

/** A class to translate between the sfz opcodes in their string representation and their enum 
values and some additional related "translations" such as the mapping of opcodes to the type of
signal processor to which they apply, etc. Essentially, it does many of the conversions and 
mappings that are required within the sampler engine for sfz-parsing or saving the state to an
sfz-file etc.

We need only one such object for the whole lifetime of the sampler engine and we need convenient 
access to it from various of its subsystems, so we implement a variation of the singleton pattern 
here.  */

class SfzCodeBook
{

public:

  //-----------------------------------------------------------------------------------------------
  // \name Lifetime. We implement a variation of the singleton pattern.

  /** Creates the sole instance of the class. In contrast to the textbook version of the singleton
  pattern where the object is implicitly created on first use and never deleted, we give explicit 
  control over the singleton's creation via this function which is supposed to be called by client
  code before using the singleton for the first time. We call this in the constructor of the 
  sampler engine when the first instance is created. */
  static void createInstance();

  /** Returns a pointer to the sole instance of the class which is assumed to have been previously
  created by createInstance. If you didn't explicitly create the instance, it will implicitly be 
  created within this call and in debug builds, you'll hit an assert yelling at you to be more 
  epxlicit. */
  static SfzCodeBook* getInstance();

  /** This function is supposed to be called by client code for cleaning up, when the instance is 
  not needed anymore. The sampler engine does this when its last instance is destroyed. */
  static void deleteInstance();

  //-----------------------------------------------------------------------------------------------
  // \name Translations

  /** Returns the type of the given opcode. */
  OpcodeType getOpcodeType(Opcode op);
  // maybe remove the "get" here or add a "get" to opcodeDefaultValue to make these functions
  // consistent..maybe rename this function to just getType and opcodeDefaultValue to 
  // getDefault(Value)

  OpcodeFormat getOpcodeFormat(Opcode op);


  // Maybe these housl go into the "Inquiry" section

  float opcodeMinValue(Opcode op);

  float opcodeMaxValue(Opcode op);


  /** Returns the default value for the given opcode as floating point number. If the format of the
  value is integer or an enum value, you'll need to convert the returned value to int and then
  possibly to the enum. Some opcodes contain an index like eq2_freq. For these, you need to pass 
  that index too because the default value may actually depend on that index. For example, the 
  defaults are 50,500,5000 respectively for eq1_freq, eq2_freq, eq3_freq. If indexing is not 
  applicable to the given opcode, you should pass index = -1. */
  float opcodeDefaultValue(Opcode op, int index);


  /** Returns the default modulation mode for the given opcode, i.e. the modulation mode that is 
  used when the user doesn't use a unit suffix to the modulation depth opcode in the sfz file. */
  ModMode opcodeDefaultModMode(Opcode op);

  /** Returns a string that represents the given opcode op with given index in a format that can be
  written into an .sfz file. The string is returned by value because for those opcodes that include 
  indices (such as eq2_freq), we need to generate the string dynamically. */
  std::string opcodeToString(Opcode op, int index) const;

  /** Translates an opcode string into the corresponding enum value. Some opcodes contain an index.
  For these, the index is passed in the "index" output parameter. For opcodes without index, it 
  will be assigned to -1. */
  Opcode stringToOpcode(const std::string& str, int* index);

  /** Translates the value of an opcode represented as floating point number into a corresponding
  string that can be written to an sfz file. */
  std::string valueToString(Opcode op, float value);

  /** Translates a string that represents the value for a given opcode into its floating point 
  representation. */
  float stringToValue(Opcode op, const std::string& str);

  /** Translates a filter type enum value to the string that represents it in an sfz file. */
  const std::string& filterTypeToString(FilterType ft);

  /** Translates a string representing a filter type into its corresponding enum value. */
  FilterType stringToFilterType(const std::string& str);


  LoopMode stringToLoopMode(const std::string& str);

  std::string loopModeToString(LoopMode loopMode);


  std::string modSourceToString(OpcodeType sourceType, int sourceIndex);

  std::string modTargetToString(OpcodeType targetType, int targetIndex, Opcode targetOpcode);

  std::string modDepthToString(float depth, ModMode mode, Opcode targetOpcode);

  OpcodeType stringToModSource(const std::string& str, int* index);

  float stringToModDepth(const std::string& str, ModMode* modMode, Opcode targetOpcode);

  ModMode stringToModMode(const std::string& str);


  std::string settingToString(const PlaybackSetting& s);

  std::string modRoutingToString(const ModulationRouting& r);




  // We need also a modRoutingToString(const ModulationRouting& routing) to show it in the tree
  // for consistency we also want settingToString(const PlaybackSetting& s)


  //-----------------------------------------------------------------------------------------------
  // \name Inquiry


  bool isValidOpcode(Opcode op)
  {
    return (int)op >= 0 && (int)op < (int)opcodeEntries.size();
  }

  /** Returns true iff the given opcode applies to the sample playback source such as tune, 
  delay, offset, etc. */
  bool isPlayerSetting(Opcode op)
  {
    OpcodeType type = getOpcodeType(op);
    return type == OpcodeType::SamplePlayer;
  }
  // todo: make static (maybe)


  /** Returns true iff the given opcode type applies to an effect in the effect chain like e.g.
  Amplifier, Filter, Equalizer, etc. */
  static bool isEffectSetting(OpcodeType type)
  {
    using OT = OpcodeType;
    return type > OT::_TagEffectsStart && type < OT::_TagEffectsEnd;
    // ToDo: maybe use a helper function:
    // return isStrictlyBetween(type, OT::_TagEffectsStart, OT::_TagEffectsEnd);
  }
  // maybe rename to isEffect


  /** Returns true iff the given opcode applies to an effect in the effect chain like e.g. volume, 
  cutoff, eq2_freq, etc. but not pitch_keycenter or tune. */
  bool isEffectSetting(Opcode op)
  {
    return isEffectSetting(getOpcodeType(op));

    /*
    using OT = OpcodeType;
    OT type = getOpcodeType(op);
    return type > OT::_TagEffectsStart && type < OT::_TagEffectsEnd;
    */
    // maybe use a helper function :
    // return isStrictlyBetween(type, OT::_TagEffectsStart, OT::_TagEffectsEnd);
  }
  // todo: make static

  static bool isModSourceSetting(OpcodeType type)
  {
    using OT = OpcodeType;
    return type > OT::_TagFixedModulatorsStart && type < OT::_TagFreeModulatorsEnd;
    // It is not a bug that in the > comparison we compeare against the "..Fixed.." tag and in 
    // the < comparision against the "..Free.." tag. We want to catch all ModSource settings here,
    // regardless whether its a fixed or free source.
  }

  bool isModSourceSetting(Opcode op) { return isModSourceSetting(getOpcodeType(op)); }
  // todo: needs test, make static


  static bool isModRoutingSetting(OpcodeType type)
  {
    using OT = OpcodeType;
    return type > OT::_TagModRoutingStart && type < OT::_TagModRoutingEnd;
  }

  bool isModRoutingSetting(Opcode op) { return isModRoutingSetting(getOpcodeType(op)); }

  // todo: isModulationSetting, isModulatorSetting  (modulatiON settings define mod-connections,
  // modulatOR settings parameters of the modulators), hasIndex..but this may be expensive to 
  // implement - we may go through the opcode string and try to find an 'N'. Not too bad for 
  // non-realtime code but perhaps not good to call in realtime. Or, maybe, if we need that at 
  // realtime add another field to the OpcodeEntry record. maybe a set of flags

  //-----------------------------------------------------------------------------------------------
  // \name Text search utilities

  static bool hasImplicitFirstGroup(const std::string& code);

  /** Finds the character index within the given sfzCode of the opening angle-bracket '<' of the 
  "<group>" header for the group with given index and assigns it to the output variable "startIndex".
  The endIndex is assigned to the last character before the '<' of the next group or the last 
  character in the string (if the group is the last one and no other group follows). */
  static void findGroup(const std::string& sfzCode, int groupIndex, int* startIndex, int *endIndex);
  // maybe rename startIndex/endIndex to startPos/endPos

  /** Similar to findSfzGroup but for the <region> header and with an additional specification of the
  search range within the code. The intention is that client code first figures out, where a given 
  group definition starts and ends and then searches for a region *within* the group definition, i.e.
  the "startIndex", "endIndex" outputs of findSfzGroup are suitable as inputs for "searchStart" and 
  "searchEnd" in a subsequent findSfzRegion call. */
  static void findRegion(const std::string& sfzCode, int regionIndex, int searchStart, int searchEnd,
    int* startIndex, int *endIndex);

  /** Similar to findSfzGroup and findSfzRegion but finds a single opcode definition. For example, if
  you want to find the substring cutoff2=765, you would pass Opcode::cutoff and 2 for the "opcode" 
  and "opcodeIndex" parameters and the function would assign "startIndex" to the character position 
  of the 'c' and "endIndex" to the character position of the '5'. If there are multiple definitions 
  of cutoff2 within the given search range, it will return the position of the last one. This is 
  because the last one is the one that will actually be used (re-assigning an opcode in sfz code 
  works like re-assigning variables in program code). */
  static void findOpcode(const std::string& sfzCode, Opcode opcode, int opcodeIndex, 
    int searchStart, int searchEnd, int* startIndex, int *endIndex);




protected:

  // Here, we may keep some "dictionary"-like data-structures to facilitate fast translations.
  // For the opcode-to-anything translations we want an O(1) complexity with small factors. It 
  // should essentially behave like a simple array access. The same goes for things like filter
  // types or basically all enum-types. These things may be called in the audio-thread when 
  // note-events are received. For string-to-anything conversions, we have less strict 
  // performance requirements because they happen only on preset load. However, we probably also 
  // want at most O(log(N)), i.e. maybe a lexicographically sorted array with binary search 
  // should be used or maybe std::map is suitable which should give expected O(1) at the cost
  // of some memory overhead - it's basically a hash-table (i think). I'm not yet sure how to 
  // implement it best. ...tbc...

  /** Adds a record for an sfz-opcode together with some information about the opcode to our 
  database to make it available for later lookup. */
  void addOpcode(Opcode op, OpcodeFormat type, const std::string& sfzStr, 
    float minValue, float maxValue, float defaultValue, 
    OpcodeType opType, OpcodeUnit unit, OpcodeSpec spec);

  /** Adds a filter type enum index with its associated sfz-string to our lookup table for later
  lookup. */
  void addFilterType(FilterType type, const std::string& sfzStr);


  // I think, these member function could actually be static - if so, make them so:

  /** Returns true, if the opcode is related to the filter, i.e. is cutoff, fil_type, resonance, 
  etc. These opcodes need some special rules for parsing because in sfz, only 1 filter exists which
  doesn't have any index...tbc... */
  bool isFilterRelated(Opcode op) const;

  /** From a string representing an indexed opcode such as "eq2_freq", it extracts the index (here, 
  the 2) and returns it. It also modifies the passed string to replace the "2" with the placeholder
  "N". */
  int getIndexAndReplaceByN(std::string& str) const;

  /** For historical reasons, certain opcode strings such as "cutoff" or "fil_type" are actually 
  supposed to mean "cutoff1" or "fil1_type" in the context of this extended implementation of sfz. 
  This function modifies the string in such cases to spell out the implicit index. */
  void makeImplicitIndexExplicit(std::string& str) const;

  /** For those opcodes that should support an implicit index 1, this functions deletes the 
  explicitly spelled out "1", if present. */
  void makeExplicitIndexImplicit(std::string& str) const;



  /** Structure for one record in our little database of sfz opcodes. */
  struct OpcodeEntry
  {                           // Example    ...maybe use width opcode as example
    Opcode       op;          // Opcode::volume
    OpcodeFormat format;      // OpcodeFormat::Float
    std::string  text;        // "volume"
    float        minVal;      // -144
    float        maxVal;      //  +12
    float        defVal;      //    0 default value
    //float        neutVal;   //    0 neutral value (some sfz defaults, like width, seem not to be)
    OpcodeType      dsp;         // OpcodeType::Amplifier
    OpcodeUnit   unit;        // OpcodeUnit::Decibels
    OpcodeSpec   spec;        // OpcodeSpec::Sfz_1
    //std::string  comment;   // "Maximum in sfz spec is 6dB. We extend it to +12 dB."
  };
  std::vector<OpcodeEntry> opcodeEntries; /**< Our lookup table of records. */

  /** Simple structure for a record to associate filter type sfz-strings with their enum 
  values. */
  struct FilterTypeEntry
  {
    FilterType  typeId;   // FilterType::lpf_1p
    std::string sfzStr;  // "lpf_1p"
  };
  std::vector<FilterTypeEntry> filterTypeEntries;  /** Lookup table for the filter types. */


  std::string dummyString;
  /**< The string-ref returning functions return a reference to this, if they do not find a 
  suitable actual string to return. 
  \todo: is this still needed? I think, we don't return string-refs anymore but actual string
  objects. may be obsolete - if so, delete */


  static SfzCodeBook* instance;
  /** Sole instance of the class.  */


  SfzCodeBook();
  /**< Constructor is protected due to singleton pattern. */

};





}}      // namespaces
#endif  // #ifndef rosic_SfzCodeBook_h