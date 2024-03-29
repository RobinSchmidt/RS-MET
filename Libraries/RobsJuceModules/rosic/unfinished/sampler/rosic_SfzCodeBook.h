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
  FilEnvDelay, FilEnvStart, FilEnvAttack, FilEnvHold, FilEnvDecay, FilEnvSustain,
  FilEnvRelease, FilEnvDepth, FilEnvVel2Delay, FilEnvVel2Attack, FilEnvVel2Hold, 
  FilEnvVel2Decay, FilEnvVel2Sustain, FilEnvVel2Release, FilEnvVel2Depth,

  // Filter LFO:
  FilLfoDelay, FilLfoFade, FilLfoFreq, FilLfoDepth, FilLfoDepthCtrlN, 
  FilLfoDepthChanAft, FilLfoDepthPolyAft, FilLfoFreqCtrlN, FilLfoFreqChanAft, FilLfoFreqPolyAft,

  // Amplifier:
  /*Volume, Pan, Width, Position, */  // replaced by volumeN, etc.
  /*AmpKeyTrack, AmpKeyCenter, AmpVelTrack,*/ AmpVelCurveN, AmpRandom,
  RelTrigDecay, Output, GainCtrlN,
  FadeInLoKey, FadeInHiKey, FadeOutLoKey, FadeOutHiKey, FadeCurveKey,
  FadeInLoVel, FadeInHiVel, FadeOutLoVel, FadeOutHiVel, FadeCurveVel,
  FadeInLoCtrlN, FadeInHiCtrlN, FadeOutLoCtrlN, FadeOutHiCtrlN, FadeCurveCtrl,

  // new amp-opcodes - will make the above ones obsolete:
  volumeN, panN, widthN, positionN, 
  ampN_keytrack, ampN_keycenter, ampN_veltrack,
  // ampN_velcurve_X /* ? */, ampN_random,

  // Amplifier Envelope:
  AmpEnvDelay, AmpEnvStart, AmpEnvAttack, AmpEnvHold, AmpEnvDecay, AmpEnvSustain,
  AmpEnvRelease, /*AmpEnvDepth,*/ AmpEnvVel2Delay, AmpEnvVel2Attack, AmpEnvVel2Hold, 
  AmpEnvVel2Decay, AmpEnvVel2Sustain, AmpEnvVel2Release, /*AmpEnvVel2Depth,*/
  AmpEnvDelayCtrlN, AmpEnvStartCtrlN, AmpEnvAttackCtrlN, AmpEnvHoldCtrlN, AmpEnvDecayCtrlN, 
  AmpEnvSustainCtrlN, AmpEnvReleaseCtrlN, 
  // depth-parameters do not exist in sfz for amp-env. wouldn't make much sense, i guess

    // Amplifier LFO:
  AmpLfoDelay, AmpLfoFade, AmpLfoFreq, AmpLfoDepth, AmpLfoDepthCtrlN, 
  AmpLfoDepthChanAft, AmpLfoDepthPolyAft, AmpLfoFreqCtrlN, AmpLfoFreqChanAft, AmpLfoFreqPolyAft,

  // Equal�zer:
  eqN_freq, eqN_freqccX, eqN_vel2freq, eqN_bw, eqN_bwccX, eqN_gain, eqN_gainccX, eqN_vel2gain,

  // Effects:
  effN,         // effect send levels in percent (eff1: reverb, eff2: chorus)

  // SFZ 2.0:
  // fil2_type, cutoff2, etc. ...generalize to filN_type, cutoffN, etc. - maybe for translations
  // and lookup, we'll need a special rule for filter-related parameters that leaves out the 1,
  // i.e. fil1_type, cutoff1, etc. should not be produced or accepted - instead, we'll use 
  // fil_type, cutoff, etc. because that's what sfz 1 has. it doesn't have a concept of multiple
  // filters

  // egN_timeX, egN_levelX
  // lfoN_freq, lfoN_delay, lfoN_fade, lfoN_phase, lfoN_wave, lfoN_volume

  // ARIA:
  PanLaw,




  // My own extensions - preliminary, for experimentation. Before defining extensions, check what 
  // is already there in SFZ 2.0 and in other sfz engines with extensions. Try to be compatible 
  // with the  largest possible range of other engines. Discuss extensions on KVR before 
  // implementing them in production code:
    
  distortN_shape, distortN_drive, distortN_dc, 
  //distortN_shape, distortN_drive, distortN_dc, // DistGain...may be redundant with Volume


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


  numTypes
};

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

enum class OpcodeFormat  // maybe rename to OpcodeFormat
{
  Unknown = 0,

  Boolean,
  Natural,  // a.k.a. unsigned int
  Integer,
  Float,
  String,    // rename to Text

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
  Aria,
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

enum class DspType  // rename to DspType or ProcessorType or OpcodeType
{
  Unknown,

  SamplePlayer,

  // The modulators:
  // AmpEnv, FilterEnv, PitchEnv, AmpLFO, ...

  // The actual DSP processors:
  _TagDspsStart,
  Amplifier,
  Filter,
  Equalizer,
  WaveShaper,
  _TagDspsEnd,

  // Opcodes for controlling key- and velocity tracking:
  //Tracking,           // no - we don't need this - the tracking parameters bleong into the same
  //_TagTackingStart, // category as their target parameters
  //_TagTrackingEnd,



  NumDsps
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

  /** Returns the type of signal processor to which the given opcdoe applies. */
  DspType opcodeToProcessor(Opcode op);

  /** Returns the default value for the given opcode as floating point number. If the format of the
  value is integer or an enum value, you'll need to convert the returned value to int and then
  possibly to the enum. Some opcodes contain an index like eq2_freq. For these, you need to pass 
  that index too because the default value may actually depend on that index. For example, the 
  defaults are 50,500,5000 respectively for eq1_freq, eq2_freq, eq3_freq. If indexing is not 
  applicable to the given opcode, you should pass index = -1. */
  float opcodeDefaultValue(Opcode op, int index);

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



  //-----------------------------------------------------------------------------------------------
  // \name Inquiry

  /** Returns true iff the given opcode applies to a DSP in the DSP chain like volume, cutoff, 
  etc. but not pitch_keycenter or tune. */
  bool isDspSetting(Opcode op)
  {
    DspType dspType = opcodeToProcessor(op);
    return dspType > DspType::_TagDspsStart && dspType < DspType::_TagDspsEnd;
  }

  /** Returns true iff the given opcode applies to the sample playback source such as tune, 
  delay, offset, etc. */
  bool isPlayerSetting(Opcode op)
  {
    DspType dspType = opcodeToProcessor(op);
    return dspType == DspType::SamplePlayer;
  }

  // todo: isModulationSetting, isModulatorSetting  (modulatiON settings define mod-connections,
  // modulatOR settings parameters of the modulators), hasIndex..but this may be expensive to 
  // implement - we may go through the opcode string and try to find an 'N'. Not too bad for 
  // non-realtime code but perhaps not good to call in realtime. Or, maybe, if we need that at 
  // realtime add another field to the OpcodeEntry record. maybe a set of flags


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
    DspType dspType, OpcodeUnit unit, OpcodeSpec spec);

  /** Adds a filter type enum index with its associated sfz-string to our lookup� table for later
  lookup. */
  void addFilterType(FilterType type, const std::string& sfzStr);


  // these could actually be static:

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
    float        defVal;      //    0
    //float        neutVal;   //    0 neutral value (some sfz defaults, like width, seem not to be)
    DspType      dsp;         // DspType::Amplifier
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
  objects */


  static SfzCodeBook* instance;
  /** Sole instance of the class.  */


  SfzCodeBook();
  /**< Constructor is protected due to singleton pattern. */

};


}}      // namespaces
#endif  // #ifndef rosic_SfzCodeBook_h