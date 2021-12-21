#ifndef rosic_SfzOpcodeTranslator_h
#define rosic_SfzOpcodeTranslator_h
namespace rosic { namespace Sampler {

//-------------------------------------------------------------------------------------------------
/** Enumeration of possible types of settings. These types correspond to the opcodes defined
in the sfz specification. */

enum class Opcode
{
  Unknown = 0,

  // Input controls:
  LoKey, HiKey, LoVel, HiVel,
  ControllerRangeLo, ControllerRangeHi, PitchWheelRange,  // 

  // Muted: convenient to switch regions or groups off wihthout removing them - check if 
  // sfz has such a thing

  // Pitch:
  PitchKeyCenter, Transpose, Tune,

  // Amplitude:
  Volume, Pan, PanRule,
  AmpEnvAttack, AmpEnvDecay, AmpEnvSustain, AmpEnvRelease,
  // ToDo: Width, Position

  // Player:
  Delay, Offset,
  // ToDo: loop-stuff

  // Filter:
  FilterType, FilterCutoff, FilterResonance,

  // Some of my own extensions
  // Distortion:
  DistShape, DistDrive,  // DistGain...may be redundant with Volume

  NumTypes 
  //...tbc...
};
// The underlying integers for the opcodes in this enum must start at zero and increment by one.
// Dont do something like FilterCutoff = 1000 etc. They must be usable as array indices and we dont
// want to waste space

// maybe don't capitalize first letter - make it conistent with other (newer) enums in the 
// library. see community stadards recommendations...


//-------------------------------------------------------------------------------------------------
/** Enumeration of the different data formats of the values of the opcode. */

enum class OpcodeType  // maybe rename to OpcodeFormat
{
  Unknown = 0,

  Float,
  Integer,
  Boolean,
  String,

  NumTypes
};

//-------------------------------------------------------------------------------------------------
/** Enumeration of the physical units that apply to the value of the opcode. */

enum class OpcodeUnit
{
  Unknown = 0,

  Hertz,
  Semitones,
  Cents,
  Decibels,
  Seconds,
  Milliseconds,
  Beats,
  BeatsOrSeconds,
  Text,
  MidiKey,           // i think, it can be a number or text
  Percent,
  Index,
  RawInt,
  RawFloat,

  NumUnits
};
//   // ToDo: have a field for the unit: Hz, cents, semitones, noteNumber, dB, sec, beats, .

//-------------------------------------------------------------------------------------------------
/** Enumeration of the different specifications in which the respective opcode was defined. */

enum class OpcodeSpec
{
  Unknown,

  Sfz,
  Sfz2,
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

enum class SignalProcessorType  // rename to DspType or ProcessorType
{
  Unknown,

  SamplePlayer,

  // The modulators:
  // AmpEnv, FilterEnv, PitchEnv, AmpLFO, ...

  // The actual DSP processors:
  Filter,
  WaveShaper,

  NumDsps
};

//-------------------------------------------------------------------------------------------------
/** Enumeration of the different filter types that are available. */

enum class FilterType
{
  off = 0, lp_6, lp_12, hp_6, hp_12, bp_6_6, br_6_6,

  numFilterTypes
};

//-------------------------------------------------------------------------------------------------

enum class PanRule
{
  linear, sinCos,

  numPanRules
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

//=================================================================================================

/** A class to translate between the sfz opcodes in their string representation and their enum 
values and some additional related "translations" such as the mapping of opcodes to the type of
signal processor to which they apply, etc. Essentially, it does many of the conversions and 
mappings that are required within the sampler engine. An object of this class shall be created 
once when the sampler is created - maybe as global object or as singleton. The creating may be
somewhat expensive, so it should be done only once and then the object should remain available
for the whole lifetime of the sampler engine....tbc... */

class SfzOpcodeTranslator // maybe rename to SfzDictionary, SfzTranslator
{

public:

  SfzOpcodeTranslator();


  const std::string& opcodeToString(Opcode op);
  Opcode stringToOpcode(const std::string& str);
  SignalProcessorType opcodeToProcessor(Opcode op);


  // ToDo:
  //const std::string& filterTypeToString(FilterType ft);
  //FilterType stringToFilterType(const std::string& str);
  // add documentation

  // ...etc.

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


  /** Adds an opcode to our database to make it available for later lookup. */
  void addOpcode(Opcode op, OpcodeType type, const std::string& sfzStr, 
    float minValue, float maxValue, float defaultValue, SignalProcessorType dspType);

  /** Structure for one record in our little database or lookup table. */
  struct OpcodeEntry
  {                       // Examples
    Opcode       op;      // Opcode::Cutoff
    OpcodeType   type;    // OpcodeType::Float
    std::string  str;     // cutoff
    float        minVal;  // 20?
    float        maxVal;  // 20000?
    float        defVal;  // 1000?
    SignalProcessorType dsp;
  };
  std::vector<OpcodeEntry> opcodeEntries; /**< Our lookup table of records. */


  std::string dummyString;
  /**< The string-ref returning functions return a reference to this, if they do not find a 
  suitable actual string to return. */

};


}}      // namespaces
#endif  // #ifndef rosic_SfzOpcodeTranslator_h