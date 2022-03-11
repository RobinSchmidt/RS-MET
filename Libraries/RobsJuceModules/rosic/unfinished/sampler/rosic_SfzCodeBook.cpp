namespace rosic { namespace Sampler {

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

  OF Nat  = OF::Natural;
  OF Int  = OF::Integer;
  OF Flt  = OF::Float;
  OF Txt  = OF::String;

  OS Sfz1  = OS::Sfz_1;
  OS Sfz1e = OS::Sfz_1_E;

  // Player response constraints (aka "Input Control" in the sfz doc):
  SP dsp = OpcodeType::SamplePlayer;
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
  add(OC::volumeN,   Flt, "volumeN",   -144.f,   +6.f,   0.f, dsp, OU::Decibels, Sfz1e);
  add(OC::panN,      Flt, "panN",      -100.f, +100.f,   0.f, dsp, OU::RawFloat, Sfz1e);
  add(OC::widthN,    Flt, "widthN",    -100.f, +100.f, 100.f, dsp, OU::Percent,  Sfz1e);
  add(OC::positionN, Flt, "positionN", -100.f, +100.f,   0.f, dsp, OU::Percent,  Sfz1e);
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

  dsp = OpcodeType::AmpLfo;
  add(OC::amplfo_freq,  Flt, "amplfo_freq",    0.0f,  20.f, 0.f, dsp, OU::Hertz,    Sfz1e);
  add(OC::amplfo_depth, Flt, "amplfo_depth", -10.0f, +10.f, 0.f, dsp, OU::Decibels, Sfz1e);
  // Or should the amplfo_depth be of type ModulationRouting, ModConnection?




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
  // https://www.linuxsampler.org/sfz/
  // http://ariaengine.com/forums/index.php?p=/discussion/4389/arias-custom-opcodes/p1
  // ...it has 6-pole filters! :-O can we realize that with the current filter implementation 
  // without increasing its memory footprint? maybe using 3 equal biquads in DF2 or TDF1?
  // ...but maybe it would be a better idea to not use one class that does all filter types but
  // instead have for each filter topology an extra class. This will reduce the memory footprint
  // when only simple filters are used in a patch and opens the possibility to later include really
  // fancy filters without blowing up the memory footprint of patches which use only simple 
  // filters


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

OpcodeType SfzCodeBook::opcodeToProcessor(Opcode op)
{
  if((int)op < 0 || (int)op >= (int)opcodeEntries.size()) {
    RAPT::rsError("Unknown opcode in SfzCodeBook::opcodeToProcessor");
    return OpcodeType::Unknown; 
  }
  return opcodeEntries[(int)op].dsp;
}

float SfzCodeBook::opcodeDefaultValue(Opcode op, int index)
{
  if((int)op < 0 || (int)op >= (int)opcodeEntries.size()) {
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

std::string SfzCodeBook::opcodeToString(Opcode op, int index) const
{
  if((int)op < 0 || (int)op >= (int)opcodeEntries.size()) {
    RAPT::rsError("Unknown opcode in SfzCodeBook::opcodeToString");
    return dummyString; }
  RAPT::rsAssert(index > 0 || index == -1, "Invalid index in SfzCodeBook::opcodeToString");
  // Either index is actually a valid index in which case must be a positive natural number or 
  // indexing doesn't apply to the given opcode in which case we use -1 as code to indicate this
  // situation.

  if(index != -1) {                              // if we are dealing with an indexed opcode...
    std::string s = opcodeEntries[(int)op].text; //   retrieve opcode template (e.g. eqN_freq)
    rsReplace(s, "N", std::to_string(index));    //   replace placeholder "N" with actual index
    makeExplicitIndexImplicit(s);                //   turn e.g. fil1_type into fil_type
    return s; }
  else
    return opcodeEntries[(int)op].text;
}

//-------------------------------------------------------------------------------------------------
// These helpers should be moved to rosic_String.h/cpp:

/** Returns true, iff character c is a digit, i.e. one of 0,1,2,3,4,5,6,7,8,9. */
bool isDigit(char c) { return c >= '0' && c <= '9'; }

int charToInt(char c)
{
  RAPT::rsAssert(isDigit(c));
  return c - '0';
}

/** Returns the index of the first digit that occurs in the string. For example, if 
str = "Freq=3520Hz", it would return 5, i.e. the index of the digit 3. If no digit is found in the
string, it returns -1. */
int firstDigit(const std::string& str)
{
  for(int i = 0; i < (int) str.size(); i++)
    if(isDigit(str[i]))
      return i;
  return -1;
}

/** Returns the index of the last digit in a sequence of digits (i.e. in a number) starting at the 
given startIndex. */
int lastDigitInSeq(const std::string& str, int startIndex)
{
  RAPT::rsAssert(isDigit(str[startIndex]));
  int i = startIndex;
  while(i < ((int)str.size())-1) {
    if(!isDigit(str[i+1]))
      return i;
    i++; }
  return i;
}
int parseNaturalNumber(const std::string& str, int startIndex, int endIndex)
{
  int num    = 0;
  int scaler = 1;
  for(int i = endIndex; i >= startIndex; i--)
  {
    num += scaler * charToInt(str[i]);
    scaler *= 10;
  }
  return num;
}
// needs unit tests

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
  RAPT::rsError("Unknown opcode in SfzCodeBook::stringToOpcode");
  return Opcode::Unknown;

  // This lookup has currently linear complexity in the number of opcodes. Maybe bring this down 
  // to at most O(log(N)) by maintaining a map of indices into the opcodeEntries array that is 
  // sorted lexicographically according to the opcode string. I'm not sure, if it's worth it 
  // though. This is called only on patch loading and maybe it's fast enough as is. We'll see.
  // Another option could be to use a hash-table as in std::map.

  // ToDo:
  // -We need a way to return the index, too!
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
}

std::string SfzCodeBook::loopModeToString(LoopMode lm)
{
  if(lm == LoopMode::no_loop)         return "no_loop";
  if(lm == LoopMode::loop_continuous) return "loop_continuous";
  if(lm == LoopMode::one_shot)        return "one_shot";
  return "unknown";
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
  case Opcode::filN_type: return (float) stringToFilterType(str);
  case Opcode::LoopMode:  return (float) stringToLoopMode(str);
  default: return std::stof(str);
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

SfzCodeBook* SfzCodeBook::getInstance()
{
  RAPT::rsAssert(instance != nullptr);
  // Client code is supposed to explicitly create the singleton instance using createInstance() 
  // before using it. It should also clean up by calling deleteInstance(), when the object is not 
  // needed anymore. We need this explicit lifetime management (in particular, the clean up) of 
  // the singleton to prevent false positives from the memory leak checker. Well, it's actually
  // a valid positive - the (GoF) textbook version of the pattern doesn't clean up.

  if(instance == nullptr)  // ...yeah, ok - just in case...but it's really cleaner to do an 
    createInstance();      // explicit creation somewhere before usage.
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


}}

/*

ToDo:
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

https://sfzformat.com/software/players/  players (also open source)
https://plugins4free.com/plugin/217/   sfz by rgcaudio

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


*/