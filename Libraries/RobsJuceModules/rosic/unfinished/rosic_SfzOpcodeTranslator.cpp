namespace rosic { namespace Sampler {

SfzOpcodeTranslator* SfzOpcodeTranslator::instance = nullptr;

SfzOpcodeTranslator::SfzOpcodeTranslator()
{
  // On construction, we build our database (maybe factor out):
  using OC = Opcode;
  using OF = OpcodeFormat;
  using SP = DspType;
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

  OS Sfz1 = OS::Sfz_1;

  // Player response constraints (aka "Input Control" in the sfz doc):
  SP dsp = DspType::SamplePlayer;
  add(OC::LoKey, Nat, "lokey", 0, 127,   0, dsp, OU::MidiKey, Sfz1);
  add(OC::HiKey, Nat, "hikey", 0, 127, 127, dsp, OU::MidiKey, Sfz1);
  add(OC::LoVel, Nat, "lovel", 0, 127,   0, dsp, OU::RawInt,  Sfz1);
  add(OC::HiVel, Nat, "hivel", 0, 127, 127, dsp, OU::RawInt,  Sfz1);

  // Player playback settings:
  add(OC::Delay,  Flt, "delay",  0,        100, 0, dsp, OU::Seconds, Sfz1);
  add(OC::Offset, Nat, "offset", 0, 4294967296, 0, dsp, OU::Samples, Sfz1);
  add(OC::Count,  Nat, "count",  0, 4294967296, 0, dsp, OU::RawInt,  Sfz1);

  add(OC::LoopMode, Txt, "loop_mode",  (float)LoopMode::Unknown + 1.f, 
    (float)LoopMode::numModes - 1.f, (float)LoopMode::no_loop, dsp, OU::Text, Sfz1);
  // In sfz, the default value depends on whether or not the sample file itself defines a loop,
  // if yes: continuous, if not: no loop

  add(OC::LoopStart, Nat, "loop_start", 0, 4294967296, 0, dsp, OU::Samples, Sfz1);
  // If not specified, the defined value in the sample will be used, if any, zero otherwise.

  add(OC::LoopEnd,   Nat, "loop_end",   0, 4294967296, 0, dsp, OU::Samples, Sfz1);
  // If not specified, the defined value in the sample will be used, if any.

  // Player pitch:
  add(OC::Transpose,      Int, "transpose",       -127, 127,  0, dsp, OU::Semitones, Sfz1);
  add(OC::Tune,           Int, "tune",            -100, 100,  0, dsp, OU::Cents,     Sfz1);
  add(OC::PitchKeyCenter, Nat, "pitch_keycenter",    0, 127, 60, dsp, OU::MidiKey,   Sfz1); 
  // For pitch_keycenter, the spec actually says, -127..127 for the range but also C-1...G9 and 
  // C-1 would map to 0. So, that's probably a mistake in the sfz documentation and they actually 
  // mean 0..127. I assume this here but should check what other implementation do. In sfz files, 
  // keycenter can be given numerically or textually as e.g. c#2, so the parser should support 
  // midi-note to number conversion.


  // Filter:
  dsp = DspType::Filter;
  add(OC::FilType, Txt, "fil_type", (float)FilterType::Unknown + 1.f, 
    (float)FilterType::numFilterTypes - 1.f, (float)FilterType::lp_12, dsp, OU::Text, Sfz1); 
  // sfz default is lpf_2p - maybe rename our enum values to be consistent with sfz

  add(OC::Cutoff, Flt, "cutoff", 0.f, 22050.f, 22050.f, dsp, OU::Hertz, Sfz1);
  // Range is 0..fs/2, default is: filter disabled, so perhaps, the default should depend on the 
  // selected type: fs/2 for a lowpass, 0 for a highpass - figure out what sfz+ does

  add(OC::Resonance, Flt, "resonance", 0.0f, +40.f, 0.f, dsp, OU::Decibels, Sfz1);
  // In sfz, the resonance is adjusted in terms of the resonance gain in dB. If zero dB is the 
  // minimum, does that mean there is always some resonance? Because without resonance, the gain
  // at cutoff would be -3.01 dB, I think. Does the resonance parameter give the gain at the cutoff
  // freq or at the peak? -> figure out experimentally or by looking at other implementations. 
  // Maybe I should derive a formula for Q in terms of the resonance gain. The formula for the peak
  // gain would probably be more complicated.

  // Player amplifier:
  dsp = DspType::SamplePlayer;  // use Amplifier later
  add(OC::Volume, Flt, "volume", -144.f,   +6.f, 0.f, dsp, OU::Decibels, Sfz1);
  add(OC::Pan,    Flt, "pan",    -100.f, +100.f, 0.f, dsp, OU::RawFloat, Sfz1);



  // This is very very preliminary - don't use it yet to define actual instruments - it's going
  // to change:
  dsp = DspType::WaveShaper;
  OS RsMet = OS::RsMet;
  add(OC::DistShape, Nat, "dist_shape", 0.f, 0.f, 0.f, dsp, OU::RawInt,   RsMet);
  add(OC::DistDrive, Flt, "dist_drive", 0.0, 8.0, 1.0, dsp, OU::RawFloat, RsMet);




  // ToDo: 
  // -PanLaw
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

  // Filter types:
  // SFZ 1: lpf_1p, hpf_1p, lpf_2p, hpf_2p, bpf_2p, brf_2p
  // SFZ 2: lpf_4p, hpr_4p, lpf_6p, hpf_6p, bpf_1p, brf_1p, apf_1p, pkf_2p, lpf_2p_sv, hpf_2p_sv, 
  //        bpf_2p_sv, brf_2p_sv, comb, pink.
  // https://sfzformat.com/legacy/
  // https://www.linuxsampler.org/sfz/
  // http://ariaengine.com/forums/index.php?p=/discussion/4389/arias-custom-opcodes/p1
  // ...it has 6-pole filters! :-O can we realize that with the current filter implementation 
  // without increasing its memory footprint? maybe using 3 equal biquads in DF2 or TDF1?


  int dummy = 0;
}

template<class T>
inline void rsEnsureSize(std::vector<T>& v, size_t s)
{
  if(v.size() < s)
    v.resize(s);
} // maybe move to rapt
void SfzOpcodeTranslator::addOpcode(Opcode op, OpcodeFormat type, const std::string& sfzStr,
  float minVal, float maxVal, float defVal, DspType dspType, OpcodeUnit unit, OpcodeSpec spec)
{
  int i = (int)op;
  rsEnsureSize(opcodeEntries, size_t(i+1));
  opcodeEntries[i] = 
    OpcodeEntry({ op, type, sfzStr, minVal, maxVal, defVal, dspType, unit, spec });
  // Actually, storing the "op" is redundant because it's implicitly given by the array index, so
  // maybe remove that field...but maybe it's useful in other contexts
}

const std::string& SfzOpcodeTranslator::opcodeToString(Opcode op) const
{
  if((int)op < 0 || (int)op >= (int)opcodeEntries.size()) {
    RAPT::rsError("Unknown opcode in SfzOpcodeTranslator::opcodeToString");
    return dummyString; 
  }
  return opcodeEntries[(int)op].text;
}
// needs test

Opcode SfzOpcodeTranslator::stringToOpcode(const std::string& str)
{
  for(int i = 0; i < opcodeEntries.size(); i++)
    if(opcodeEntries[i].text == str)
      return opcodeEntries[i].op;    // op should be equal to i
  RAPT::rsError("Unknown opcode in SfzOpcodeTranslator::stringToOpcode");
  return Opcode::Unknown;

  // This has currently linear complexity in the number of opcodes. Maybe bring this down to at
  // most O(log(N)) by maintaining a map of indices into the opcodeEntries array that is sorted
  // lexicographically according to the opcode string. I'm not sure, if it's worth it though.
  // This is called only on patch loading and maybe it's fast enough as is. We'll see.
}

DspType SfzOpcodeTranslator::opcodeToProcessor(Opcode op)
{
  if((int)op < 0 || (int)op >= (int)opcodeEntries.size()) {
    RAPT::rsError("Unknown opcode in SfzOpcodeTranslator::opcodeToProcessor");
    return DspType::Unknown; 
  }
  return opcodeEntries[(int)op].dsp;
}

float SfzOpcodeTranslator::opcodeDefaultValue(Opcode op)
{
  if((int)op < 0 || (int)op >= (int)opcodeEntries.size()) {
    RAPT::rsError("Unknown opcode in SfzOpcodeTranslator::opcodeDefaultValue");
    return 0.f;
  }
  return opcodeEntries[(int)op].defVal;
}

SfzOpcodeTranslator* SfzOpcodeTranslator::getInstance()
{
  RAPT::rsAssert(instance != nullptr);
  // Client code is supposed to explicitly create the singleton instance using createInstance() 
  // before using it. It should also clean up by calling deleteInstance(), when the object is not 
  // needed anymore. We need this explicit lifetime management (in particular, the clean up) of 
  // the singleton to prevent false positives from the memory leak checker...

  if(instance == nullptr)  // ...yeah, ok - just in case...but it's really cleaner to do an 
    createInstance();      // explicit creation somewhere before usage.
  return instance;
}

void SfzOpcodeTranslator::createInstance()
{
  RAPT::rsAssert(instance == nullptr);
  // Don't create a new instance before deleting the old one. That's a memory leak!

  instance = new SfzOpcodeTranslator;
}

void SfzOpcodeTranslator::deleteInstance()
{
  delete instance;
  instance = nullptr;
}


}}

/*

ToDo:
-Turn it into a Singleton: 
 -add static member of type SfzOpcodeTranslator*, init it to nullptr -> done
 -add a static getInstance() method returning a pointer to our SfzOpcodeTranslator, creating the
  object first, if not yet done -> done
 -provide a cleanup function that deallocates the object -> done
 -create may be called from the constructor and cleanup from the destructor of rsSamplerEngine.
  ...hmm...but that's not good when multiple instance of the sampler engine are openened -> maybe
  the sampler engine shopuld have an instance counter and only the last of them deallocates?
 -make constructor and assignment operator protected
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


SFZ - Resources:
https://sfzformat.com/legacy/   opcode reference
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