namespace rosic { namespace Sampler {

SfzOpcodeTranslator::SfzOpcodeTranslator()
{
  // On construction, we build our database (maybe factor out):
  using OC = Opcode;
  using OF = OpcodeFormat;
  using SP = DspType;
  using OU = OpcodeUnit;
  using OS = OpcodeSpec;
  opcodeEntries.resize((int)Opcode::NumTypes);

  // Sample playback:
  addOpcode(OC::LoKey, OF::Integer, "lokey", 0, 127,   0, SP::SamplePlayer, OU::MidiKey, OS::Sfz_1);
  addOpcode(OC::HiKey, OF::Integer, "hikey", 0, 127, 127, SP::SamplePlayer, OU::MidiKey, OS::Sfz_1);
  addOpcode(OC::LoVel, OF::Integer, "lovel", 0, 127,   0, SP::SamplePlayer, OU::RawInt,  OS::Sfz_1);
  addOpcode(OC::HiVel, OF::Integer, "hivel", 0, 127, 127, SP::SamplePlayer, OU::RawInt,  OS::Sfz_1);

  // Pitch:
  addOpcode(OC::PitchKeyCenter, OF::Integer, "pitch_keycenter", -127, 127, 60, SP::SamplePlayer, OU::MidiKey,   OS::Sfz_1);
  addOpcode(OC::Transpose,      OF::Integer, "transpose",       -127, 127,  0, SP::SamplePlayer, OU::Semitones, OS::Sfz_1);
  addOpcode(OC::Tune,           OF::Integer, "tune",            -100, 100,  0, SP::SamplePlayer, OU::Cents,     OS::Sfz_1);

  // Amplitude:

  // Filter:
  addOpcode(OC::FilterCutoff, OF::Float, "cutoff", 20.f, 20000.f, 1000.f, SP::Filter, OU::Hertz, OS::Sfz_1);
  // verify min/max/def (i made them up!)


  // ToDo: try to make the calls shorter by:
  // -define a lambda-function add to use in place of addOpcode
  // -using abbreviations for OF::Integer, etc.
  // -using abbreviations for SP::SamplePlayer etc. - if we collect the calls belonging to the same
  //  DspType we can use variable dsp for that which we re-assign before the next round of calls

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

const std::string& SfzOpcodeTranslator::opcodeToString(Opcode op)
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

}}

/*

ToDo:
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

*/