namespace rosic { namespace Sampler {

SfzOpcodeTranslator::SfzOpcodeTranslator()
{
  // On construction, we build our database (maybe factor out):
  using OC = Opcode;
  using OT = OpcodeType;
  using SP = SignalProcessorType;
  opcodeEntries.resize((int)Opcode::NumTypes);

  // Sample playback:
  addOpcode(OC::LoKey, OT::Integer, "lokey", 0, 127,   0, SP::SamplePlayer);
  addOpcode(OC::HiKey, OT::Integer, "hikey", 0, 127, 127, SP::SamplePlayer);
  addOpcode(OC::LoVel, OT::Integer, "lovel", 0, 127,   0, SP::SamplePlayer);
  addOpcode(OC::HiVel, OT::Integer, "hivel", 0, 127, 127, SP::SamplePlayer);

  // Pitch:


  // Amplitude:



  // Filter:
  addOpcode(OC::FilterCutoff, OT::Float, "cutoff", 20.f, 20000.f, 1000.f, SP::Filter);
  // verify min/max/def (i made them up!)

  int dummy = 0;
}

template<class T>
inline void rsEnsureSize(std::vector<T>& v, size_t s)
{
  if(v.size() < s)
    v.resize(s);
} // maybe move to rapt
void SfzOpcodeTranslator::addOpcode(Opcode op, OpcodeType type, const std::string& sfzStr,
  float minVal, float maxVal, float defVal, SignalProcessorType dspType)
{
  int i = (int)op;
  rsEnsureSize(opcodeEntries, size_t(i+1));
  opcodeEntries[i] = OpcodeEntry({ op, type, sfzStr, minVal, maxVal, defVal, dspType });
  // Actually, storing the "op" is redundant because it's implicitly given by the array index, so
  // maybe remove that field...but maybe it's useful in other contexts
}

const std::string& SfzOpcodeTranslator::opcodeToString(Opcode op)
{
  if((int)op < 0 || (int)op >= (int)opcodeEntries.size()) {
    RAPT::rsError("Unknown opcode in SfzOpcodeTranslator::opcodeToString");
    return dummyString; 
  }
  return opcodeEntries[(int)op].str;
}
// needs test

Opcode SfzOpcodeTranslator::stringToOpcode(const std::string& str)
{
  for(int i = 0; i < opcodeEntries.size(); i++)
    if(opcodeEntries[i].str == str)
      return opcodeEntries[i].op;    // op should be equal to i
  RAPT::rsError("Unknown opcode in SfzOpcodeTranslator::stringToOpcode");
  return Opcode::Unknown;

  // This has currently linear complexity in the number of opcodes. Maybe bring this down to at
  // most O(log(N)) by maintaining a map of indices into the opcodeEntries array that is sorted
  // lexicographically according to the opcode string. I'm not sure, if it's worth it though.
  // This is called only on patch loading and maybe it's fast enough as is. We'll see.
}

SignalProcessorType SfzOpcodeTranslator::opcodeToProcessor(Opcode op)
{
  if((int)op < 0 || (int)op >= (int)opcodeEntries.size()) {
    RAPT::rsError("Unknown opcode in SfzOpcodeTranslator::opcodeToProcessor");
    return SignalProcessorType::Unknown; 
  }
  return opcodeEntries[(int)op].dsp;
}

}}

/*

ToDo:
-Write a unit test that loops through all opcodes and translates them to a string and back. It 
 should also ensure that the op member of the record is equal to the array index. Wrtie similar
 tests for back-and-forth conversion of the string-type parameters.
-Maybe write a little program that generates a string which lists all opcodes with their 
 properties, perhaps sorted by name, dsp-type, etc. This could be used to generate a little 
 textfile as reference manual which would be a convenient thing to have. Maybe it could even be 
 member function here: generateReferenceManual() or something. It could perhaps take a couple of
 parameters to select formatting and sorting options

*/