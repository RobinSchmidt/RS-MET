namespace rosic { namespace Sampler {

SfzOpcodeTranslator::SfzOpcodeTranslator()
{
  opcodeEntries.reserve((int)Opcode::NumTypes);

  using OC = Opcode;
  using OT = OpcodeType;
  using SP = SignalProcessorType;

  addOpcode(OC::FilterCutoff, OT::Float, "cutoff", 20.f, 20000.f, 1000.f, SP::Filter);
  // verify min/max/def (i made them up!)



}

void SfzOpcodeTranslator::addOpcode(Opcode op, OpcodeType type, const std::string& sfzStr,
  float minVal, float maxVal, float defVal, SignalProcessorType dspType)
{
  // ah - that will lead to the opcodes being stored in the order we add them. What we should do 
  // instead is put it at the place determined by (int)op. We need may to resize the array first

  opcodeEntries.push_back(
    OpcodeEntry({ op, type, sfzStr, minVal, maxVal, defVal, dspType }));

  int dummy = 0;
}

const std::string& SfzOpcodeTranslator::opcodeToString(Opcode op)
{


  RAPT::rsError("Unknown opcode in SfzOpcodeTranslator::opcodeToString");
  return std::string();
}

Opcode SfzOpcodeTranslator::stringToOpcode(const std::string& str)
{



  RAPT::rsError("Unknown opcode in SfzOpcodeTranslator::stringToOpcode");
  return Opcode::Unknown;
}


}}