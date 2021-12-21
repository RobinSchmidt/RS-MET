namespace rosic { namespace Sampler {

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