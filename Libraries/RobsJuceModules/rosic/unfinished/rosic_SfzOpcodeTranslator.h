#ifndef rosic_SfzOpcodeTranslator_h
#define rosic_SfzOpcodeTranslator_h
namespace rosic { namespace Sampler {




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

  // ToDo:
  //const std::string& opcodeToString(Opcode op);
  //Opcode stringToOpcode(const std::string& str);

  //SignalProcessorType opcodeToProcessor(Opcode op);

  //const std::string& filterTypeToString(FilterType ft);
  //FilterType stringToFilterType(const std::string& str);

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

};


}}      // namespaces
#endif  // #ifndef rosic_SfzOpcodeTranslator_h