#ifndef rosic_SamplerTools_h
#define rosic_SamplerTools_h

namespace rosic {
namespace Sampler {

//=================================================================================================
// Potentially generally useful utility classes and functions that are used in the sampler engine.
// The code here may later go into some more general place somewhere in the RAPT or rosic library.

/** Return codes for the setup functions. We use encodings as negative integers so we can use
them also for functions which use positive integers as valid return values. */
enum rsReturnCode
{
  success        =  -1,  //< Operation completed successfully. 
  nothingToDo    =  -2,  //< There was nothing to actually do. State was already as desired.
  memAllocFail   =  -3,  //< Memory allocation failure.
  invalidIndex   =  -4,  //< An invalid (i.e. out-of-range) index was passed.
  overload       =  -5,  //< Not enough resources available to complete action.
  notFound       =  -6,  //< A region, group, sample or whatever was not found.
  fileLoadError  =  -7,  //< A file could not be loaded (reasons: not found or failed alloc).
  notImplemented =  -8,  //< Feature not yet implemented (relevant during development).
  failed         =  -9,  //< General failure report wihtout further specification
  malformed      = -10   //< A malformed string or object was encountered.
};
// todo: make it an enum class, maybe include also return codes for inquiry functions such as for
// "unknown", etc. ...but maybe that's no good idea when we want to use it for functions which
// need to return valid integers (like, for numChannels, etc. - we could use negative numbers to
// encode such things)
// -maybe rename "success" to "completed" or "done" because "success" has actually a more general 
//   meaning: "nothingToDo" is also a kind of "success" (or maybe "workDone" or "workCompleted"
//  ...but maybe it's not a good idea to distinguish between different kinds of "success". At the 
//  call site, we want to do tests like if(returnCode == ReturnCode::success) which should pass 
//  even when there was nothing to do
// -rename the layerOverload to a general "overload" or ressourcesUsedUp or something - to make the
//  enum more genrally useful...maybe just overload
// -maybe move it out of the Sampler sub-namespace - it may be more generally useful
// -other possibly useful codes: unavailable, denied, malformed
//  maybe invalidIndex should be a more general invalidArgument - it could be an out-of-range index
//  or a malformed string or whatever
// -maybe the fileLoadError should be replaced by more specific errors: fileNotFound, 
//  fileTooLarge, fileCorrupted ...TooLarge is returned when the mem-alloc for the buffer failed,
//  Corrupted when the sfz parser couldn't parse the content. Then we can display more specific 
//  errors on the GUI...maybe the "corrupted" should be even more specific: inform, which opcodes
//  were to blame, etc.

/** Enumeration of the available modulation modes. */
enum class ModMode  // rename to ModulationMode, maybe move out of Sampler sub-namespace
{
  // Used in sfz:
  absolute,          // for ampeg
  cents,             // for fileg, pitcheg

  // My own additions:
  relative,          // not yet used
  //exponential,
  //multiplicative,
  //percent_of_nominal,  // maybe rename to percent_relative
  //percent_absolute,    // used for ampeg where nominal value is 0
  //percent_of_max,
  //raw,

  unknown
};


//=================================================================================================

/** A class for representing musical events such as note-on/off etc. Think of it as a class to
represent MIDI events but with some of its anachronistic restrictions lifted, such as the abysmal
resolution of values (typically 7 or 14 bit). The template parameter T is supposed to be either
float or double. For easy conversion and compatibility with MIDI, we still follow the (now
historical) convention that values are in the range 0..127, but now with much higher resolution
due to the floating point representation. So, in a nutshell, this is a class for MIDI events but
with higher resolution for all the values. */

template<class T>
class rsMusicalEvent
{

public:

  enum class Type
  {
    noteOn,
    noteOff,
    controlChange,
    pitchWheel,
    reset
    // ...tbc...
  };

  using uint = RAPT::rsUint32;

  rsMusicalEvent(Type eventType, T value1, T value2, uint sampleTime = 0)
    : type(eventType), val1(value1), val2(value2), time(sampleTime) {}

  Type getType() const { return type; }

  T getValue1() const { return val1; }

  T getValue2() const { return val2; }

protected:

  Type type;  // e.g. noteOn/Off, controlChange, pitchWheel
  T    val1;  // e.g. key, controller number, pitchWheelMSB
  T    val2;  // e.g. velocity, controller value, pitchWheelLSB
  uint time;  // time-stamp of the event in samples, usually relative to the buffer start

};
// maybe move the class elsewhere for more general use - maybe it should go into rapt due to the
// templatized nature

//=================================================================================================

/** Implements a pool of objects from which items can be grabbed and later returned. To "grab" an
item means to retrieve a pointer to it using grabItem(). The caller can then use the item fro as 
long as it needs it and when it is finished, it returns the item into the pool by calling
repositItem(T*) with the pointer. This will make the item available again such that it may be 
handed out again to other clients (or to the same client again). */

template<class T>
class rsObjectPool
{

public:

  /** Initializes the pool by allocating the given number of items. */
  void init(int numItems);

  /** Returns a pointer to an item or a nullptr if no items are available anymore. */
  T* grabItem();

  /** Reposits an item into the pool such that it becomes available again. Returns the index in 
  our array where the item is located or -1, if the item could not be found. The latter condition
  indicates an error in the client code: you may have tried to return an item that you did 
  previously grabbed from the pool or you may have re-alloacted the pool while there were still
  items  in use from the old one. Typically, it is irrelevant for the client, at which index the
  object lives, but the error detection may be useful. */
  int repositItem(T* item);

  /** Like repositItem(T*) but accepts also pointers that were not declared as pointing to T. 
  This version is needed when the client uses a pointer to a baseclass of T. */
  int repositItem(void* item) { return repositItem((T*)item); }
  // todo: maybe templalize it on some type B and make an assertion that B is baseclass of T, if 
  // that is possible. That would be a bit more restrictive than accepting void and perhaps
  // catch more errors at compile time


  int getNumItems() const { return numUsed + numIdle; }

  int getNumUsedItems() const { return numUsed; }

  int getNumIdleItems() const { return numIdle; }



  /** Checks, if the class invariants are satisfied. */
  bool isInConsistentState();


protected:

  std::vector<T> items; 
  /**< Array of the actual items. We store them here as direct objects and hand out pointers
  to them to the client. 
  \todo: Maybe it should later use rsRetainedArray (a class currently in development in the 
  Prototypes) to enable dynamic growth of the pool without invalidation of the pointers at the 
  cost of discontiguous object storage (a cost which is only incurred when growth is actually 
  happening - the first pre-allocated chunk, the size of which can be chosen at initialization, 
  will be stored contiguously). Currently, dynamic growth is not supported. */

  std::vector<char> used;  
  /**< Parallel array of flags indicating, if item[i] is currently used by anyone. If 0 (false),
  the item is avaible and can be handed out. If true, it can't **/

  int numUsed = 0;  /**< Number of used items */
  int numIdle = 0;  /**< Number of items available for handing out. */

};
// todo: 
// -write unit tests
// -maybe switch std::vector to rsRetainedArray
// -in this implementation, the items are assumed to be indistinguishable which means it doesn't
//  matter which object is handed out to the client - they are all the same
// -maybe don't use a std::vector for the items but rather a raw pointer
// -maybe introduce another level on indirection to enable grab/reposit in constant or at least
//  log(N) time: always grab objects from the end of a list of available objects, when objects
//  are reposited, they are appended to the end...hmm...i have to think about it some more - we 
//  want to avoid having to iterate through the whole used[] array to find an idle item and we also
//  want to avoid iterating through the list on reposit to find the location where we need to set
//  the used flag back to false

template<class T>
void rsObjectPool<T>::init(int numItems)
{
  RAPT::rsAssert(numUsed == 0 && isInConsistentState());
  items.resize(numItems);
  used.resize(numItems);  // will be initialized to zero (that's important!)
  numUsed = 0;
  numIdle = (int) used.size();
}

template<class T>
T* rsObjectPool<T>::grabItem()
{
  for(size_t i = 0; i < used.size(); i++) {
    if(!used[i]) {
      used[i] = 1;
      numUsed++;
      numIdle--;
      return &(items[i]); }}
  return nullptr;
}

template<class T>
int rsObjectPool<T>::repositItem(T* item)
{
  for(size_t i = 0; i < items.size(); i++) {
    if(item == &(items[i]) ) {
      used[i] = 0;
      numUsed--;
      numIdle++;
      return (int)i; }}
  return -1;
}
// -grab/reposit have currently linear complexity in the number of stored items. ToDo: try to 
//  reduce it to log or even const - keep this class then as prototype
// -maybe use range-based loops to later make the switch to rsRetainedArray easier

template<class T>
bool rsObjectPool<T>::isInConsistentState()
{
  bool ok = items.size() == used.size();
  int count = 0;
  for(const auto& it : used) {
    if( it != 0 )
      count++; }
  ok &= count == numUsed;
  ok &= numIdle == (int) used.size() - numUsed;
  return ok;
}

template<class T>
void rsSetupPointers(std::vector<T>& objectArray, std::vector<T*>& pointerArray)
{
  pointerArray.resize(objectArray.size());
  for(size_t i = 0; i < objectArray.size(); i++)
    pointerArray[i] = &objectArray[i];
}
// move to library, maybe use it in the production code, too

























//=================================================================================================
// String processing:

// Some helpers for string parsing:
void rsRemoveLineComments(std::string& str, char commentStart);
void rsReplaceCharacter(std::string& str, char oldChar, char newChar);
void rsRemoveRepeats(std::string& s, char c);

/** Returns true, iff character c is a digit, i.e. one of 0,1,2,3,4,5,6,7,8,9. */
inline bool isDigit(char c) { return c >= '0' && c <= '9'; }

inline int charToInt(char c)
{
  RAPT::rsAssert(isDigit(c));
  return c - '0';
}

/** Returns the index of the first digit that occurs in the string. For example, if 
str = "Freq=3520Hz", it would return 5, i.e. the index of the digit 3. If no digit is found in the
string, it returns -1. */
int firstDigit(const std::string& str);

/** Returns the index of the last digit in a sequence of digits (i.e. in a number) starting at the 
given startIndex. */
int lastDigitInSeq(const std::string& str, int startIndex);  // remove the "InSeq" from name

/** Returns the index of the last character in the string that is not a digit. For example, in the 
string lfo13, it would return 2, the index of the letter o. */
int lastNonDigit(const std::string& str, int startIndex = 0);

/** Returns the index at which the unit suffix starts. If the is no suffix, it returns 
str.length()-1. If the whole string is a suffix, it returns 0. */
int suffixStart(const std::string& str);

int parseNaturalNumber(const std::string& str, int startIndex, int endIndex);

/** Tries to parse the given string as floating point number and returns the result. If parsing 
fails, i.e. the given string does not represent a valid floating point number, it will return 
positive NaN. */
float rsStringToFloat(const std::string& str);

/** Converts the given floating point number into a string. */
inline std::string rsFloatToString(float v) { return std::to_string(v); }

/** Returns true, iff the given string str starts with the given possible initial section. For 
example, rsStartsWith("abcd", "ab") returns true because "ab" is indeed an initial section of 
"abcd". */
bool rsStartsWith(const std::string& str, const std::string& possibleInitialSection);
// needs unit tests


//=================================================================================================
// Container tools:

/*
template<class T> // Grow v by given amount
inline void rsGrow(std::vector<T>& v, size_t amount = 1)
{
  v.resize(v.size() + amount);
}
template<class T> // Shrink v by given amount
inline void rsShrink(std::vector<T>& v, size_t amount = 1)
{
  RAPT::rsAssert(v.size() >= amount);
  v.resize(v.size() - amount);
}
template<class T> // Pointer to last element in v
inline T* rsLastPointer(std::vector<T>& v)
{
  RAPT::rsAssert(!v.empty());
  return &v[v.size()-1];
}
template<class T> // Pointer to last element, shrink by 1
inline T* rsGetLastPtrAndShrink(std::vector<T>& v)
{
  if(v.empty())
    return nullptr;
  T* p = rsLastPointer(v);
  rsShrink(v);
  return p;
}
*/
// Not used anymore and I'm not sure, if they are useful enough in general to include in the 
// library.


}
}

#endif