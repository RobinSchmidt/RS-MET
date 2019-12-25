#ifndef romos_ModuleTypeRegistry_h
#define romos_ModuleTypeRegistry_h

namespace romos
{

/** A data structure to store information about the various available module types. Conatins also
(pointer to) a factory function "createModule" that is supposed to be used to instantiate modules
of the respective type. */

class ModuleTypeInfo
{
public:
  int id = -1;  // id must be assigned by the type registry object when a type is registered
  Module* (*createModule)() = nullptr; // returns a pointer to a new (subclass of) Module object
  bool hasEditor   = true;  // maybe remove
  //bool hasTreeNode = true;   // should be shown as a node in the tree
  bool hasHeader   = true;   // show a header at the top of the block (most do, so default to true)
  std::string shortName, fullName, description;
  std::vector<std::string> inputShortNames, inputFullNames, inputDescriptions;
  std::vector<std::string> outputShortNames, outputFullNames, outputDescriptions;
  std::string category; // may be a path with "." as delimiter

  /** Used to store information about one of the inputs - the short name is what appears on the 
  pins in the structure view. The long name and description (are supposed to be) used for a more
  verbose information screen. */
  void addInputPinInfo(const char* shortName, const char* fullName = "",
    const char* description = "");

  /** Same as addInputPinInfo but for an output pin. */
  void addOutputPinInfo(const char* shortName, const char* fullName = "",
    const char* description = "");
};

//=================================================================================================

class TopLevelModule;

/** This class is used to store information about the available module types by maintining an array 
of (pointers to) ModuleTypeInfo objects. This structure also contains a function pointer to
a creator/factory function that can be used to instantiate objects of the respective atomic module
type. To access the information and/or invoke the factory functions, there is a global object
moduleFactory of class ModuleFactory in the romos namespace. */

class ModuleFactory
{

public:

  /** Constructor. Pre-populates the registry with standard modules. */
  ModuleFactory();

  /** Destructor. Cleans up the memory. */
  ~ModuleFactory();

  //-----------------------------------------------------------------------------------------------
  // \name Creation/deletion of module instances

  /** Creates a module with given (long) name. */

  /** Creates a module (using the new operator). The kind of the module (i.e. which subclass) will 
  be determined by "fullTypeName". It should be one the types that were previously registered
  via registerModuleType otherwise an error will occur. After creation, it will call the 
  initialize() method of the module and then do some other calls that allocate memory, reset the 
  state and some other stuff that is common to the initialization of all modules (these additional 
  calls avoid code-duplication in the respective initialize() methods of the module subclass). */
  Module* createModule(const std::string& fullTypeName, const std::string& name = "", 
    int x = 0, int y = 0, bool polyphonic = false);

  /** Creates a module with given unique identifier. */
  Module* createModule(int id, const std::string& name = "", int x = 0, int y = 0, 
    bool polyphonic = false);

  /** Creates a top-level module. In Liberty, this is called just once on start-up. */
  TopLevelModule* createTopLevelModule(const std::string& name = "TopLevelModule", 
    int x = 0, int y = 0, bool polyphonic = false) const;

  /** Deletes the passed module by first calling its cleanUp() method and then using the 
  delete-operator. */
  void deleteModule(romos::Module* moduleToDelete);


  // from ModuleFactory:
  //static romos::Module* createModule(int typeIdentifier, rosic::rsString name = rosic::rsString(), 
  //  int x = 0, int y = 0, bool polyphonic = false);

  // todo: copy/edit comments from ModuleFactory

  //-----------------------------------------------------------------------------------------------
  // \name Inquiry

  /** Returns the unquie id for the module with given (long) name. Returns -1, if the module 
  type with given name was never registered. */
  int getModuleId(const std::string& fullTypeName);
  // rename to getModuleTypeId

  /** Returns a pointer the full type info object given the name of the type (if the name is
  unknown, it returns a nullptr). */
  ModuleTypeInfo* getModuleTypeInfo(const std::string& fullTypeName);

  /** Returns the short name (taht appears at the top of the blocks) give the full name of a 
  module type */
  std::string getShortTypeName(const std::string& fullTypeName)
  {
    return getModuleTypeInfo(fullTypeName)->shortName;
  }

  /** Returns a pointer the type info object given the index/id of the type (if the index is
  out of range, it returns a nullptr). */
  ModuleTypeInfo* getModuleTypeInfo(size_t id);

  /** Preliminary - todo: disallow editing for certain simple modules such as adders, unit-delays,
  etc - see old implementation. */
  bool isModuleNameEditable(int id) { return true; } 
  // maybe we should have a nameEditable flag in the Module baseclass that defaults to true

  bool hasModuleTypeEditor(size_t id) { return (*typeInfos)[id]->hasEditor; }

  /** Returns the number of available module types. */
  size_t getNumModuleTypes() { return typeInfos->size(); }

  /** Returns true, iff a module of the given type exists (i.e. was registered). */
  inline bool doesTypeExist(const std::string& fullName)
  {
    return getModuleId(fullName) != -1;
  }

  //-----------------------------------------------------------------------------------------------
  // \name Registration of module types

  /** Registers the module with given type info. You need to pass a pointer that you may create via
  new and forget about it. This object takes over responsibility for deleting it */
  void registerModuleType(ModuleTypeInfo* newTypeInfoToAdd);

  /** Removes the module type with given name (if it is found in the list of registered types, 
  otherwise does nothing). */
  void removeModuleType(const std::string& fullTypeName);

  /** Registers all the standard module types that are commonly used. */
  void registerStandardModules();

  /** Registers a bunch of programatically pre-built containers. */
  void registerPreBuiltContainers();

  /** Cleans up the memory, i.e. the registered type-info objects that were passed as pointers. */
  void clearRegisteredTypes();

protected:

  void ensureTypeInfoArrayAllocated();

  void setupModule(romos::Module* module, const std::string& name, int x, int y, 
    bool polyphonic) const;

  //std::vector<ModuleTypeInfo*> typeInfos;          // causes false memleak error triggers
  //rosic::rsArrayTools<ModuleTypeInfo*> typeInfos;       // dito
  std::vector<ModuleTypeInfo*>* typeInfos = nullptr; // using a pointer-to vector avoids this

};

extern ModuleFactory moduleFactory;  // declaration of the global object















//=================================================================================================

/*
Module-Types to do:

for sure:

regular:

Difference:
y[n] = x[n] - x[n-1]

Clip:
if( in < min )
  out = min;
else if( in > max )
  out = max;
else
  out = in;

CenterClip:
if( in < min || in > max )
  out = in;
else
  out = centerValue;  // could be 0.5*(min+max) but maybe also set from outside

Select:
if( in1 >= threshold )
  out = in2;
else
  out = in3;

  // nah, make a multi-selector with variable number of ins and switch on round(in1)

Min, Max, Floor, Ceil, Round, Abs, Angle(x, y), Length(x, y), ToRA(x, y) (radius/angle), 
ToXY(r, a), Equal, InEqual, Greater, GreaterOrEqual, Less, LessOrEqual, And, Or, Not, Nand, Nor, 
Xor

associated with per-voice values (maybe have a PerVoiceValue Module):

associated with global values (maybe have a GlobalValue Module):
SampleRate, Tempo, Time

detection:
-zero-crossing finder (outputs a trigger signal and an offset to indicate subsample location of the 
 crossing)
-maximum-finder: x[n-1] > x[n] && x[n-1] > x[n-2] - output a trigger and an offset (found by 
 parabolic interpolation)
 -min-finder likewise
-Schmitt-Trigger (see CSound-book, page 568)

maybe:
SendAudio, ReceiveAudio, SendEvent, ReceiveEvent

UpSample, DownSample
-when connecting modules with different oversampling factors, automatic conversion is done by zero 
 insertion or discarding

*/

}

#endif
