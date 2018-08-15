#ifndef romos_ModuleFactory_h
#define romos_ModuleFactory_h

namespace romos
{

#ifdef RS_BUILD_OLD_MODULE_FACTORY

/** This class is used to handle the creation and destruction of modules. We use a factory in 
order to avoid duplication of code which otherwise would be necessary in the constructors of all
the individual modules. With a factory, we can write the common code once and delegate the 
subclass-specific things into an initialize method. A baseclass constructor won't cut it in this
case because a baseclass constructor will always be executed before the subclass constructor - and
some common things must be executed after the non-common things. Also, each module subclass must 
allocate memory (via a call to allocateMemory) - but we can't call allocateMemory in the baseclass
constructor because allocateMemory has to be virtual (each module class must allocate a different 
amount of memory). Even if it's not overriden and the baseclass implementation is used, it must be
called after the subclass has initialized the numInputs/numOutputs members. All these 
considerations justify to use a factory for module creation.

The deletion through a factory is done because...tbc */

class ModuleFactory
{

public:

  //-----------------------------------------------------------------------------------------------
  // module creation and deletion:

  /** Creates a module (using the new operator). The kind of the module (i.e. which subclass) will 
  be determined by "typeIdentifier". After creation, it will call the initialize() method of the 
  module and then do some other calls that allocate memory, reset the state and some other stuff 
  that is common to the initialization of all modules (these additional calls avoid 
  code-duplication in the respective initialize() methods of the module subclass). */
  static romos::Module* createModule(int typeIdentifier, rosic::rsString name = rosic::rsString(), 
    int x = 0, int y = 0, bool polyphonic = false);

  /** Deletes the passed module by first calling its cleanUp() method and then using the 
  delete-operator. */
  static void deleteModule(romos::Module *moduleToDelete);

};

// todo: have a datastructure ModuleTypeInfo containing
// a unique integer ID, a name (string), a creator function and (maybe, if necessary, a deletor
// function.
// to register modules, we would the just do
// moduleFactory.registerModuleType(PHASOR, "Phasor", createPhasor)
// maybe instead of createPhasor, &(new Phasor) could work? ..try it - would avoid the boilerplate
// function createPhasor which would otherwise have to be defined
// classes ModuleTypeRegistry and ModuleFactory can be merged into one single class (ModuleFactory)

#endif // RS_BUILD_OLD_MODULE_FACTORY

}

#endif
