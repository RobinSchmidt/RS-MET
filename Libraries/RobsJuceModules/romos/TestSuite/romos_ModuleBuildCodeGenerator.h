#ifndef romos_ModuleBuildCodeGenerator_h
#define romos_ModuleBuildCodeGenerator_h

//#include "../framework/romos_ContainerModule.h"

namespace rsTestRomos
{

  /**  

  Given a Module object (which might be a container with complex structure), this class creates C++ code that can be used to build the 
  module.

  */

  class ModuleBuildCodeGenerator
  {

  public:

    //-------------------------------------------------------------------------------------------------------------------------------------
    // code generation:

    /** Returns the string representing the code that can be used to build the module. */
    static rosic::rsString getCodeForModule(romos::Module *module);

    //=====================================================================================================================================

  protected:

    /** creates a suitable variable-name to be used inside the code to refer to the module. For most modules, this will be just the module's
    name with spaces removed and first letter converted to lowercase. But some atomic modules require special treatment, because their names
    are mere arithmetic operators ("+", "-", "*", etc.). */
    static rosic::rsString makeModuleVariableName(romos::Module *module);

    /** Returns the number of modules of the same type that come before the passed module within the parent. */
    static int getNumOfSameModulesBefore(romos::Module *module);

    /** Creates the variable names to be used for the child-modules in the generated code. */
    static rosic::rsDynamicArray<rosic::rsString> createVariableNames(romos::ContainerModule *container);

    /** Creates the module names to be used for the child-modules in the generated code, including a comma and possibly padding whitespace
    to be used in the function call. */
    static rosic::rsDynamicArray<rosic::rsString> createModuleNames(romos::ContainerModule *container);

    /** Creates the strings that are used to infer the type-indentifiers for the modules in the generated code. */
    static rosic::rsDynamicArray<rosic::rsString> createTypeRetrievalStrings(romos::ContainerModule *container);

    /** Finds the longest string in the passed array and pads all shorter strings with whitespace such that all strings in the array have 
    the same length. */
    static void padShortStringsWithSpace(rosic::rsDynamicArray<rosic::rsString> &stringsToPad);

    /** Given an array with names that are associated with child-modules of the passed ContainerModule (variable-names to appear in the 
    generated code, for example), this function returns the maximum lengths of the names for atomic and container modules separately. This
    information can be used to create proper whitespace padding  for code alignment. */
    static void getMaxNameLengths(rosic::rsDynamicArray<rosic::rsString> &names, romos::ContainerModule *container, 
      int &maxAtomicNameLength, int &maxContainerNameLength);

    /** Creates a string containing a stretch pf whitespaces that can be used to append (or prepend) to the passed string to make it
    reach the desired length. */
    static rosic::rsString createPadding(rosic::rsString stringToPad, int desiredLength);

    /** Creates and returns either the string "true " or the string "false", switching on the input parameter. */
    static rosic::rsString createTrueFalseString(bool trueOrFalse);

    /** Creates and returns either the string "true,  " or the string "false, ", switching on the input parameter. */
    static rosic::rsString createTrueFalseStringWithCommaAndSpace(bool trueOrFalse);

  };

} 

#endif
