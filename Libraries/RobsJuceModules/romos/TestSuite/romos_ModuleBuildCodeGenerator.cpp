#include "romos_ModuleBuildCodeGenerator.h"
using namespace romos;
 
//-----------------------------------------------------------------------------------------------------------------------------------------
// code generation:

#define S(x) rosic::rsString(x)

rosic::rsString ModuleBuildCodeGenerator::getCodeForModule(romos::Module *module)
{
  // catch some special conditions:
  if( module == NULL )
    return rosic::rsString();
  if( !module->isContainerModule() )
    return rosic::rsString("Atomic module - nothing to do");

  // declare some local variables fo later use:
  romos::ContainerModule *container = dynamic_cast<romos::ContainerModule*> (module);
  romos::Module *child;
  rosic::rsString code;
  rosic::rsString indent = S("  ");
  rosic::rsString tmp, padding, padding2;

  // get variable names maximum name lengths (for alingment):
  rosic::rsDynamicArray<rosic::rsString> variableNames  = createVariableNames(container);
  rosic::rsDynamicArray<rosic::rsString> moduleNames    = createModuleNames(container);
  rosic::rsDynamicArray<rosic::rsString> typeRetrievals = createTypeRetrievalStrings(container);
  int maxAtomicVariableNameLength, maxAtomicModuleNameLength, maxContainerVariableNameLength, maxContainerModuleNameLength, 
    maxVariableNameLength;
  getMaxNameLengths(variableNames, container, maxAtomicVariableNameLength, maxContainerVariableNameLength);
  getMaxNameLengths(moduleNames,   container, maxAtomicModuleNameLength,   maxContainerModuleNameLength);
  maxVariableNameLength = RAPT::rsMax(maxAtomicVariableNameLength, maxContainerVariableNameLength);

  // create the function entry code:
  code += S("romos::Module* TestModuleBuilder::create") + S(module->getName()) + 
          S("(const rosic::rsString &name, int x, int y, bool polyphonic)\n");
  code += S("{\n");
  code += indent + S("ContainerModule *module = (ContainerModule*) ")
                 + S("ModuleFactory::createModule(ModuleTypeRegistry::CONTAINER, name, x, y, polyphonic);\n");
  code += indent + S("\n");

  // 1st pass through the child-modules: create code for creation of atomic child modules:
  for(unsigned int i = 0; i < container->getNumChildModules(); i++)
  {
    child = container->getChildModule(i);
    if( child->isContainerModule() )
      continue;
    //tmp = ModuleTypeRegistry::getSoleInstance()->getModuleTypeStringFromIdentifier(child->getTypeIdentifierOld());
    tmp = child->getTypeName(); // new - correct?

    code    += indent;
    padding  = createPadding(variableNames[i], maxAtomicVariableNameLength);
    code    += S("romos::Module *") + variableNames[i] + padding + S(" = module->addChildModule(") + typeRetrievals[i] + S(" ");
    padding  = createPadding(moduleNames[i], maxAtomicModuleNameLength);
    code    += moduleNames[i] + rosic::rsString(", ") + padding;
    code    += rosic::rsString::fromIntWithLeadingSpaces(child->getPositionX(), 2) + S(", "); 
    code    += rosic::rsString::fromIntWithLeadingSpaces(child->getPositionY(), 2) + S(", ");  
    code    += createTrueFalseStringWithCommaAndSpace(child->isPolyphonic());
    code    += rosic::rsString("false);");
    code    += S("\n");
  }    
  code += S("\n");

  // 2nd pass through the child-modules: create code for creation of non-atomic (container) child modules:
  bool hasChildContainers = false;
  for(unsigned int i=0; i<container->getNumChildModules(); i++)
  {
    child = container->getChildModule(i);
    if( child->isContainerModule() )
    {
      // call this function recursively and prepend the result: 
      rosic::rsString childCreationCode = getCodeForModule(child);
      code = childCreationCode + code;

      // generate code that invokes the just prepended child-creation function:
      code += indent;
      rosic::rsString name = variableNames[i].toUpperCase(0, 0);
      padding = createPadding(name, maxContainerVariableNameLength);
      code += S("romos::Module *") + variableNames[i] + padding + S(" = module->addChildModule(create") + name + S("(") + padding;
      code += S("\"") + name + S("\"") + S(", ") + padding;
      code += rosic::rsString::fromIntWithLeadingSpaces(child->getPositionX(), 2) + S(", "); 
      code += rosic::rsString::fromIntWithLeadingSpaces(child->getPositionY(), 2) + S(", ");
      code += createTrueFalseString(child->isPolyphonic()) + S("));");  ;
      code += S("\n");
      hasChildContainers = true;
    }
  }    

  if( hasChildContainers == true )
    code += S("\n");
  code += indent + S("module->sortChildModuleArray();\n");
  code += S("\n");

  // 3rd pass through the child-modules: create code for creation of connections:
  for(unsigned int i = 0; i < container->getNumChildModules(); i++)
  {
    //DEBUG_BREAK;

    child = container->getChildModule(i);
    for(unsigned int pinIndex = 0; pinIndex < child->getNumInputPins(); pinIndex++)
    {
      Module *sourceModule = child->inputPins[pinIndex].sourceModule;
      if( sourceModule != NULL )
      {
        // when the source is an output module, we actually want the hosting container:
        if( sourceModule->isOutputModule() )  
          sourceModule = sourceModule->getParentModule();

        int sourceModuleIndex = container->getIndexOfChildModule(sourceModule);
        int targetModuleIndex = i;
        int sourceOutPinIndex = child->inputPins[pinIndex].outputIndex;
        int targetInPinIndex  = pinIndex;
     
        padding  = createPadding(variableNames[sourceModuleIndex], maxVariableNameLength);
        padding2 = createPadding(variableNames[targetModuleIndex], maxVariableNameLength);

        code += indent;
        code += S("module->addAudioConnection(") + variableNames[sourceModuleIndex] + S(", ") + padding + S(sourceOutPinIndex);
        code += S(", ") + variableNames[targetModuleIndex] + S(", ") + padding2 + S(targetInPinIndex) + S(");");
        code += S("\n");
      }
    }
  }

  // create function-leave code and return the result:
  code += S("\n");
  code += indent + S("return module;\n");
  code += S("}\n");
  code += "\n\n";
  return code;
}

rosic::rsString ModuleBuildCodeGenerator::makeModuleVariableName(romos::Module *module)
{
  rosic::rsString name;
  if( module->isContainerModule() )
    name = module->getName();
  else
  {
    //name = module->getTypeNameOld();  
    name = module->getTypeName();  
    int suffix = getNumOfSameModulesBefore(module) + 1;
    name += suffix;
  }
  name = name.toLowerCase(0, 0);
  return name;
}

int ModuleBuildCodeGenerator::getNumOfSameModulesBefore(romos::Module *module)
{
  romos::ContainerModule *parent = module->getParentModule();
  if( parent == NULL )
    return 0;
  else
  {
    int count  = 0;
    for(unsigned int i=0; i<parent->getNumChildModules(); i++)
    {
      if( parent->getChildModule(i)->getTypeId() == module->getTypeId() )
      {
        if( parent->getChildModule(i) == module )
          return count;
        else
          count++;
      }
    }
    return count;
  }
}

rosic::rsDynamicArray<rosic::rsString> ModuleBuildCodeGenerator::createVariableNames(romos::ContainerModule *container)
{
  rosic::rsDynamicArray<rosic::rsString> varNames;
  for(unsigned int i = 0; i < container->getNumChildModules(); i++)
    varNames.appendElement(makeModuleVariableName(container->getChildModule(i)));
  return varNames;
}

rosic::rsDynamicArray<rosic::rsString> ModuleBuildCodeGenerator::createModuleNames(romos::ContainerModule *container)
{
  rosic::rsDynamicArray<rosic::rsString> names;
  rosic::rsString name;
  for(unsigned int i = 0; i < container->getNumChildModules(); i++)
  {
    name = rosic::rsString("\"") + container->getChildModule(i)->getName() + rosic::rsString("\"");
    names.appendElement(name);
  }
  return names;
}
   
rosic::rsDynamicArray<rosic::rsString> ModuleBuildCodeGenerator::createTypeRetrievalStrings(romos::ContainerModule *container)
{
  rosic::rsDynamicArray<rosic::rsString> typeRetrievalStrings;
  for(unsigned int i = 0; i < container->getNumChildModules(); i++)
  {
    romos::Module *child = container->getChildModule(i);

    //rosic::rsString tmp = ModuleTypeRegistry::getSoleInstance()->getModuleTypeStringFromIdentifier(child->getTypeIdentifierOld());
    rosic::rsString tmp = child->getTypeName();

    tmp = S("getTypeId(\"") + S(tmp) + S("\"),"); // old - won't work
    typeRetrievalStrings.appendElement(tmp);
  }
  padShortStringsWithSpace(typeRetrievalStrings); 
  return typeRetrievalStrings;
}

void ModuleBuildCodeGenerator::padShortStringsWithSpace(rosic::rsDynamicArray<rosic::rsString> &stringsToPad)
{
  int maxLength = 0;

  // first loop could also be factored out into a function getLengthOfLongestString
  for(int i=0; i<stringsToPad.getNumElements(); i++)
  {
    if( stringsToPad[i].getLength() > maxLength )
      maxLength = stringsToPad[i].getLength();
  }
  for(int i=0; i<stringsToPad.getNumElements(); i++)
  {
    if( stringsToPad[i].getLength() < maxLength )
      stringsToPad[i] = stringsToPad[i] + rosic::rsString::createWhiteSpace(maxLength-stringsToPad[i].getLength());
  }
}
  
void ModuleBuildCodeGenerator::getMaxNameLengths(rosic::rsDynamicArray<rosic::rsString> &names, romos::ContainerModule *container, 
                       int &maxAtomicNameLength, int &maxContainerNameLength)
{
  maxAtomicNameLength    = 0;
  maxContainerNameLength = 0;
  for(unsigned int i = 0; i < container->getNumChildModules(); i++)
  {
    romos::Module *child = container->getChildModule(i);
    if( child->isContainerModule() && names[i].getLength() > maxContainerNameLength )
      maxContainerNameLength = names[i].getLength();
    else if( names[i].getLength() > maxAtomicNameLength )
      maxAtomicNameLength = names[i].getLength();
  }
}

rosic::rsString ModuleBuildCodeGenerator::createPadding(rosic::rsString stringToPad, int desiredLength)
{
  rassert( stringToPad.getLength() <= desiredLength );
  int missingLength = RAPT::rsMax(desiredLength-stringToPad.getLength(), 0);
  return rosic::rsString::createWhiteSpace(missingLength);
}

rosic::rsString ModuleBuildCodeGenerator::createTrueFalseString(bool trueOrFalse)
{
  if( trueOrFalse == true )
    return rosic::rsString("true ");
  else
    return rosic::rsString("false");
}

rosic::rsString ModuleBuildCodeGenerator::createTrueFalseStringWithCommaAndSpace(bool trueOrFalse)
{
  if( trueOrFalse == true )
    return rosic::rsString("true,  ");
  else
    return rosic::rsString("false, ");
}
