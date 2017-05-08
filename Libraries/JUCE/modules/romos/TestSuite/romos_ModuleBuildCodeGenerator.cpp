#include "romos_ModuleBuildCodeGenerator.h"
using namespace romos;
 
//-----------------------------------------------------------------------------------------------------------------------------------------
// code generation:

#define S(x) rosic::String(x)

rosic::String ModuleBuildCodeGenerator::getCodeForModule(romos::Module *module)
{
  // catch some special conditions:
  if( module == NULL )
    return rosic::String();
  if( !module->isContainerModule() )
    return rosic::String("Atomic module - nothing to do");

  // declare some local variables fo later use:
  romos::ModuleContainer *container = dynamic_cast<romos::ModuleContainer*> (module);
  romos::Module *child;
  rosic::String code;
  rosic::String indent = S("  ");
  rosic::String tmp, padding, padding2;

  // get variable names maximum name lengths (for alingment):
  rosic::Array<rosic::String> variableNames  = createVariableNames(container);
  rosic::Array<rosic::String> moduleNames    = createModuleNames(container);
  rosic::Array<rosic::String> typeRetrievals = createTypeRetrievalStrings(container);
  int maxAtomicVariableNameLength, maxAtomicModuleNameLength, maxContainerVariableNameLength, maxContainerModuleNameLength, 
    maxVariableNameLength;
  getMaxNameLengths(variableNames, container, maxAtomicVariableNameLength, maxContainerVariableNameLength);
  getMaxNameLengths(moduleNames,   container, maxAtomicModuleNameLength,   maxContainerModuleNameLength);
  maxVariableNameLength = rmax(maxAtomicVariableNameLength, maxContainerVariableNameLength);

  // create the function entry code:
  code += S("romos::Module* TestModuleBuilder::create") + S(module->getName()) + 
          S("(const rosic::String &name, int x, int y, bool polyphonic)\n");
  code += S("{\n");
  code += indent + S("ModuleContainer *module = (ModuleContainer*) ")
                 + S("ModuleFactory::createModule(ModuleTypeRegistry::CONTAINER, name, x, y, polyphonic);\n");
  code += indent + S("\n");

  // 1st pass through the child-modules: create code for creation of atomic child modules:
  for(unsigned int i = 0; i < container->getNumChildModules(); i++)
  {
    child = container->getChildModule(i);
    if( child->isContainerModule() )
      continue;
    tmp = ModuleTypeRegistry::getSoleInstance()->getModuleTypeStringFromIdentifier(child->getTypeIdentifier());
    code    += indent;
    padding  = createPadding(variableNames[i], maxAtomicVariableNameLength);
    code    += S("romos::Module *") + variableNames[i] + padding + S(" = module->addChildModule(") + typeRetrievals[i] + S(" ");
    padding  = createPadding(moduleNames[i], maxAtomicModuleNameLength);
    code    += moduleNames[i] + rosic::String(", ") + padding;
    code    += rosic::String::fromIntWithLeadingSpaces(child->getPositionX(), 2) + S(", "); 
    code    += rosic::String::fromIntWithLeadingSpaces(child->getPositionY(), 2) + S(", ");  
    code    += createTrueFalseStringWithCommaAndSpace(child->isPolyphonic());
    code    += rosic::String("false);");
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
      rosic::String childCreationCode = getCodeForModule(child);
      code = childCreationCode + code;

      // generate code that invokes the just prepended child-creation function:
      code += indent;
      rosic::String name = variableNames[i].toUpperCase(0, 0);
      padding = createPadding(name, maxContainerVariableNameLength);
      code += S("romos::Module *") + variableNames[i] + padding + S(" = module->addChildModule(create") + name + S("(") + padding;
      code += S("\"") + name + S("\"") + S(", ") + padding;
      code += rosic::String::fromIntWithLeadingSpaces(child->getPositionX(), 2) + S(", "); 
      code += rosic::String::fromIntWithLeadingSpaces(child->getPositionY(), 2) + S(", ");
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

rosic::String ModuleBuildCodeGenerator::makeModuleVariableName(romos::Module *module)
{
  rosic::String name;
  if( module->isContainerModule() )
    name = module->getName();
  else
  {
    name = module->getTypeName();  
    int suffix = getNumOfSameModulesBefore(module) + 1;
    name += suffix;
  }
  name = name.toLowerCase(0, 0);
  return name;
}

int ModuleBuildCodeGenerator::getNumOfSameModulesBefore(romos::Module *module)
{
  romos::ModuleContainer *parent = module->getParentModule();
  if( parent == NULL )
    return 0;
  else
  {
    int count  = 0;
    for(unsigned int i=0; i<parent->getNumChildModules(); i++)
    {
      if( parent->getChildModule(i)->getTypeIdentifier() == module->getTypeIdentifier() )
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

rosic::Array<rosic::String> ModuleBuildCodeGenerator::createVariableNames(romos::ModuleContainer *container)
{
  rosic::Array<rosic::String> varNames;
  for(unsigned int i = 0; i < container->getNumChildModules(); i++)
    varNames.appendElement(makeModuleVariableName(container->getChildModule(i)));
  return varNames;
}

rosic::Array<rosic::String> ModuleBuildCodeGenerator::createModuleNames(romos::ModuleContainer *container)
{
  rosic::Array<rosic::String> names;
  rosic::String name;
  for(unsigned int i = 0; i < container->getNumChildModules(); i++)
  {
    name = rosic::String("\"") + container->getChildModule(i)->getName() + rosic::String("\"");
    names.appendElement(name);
  }
  return names;
}
   
rosic::Array<rosic::String> ModuleBuildCodeGenerator::createTypeRetrievalStrings(romos::ModuleContainer *container)
{
  rosic::Array<rosic::String> typeRetrievalStrings;
  for(unsigned int i = 0; i < container->getNumChildModules(); i++)
  {
    romos::Module *child = container->getChildModule(i);
    rosic::String tmp = ModuleTypeRegistry::getSoleInstance()->getModuleTypeStringFromIdentifier(child->getTypeIdentifier());
    tmp = S("getTypeId(\"") + S(tmp) + S("\"),");
    typeRetrievalStrings.appendElement(tmp);
  }
  padShortStringsWithSpace(typeRetrievalStrings); 
  return typeRetrievalStrings;
}

void ModuleBuildCodeGenerator::padShortStringsWithSpace(rosic::Array<rosic::String> &stringsToPad)
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
      stringsToPad[i] = stringsToPad[i] + rosic::String::createWhiteSpace(maxLength-stringsToPad[i].getLength());
  }
}
  
void ModuleBuildCodeGenerator::getMaxNameLengths(rosic::Array<rosic::String> &names, romos::ModuleContainer *container, 
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

rosic::String ModuleBuildCodeGenerator::createPadding(rosic::String stringToPad, int desiredLength)
{
  rassert( stringToPad.getLength() <= desiredLength );
  int missingLength = rmax(desiredLength-stringToPad.getLength(), 0);
  return rosic::String::createWhiteSpace(missingLength);
}

rosic::String ModuleBuildCodeGenerator::createTrueFalseString(bool trueOrFalse)
{
  if( trueOrFalse == true )
    return rosic::String("true ");
  else
    return rosic::String("false");
}

rosic::String ModuleBuildCodeGenerator::createTrueFalseStringWithCommaAndSpace(bool trueOrFalse)
{
  if( trueOrFalse == true )
    return rosic::String("true,  ");
  else
    return rosic::String("false, ");
}
