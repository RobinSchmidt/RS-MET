//#include "rosic_Module.h"
//using namespace rosic;

char* Module::getModulatableParameterName(int index)
{
  if( index >= 0 && index < modulatableParameters.getNumElements() )
    return modulatableParameters[index].getName();
  else
  {
    DEBUG_BREAK;  // parameter-index out of range
    return NULL;
  }
}






