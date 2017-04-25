#include "rosic_AngelScriptInterpreter.h"
using namespace rosic;

//-------------------------------------------------------------------------------------------------
// construction/destruction:

AngelScriptInterpreter::AngelScriptInterpreter()
{
  sampleRate = 44100.0;
}



//-------------------------------------------------------------------------------------------------
// setup:

void AngelScriptInterpreter::setSampleRate(double newSampleRate)
{
  if( newSampleRate > 0.0 )
    sampleRate = newSampleRate;
}

void AngelScriptInterpreter::initFunctionList()
{

}