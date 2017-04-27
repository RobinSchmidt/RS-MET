#include "rosic_DspScriptInterpreter.h"
#include "rosic_DspModules.h"
using namespace rosic;

//-------------------------------------------------------------------------------------------------
// construction/destruction:

DspScriptInterpreter::DspScriptInterpreter()
{
  sampleRate = 44100.0;

  // add the DSP functions to the list of functions in the inherited ExpressionEvaluator object:
  functionList.Add( new FunctionFactoryLowpass6()             );
  functionList.Add( new FunctionFactorySampleDelay()          );
  functionList.Add( new FunctionFactorySineOscillator()       );

}



//-------------------------------------------------------------------------------------------------
// setup:

void DspScriptInterpreter::setSampleRate(double newSampleRate)
{
  mutex.lock();
  if( newSampleRate > 0.0 )
    sampleRate = newSampleRate;
  ExpressionEvaluator::assignVariable("sampleRate", sampleRate);
  initFunctionList();
  mutex.unlock();
}

void DspScriptInterpreter::initFunctionList()
{
  mutex.lock();

  ExpressionEvaluator::initFunctionList();

  // add the DSP specific functions:
  /*
  functionList.Add( new FunctionFactoryLowpass6()             );
  functionList.Add( new FunctionFactorySampleDelay()          );
  functionList.Add( new FunctionFactorySineOscillator()       );
  */

  mutex.unlock();
}
