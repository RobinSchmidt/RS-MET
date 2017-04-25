#include "rosic_DspWorkbench.h"
using namespace rosic;

//-------------------------------------------------------------------------------------------------
// construction/destruction:

DspWorkbench::DspWorkbench()
{
  // init member variables and embedded objects:
  sampleRate         = 44100.0;
  oversamplingFactor = 1;
  expressionIsValid  = false;
  upsamplerL.setSubDivision(oversamplingFactor);
  upsamplerR.setSubDivision(oversamplingFactor);
  antiAliasFilterL.setSubDivision(oversamplingFactor);
  antiAliasFilterR.setSubDivision(oversamplingFactor);

  algorithm = "outL=inL;\noutR=inR;";
  setAlgorithmString("outL=inL;\noutR=inR;");

  // create and initialize variables for input and output and retrieve the addresses where they 
  // are stored:
  scriptInterpreter.assignVariable("inL",   0.0);  
  scriptInterpreter.assignVariable("inR",   0.0);
  scriptInterpreter.assignVariable("outL",  0.0);  
  scriptInterpreter.assignVariable("outR",  0.0);
  inL  = scriptInterpreter.getVariableAddress("inL");
  inR  = scriptInterpreter.getVariableAddress("inR");
  outL = scriptInterpreter.getVariableAddress("outL");
  outR = scriptInterpreter.getVariableAddress("outR");

  initializeArrays();

  //markPresetAsClean();
}

DspWorkbench::~DspWorkbench()
{

}

//-------------------------------------------------------------------------------------------------
// parameter settings:

void DspWorkbench::setSampleRate(double newSampleRate)
{
  if( newSampleRate > 0.0 )
    sampleRate = newSampleRate;
  scriptInterpreter.setSampleRate(oversamplingFactor*sampleRate);
  //evaluator.assignVariable("sampleRate", oversamplingFactor*sampleRate);  
}

void DspWorkbench::setOversamplingFactor(int newOversamplingFactor)
{
  oversamplingFactor = max(1, newOversamplingFactor);
  if( oversamplingFactor > 1 )
  {
    upsamplerL.setSubDivision(oversamplingFactor);
    upsamplerR.setSubDivision(oversamplingFactor);
    antiAliasFilterL.setSubDivision(oversamplingFactor);
    antiAliasFilterR.setSubDivision(oversamplingFactor);
  }
  scriptInterpreter.setSampleRate(oversamplingFactor*sampleRate);
  //evaluator.assignVariable("sampleRate", oversamplingFactor*sampleRate);  
  //markPresetAsDirty();
}

bool DspWorkbench::setAlgorithmString(char *newAlgorithmString)
{

  // pass the new string to the ExpressionEvaluator object and return the boolean result of this 
  // set-function to the calling function:
  //evaluatorLock.lock();
  algorithm = newAlgorithmString;
  expressionIsValid = scriptInterpreter.setExpressionString(newAlgorithmString);      
  //evaluatorLock.unlock();

  //markPresetAsDirty();
  return expressionIsValid;
}

void DspWorkbench::setParameter(int index, double newValue)
{
  if( index >= 0 && index < numParameters )
  {
    parameters[index] = newValue;
    scriptInterpreter.assignVariable(getInternalParameterName(index),  parameters[index]);
    scriptInterpreter.assignVariable(getParameterName(index),  parameters[index]);
  }
  //markPresetAsDirty();
}

void DspWorkbench::setParameterName(int index, char* newName)
{
  if( index < 0 || index >= numParameters )
    return;

  // free old and allocate new memory for the name:
  if( parameterClearTextNames[index] != NULL )
  {
    delete[] parameterClearTextNames[index];
    parameterClearTextNames[index] = NULL;
  }
  if( newName != NULL )
  {
    int newLength         = strlen(newName);
    parameterClearTextNames[index] = new char[newLength+1];
    for(int c=0; c<=newLength; c++) // the <= is valid here, because we have one more cell allocated
      parameterClearTextNames[index][c] = newName[c];
  }
}

//-------------------------------------------------------------------------------------------------
// inquiry:

int DspWorkbench::getOversamplingFactor() const
{
  return oversamplingFactor;
}

const char* DspWorkbench::getAlgorithmString()
{
  return scriptInterpreter.getExpressionString();
}


double DspWorkbench::getParameter(int index) const
{
  if( index >= 0 && index < numParameters )
    return parameters[index];
  else
  {
    DEBUG_BREAK;
    return 0.0;
  }
}

const char* DspWorkbench::getParameterName(int index) const
{
  if( index >= 0 && index < numParameters )
    return parameterClearTextNames[index];
  else
  {
    DEBUG_BREAK;
    return NULL;
  }
}

const char* DspWorkbench::getInternalParameterName(int index) const
{
  if( index >= 0 && index < numParameters )
    return parameterInternalNames[index];
  else
  {
    DEBUG_BREAK;
    return NULL;
  }
}


//-------------------------------------------------------------------------------------------------

void DspWorkbench::initializeArrays()
{
  char charBuffer[8];

  for(int p=0; p<numParameters; p++)
  {
    parameters[p]        = 0.5;

    sprintf(charBuffer, "%s%d", "Par", p);

    parameterClearTextNames[p] = NULL;
    parameterInternalNames[p]  = NULL;

    setParameterInternalName(p, charBuffer);
    setParameterName(        p, charBuffer);
    setParameter(p, 0.5); 
      // this will assign the variable name to the numeric value in the expression evaluator
  }
}

void DspWorkbench::setParameterInternalName(int index, char* newName)
{
  if( index < 0 || index >= numParameters )
    return;

  // free old and allocate new memory for the name:
  if( parameterInternalNames[index] != NULL )
  {
    delete[] parameterInternalNames[index];
    parameterInternalNames[index] = NULL;
  }
  if( newName != NULL )
  {
    int newLength         = strlen(newName);
    parameterInternalNames[index] = new char[newLength+1];
    for(int c=0; c<=newLength; c++) // the <= is valid here, because we have one more cell allocated
      parameterInternalNames[index][c] = newName[c];
  }
}





















