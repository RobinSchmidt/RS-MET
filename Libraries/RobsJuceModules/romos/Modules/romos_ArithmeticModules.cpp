//#include "romos_ArithmeticModules.h"

namespace romos // move this namespace wrapping into romos.cpp
{

void ConstantModule::initialize()
{
  initOutputPins({ "" });
  hasHeaderFlag  = false;
  value          = 0.0;
  audioOutputs   = &value;
  outFrameStride = 0;
}
void ConstantModule::processMonoFrame(Module* /*module*/, int /*voiceIndex*/)
{
  // nothing to do - the constant sits just there in the output pin
}
void ConstantModule::processPolyFrame(Module *module, int voiceIndex)                {  }
void ConstantModule::processMonoBlock(Module *module, int voiceIndex, int blockSize) {  }
void ConstantModule::processPolyBlock(Module *module, int voiceIndex, int blockSize) {  }
void ConstantModule::clearVoiceBuffer(int voiceIndex)
{
  // do nothing - overriden to avoid resetting the output to zero - the constant just remains in 
  // the output pin
}
void ConstantModule::allocateMemory()
{
  // do nothing - our "audioOutputs" member just points to a constant value and our outFrameStride 
  // is 0
}
void ConstantModule::freeMemory()
{
  // nothing to do here as well
}


// move to rosic or rapt
std::string rsDoubleToString(double value)
{
  char tmpString[64];
  sprintf(tmpString, "%lg", value);
  return std::string(tmpString);
}

void ConstantModule::setModuleName(const std::string& newName)
{
  // old - using rsString:
  //value = newName.asDouble();
  ////Module::setModuleName(rosic::rsString::fromDouble(value));
  //Module::setModuleName(rosic::rsString(value));

  try {
    value = std::stod(newName);
    //value = atof(newName.c_str());
  }
  catch(const std::invalid_argument&) {
    value = 0.0;
  }
  Module::setModuleName(rsDoubleToString(value));
  //Module::setModuleName(std::to_string(value));  // produces trailing zeros
}
CREATE_ASSIGN_FUNCTION_POINTERS(ConstantModule);

//-------------------------------------------------------------------------------------------------

void IdentityModule::initialize()
{
  initInputPins({ "" });
  initOutputPins({ "" });
  hasHeaderFlag = false;
}
INLINE void IdentityModule::process(Module *module, double *ins, double *outs, int voiceIndex)
{
  IdentityModule *identity = (IdentityModule *)module;

  outs[0] = ins[0];


  int dummy = 0;
}
CREATE_AND_ASSIGN_PROCESSING_FUNCTIONS_1(IdentityModule);

//-------------------------------------------------------------------------------------------------

void UnaryMinusModule::initialize()
{
  initInputPins({ "" });
  initOutputPins({ "" });
  hasHeaderFlag = false;
}
INLINE void UnaryMinusModule::process(Module *module, double *ins, double *outs, int voiceIndex)
{
  outs[0] = -ins[0];
}
CREATE_AND_ASSIGN_PROCESSING_FUNCTIONS_1(UnaryMinusModule);

//-------------------------------------------------------------------------------------------------

void ScalerModule::initialize()
{
  initInputPins({ "" });
  initOutputPins({ "" });
  addParameter("Multiplier", "1");
  hasHeaderFlag = false;
}
INLINE void ScalerModule::process(Module *module, double *ins, double *outs, int voiceIndex)
{
  //double p = 1; // preliminary - todo: connect to our parameter

  ScalerModule *sclr = static_cast<ScalerModule*> (module);
  outs[0] = sclr->multiplier * ins[0];
}
void ScalerModule::parameterChanged(int index)
{
  multiplier = parameters[0].value.asDouble(); // maybe rename to toDouble
}
CREATE_AND_ASSIGN_PROCESSING_FUNCTIONS_1(ScalerModule);

//-------------------------------------------------------------------------------------------------

void ReciprocalModule::initialize()
{
  initInputPins({ "" });
  initOutputPins({ "" });
  hasHeaderFlag = false;
}
INLINE void ReciprocalModule::process(Module *module, double *in, double *out, int voiceIndex)
{
  if(*in != 0.0)
    *out = 1.0 / *in;
  else
    *out = 0.0;
}
CREATE_AND_ASSIGN_PROCESSING_FUNCTIONS_1(ReciprocalModule);

//-------------------------------------------------------------------------------------------------

void AdderModule::initialize()
{
  initInputPins({ "", "" });
  initOutputPins({ "" });
  hasHeaderFlag = false;
}
INLINE void AdderModule::process(Module *module, double *in1, double *in2, double *out, 
  int voiceIndex)
{
  *out = *in1 + *in2;
}
CREATE_AND_ASSIGN_PROCESSING_FUNCTIONS_2(AdderModule);

//-------------------------------------------------------------------------------------------------

void SubtractorModule::initialize()
{
  initInputPins({ "", "" });
  initOutputPins({ "" });
  hasHeaderFlag = false;
}
INLINE void SubtractorModule::process(Module *module, double *in1, double *in2, double *out, 
  int voiceIndex)
{
  *out = *in1 - *in2;
}
CREATE_AND_ASSIGN_PROCESSING_FUNCTIONS_2(SubtractorModule);

//-------------------------------------------------------------------------------------------------

void MultiplierModule::initialize()
{
  initInputPins({ "", "" });
  initOutputPins({ "" });
  hasHeaderFlag = false;
}
INLINE void MultiplierModule::process(Module *module, double *in1, double *in2, double *out,
  int voiceIndex)
{
  *out = *in1 * *in2;
}
CREATE_AND_ASSIGN_PROCESSING_FUNCTIONS_2(MultiplierModule);

//-------------------------------------------------------------------------------------------------

void DividerModule::initialize()
{
  initInputPins({ "", "" });
  initOutputPins({ "" });
  hasHeaderFlag = false;
}
INLINE void DividerModule::process(Module *module, double *in1, double *in2, double *out, 
  int voiceIndex)
{
  if(*in2 != 0.0)
    *out = *in1 / *in2;
  else
    *out = 0.0;
}
CREATE_AND_ASSIGN_PROCESSING_FUNCTIONS_2(DividerModule);

//-------------------------------------------------------------------------------------------------

void Adder3Module::initialize()
{
  initInputPins({ "", "", "" });
  initOutputPins({ "" });
  hasHeaderFlag = false;
}
INLINE void Adder3Module::process(Module *module, double *in1, double *in2, double *in3, 
  double *out, int voiceIndex)
{
  *out = *in1 + *in2 + *in3;
}
CREATE_AND_ASSIGN_PROCESSING_FUNCTIONS_3(Adder3Module);

//-------------------------------------------------------------------------------------------------

void Adder4Module::initialize()
{
  initInputPins({ "", "", "", "" });
  initOutputPins({ "" });
  hasHeaderFlag = false;
}
INLINE void Adder4Module::process(Module *module, double *in1, double *in2, double *in3, 
  double *in4, double *out, int voiceIndex)
{
  *out = *in1 + *in2 + *in3 + *in4;  // perhaps optimizable via parentheses
}
CREATE_AND_ASSIGN_PROCESSING_FUNCTIONS_4(Adder4Module);

//-------------------------------------------------------------------------------------------------

void Adder5Module::initialize()
{
  initInputPins({ "", "", "", "", "" });
  initOutputPins({ "" });
  hasHeaderFlag = false;
}
INLINE void Adder5Module::process(Module *module, double *in1, double *in2, double *in3, 
  double *in4, double *in5, double *out, int voiceIndex)
{
  *out = *in1 + *in2 + *in3 + *in4 + *in5; // perhaps optimizable via parentheses
}
CREATE_AND_ASSIGN_PROCESSING_FUNCTIONS_5(Adder5Module);

//-------------------------------------------------------------------------------------------------

void AdderNModule::initialize()
{
  initInputPins({ "", "" });
  initOutputPins({ "" });
}
INLINE void AdderNModule::process(Module *module, double *ins, double *outs, int voiceIndex)
{
  AdderNModule *sumModule = static_cast<AdderNModule*> (module);
  outs[0] = 0.0;
  for(unsigned int i = 0; i < sumModule->getNumInputPins()-1; i++)
    outs[0] += ins[i];
}
void AdderNModule::connectInputPinTo(int inputPinIndex, Module *sourceModule, 
  int sourceOutputPinIndex)
{
  if(!isInputPinConnected(inputPinIndex))
  {
    Module::connectInputPinTo(inputPinIndex, sourceModule, sourceOutputPinIndex);
    if(inputPinIndex == getNumInputPins()-1)
      addAudioInput("");
  }
}
void AdderNModule::disconnectInputPin(int inputPinIndex)
{
  //DEBUG_BREAK;  // check code below
  Module::disconnectInputPin(inputPinIndex);
  int lastPinIndex = getNumInputPins()-1;
  while(!isInputPinConnected(lastPinIndex-1) && lastPinIndex > 1)
  {
    deleteAudioInput(lastPinIndex);
    lastPinIndex--;
  }
}
CREATE_AND_ASSIGN_PROCESSING_FUNCTIONS_N(AdderNModule);

}
