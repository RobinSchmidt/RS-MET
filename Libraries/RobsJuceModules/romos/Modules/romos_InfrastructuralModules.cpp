//#include "romos_InfrastructuralModules.h"

// why is this commented out?
/*
void AudioInputModule::initialize()
{
  initOutputPins(1, rosic::rsString());
  hasHeaderFlag = false;
}
void AudioInputModule::allocateMemory()
{
  // do nothing - the outlying container is resposible for allocation
}
void AudioInputModule::freeMemory()
{
  // do nothing - the outlying container is resposible for deallocation
}
void AudioInputModule::setAudioInputAddress(double *newAddress)
{
  DEBUG_BREAK;
  // update:
  //audioInputs = newAddress;
  //rassert( incomingAudioConnections.empty() ); // input modules should not have incoming connections
}
void AudioInputModule::setAudioOutputAddress(double *newAddress)
{
  DEBUG_BREAK;
  // update:
  //audioOutputs = newAddress;
  //std::vector<AudioConnection*> outgoingAudioConnections = getOutgoingAudioConnections();
  //for(unsigned int i = 0; i < outgoingAudioConnections.size(); i++)
  //  outgoingAudioConnections[i]->updateSourcePointer();
}
*/

//-------------------------------------------------------------------------------------------------

/*
void AudioOutputModule::initialize()
{
  initInputPins(1, rosic::rsString());
  hasHeaderFlag = false;
}
*/

void AudioOutputModule::allocateMemory()
{
  // do nothing - the outlying container is responsible for allocation
}
void AudioOutputModule::freeMemory()
{
  // do nothing - the outlying container is responsible for deallocation
}
/*
void AudioOutputModule::resetState()
{

}
*/

//-------------------------------------------------------------------------------------------------

void SystemSampleRateModule::initialize()
{
  initOutputPins({ "SampleRate" });
}
INLINE void SystemSampleRateModule::process(Module *module, double *out, int voiceIndex)
{
  *out = processingStatus.getSystemSampleRate();
}
CREATE_AND_ASSIGN_PROCESSING_FUNCTIONS_0(SystemSampleRateModule);

//-------------------------------------------------------------------------------------------------

void SystemSamplePeriodModule::initialize()
{
  initOutputPins({ "SamplePeriod" });
}
INLINE void SystemSamplePeriodModule::process(Module *module, double *out, int voiceIndex)
{
  *out = processingStatus.getSystemSamplePeriod();
}
CREATE_AND_ASSIGN_PROCESSING_FUNCTIONS_0(SystemSamplePeriodModule);

//-------------------------------------------------------------------------------------------------

void NoteGateModule::initialize()
{
  //initOutputPins({ "Gate"});  // no - messes up block
  initOutputPins({ "" });
  hasHeaderFlag = false;
}
INLINE void NoteGateModule::process(Module *module, double *out, int voiceIndex)
{
  *out = (double)voiceAllocator.isNoteOn(voiceIndex);
}
CREATE_AND_ASSIGN_PROCESSING_FUNCTIONS_0(NoteGateModule);

//-------------------------------------------------------------------------------------------------

void NoteOnTriggerModule::initialize()
{
  //initOutputPins({ "Ping"});
  initOutputPins({ "" });
  hasHeaderFlag = false;
}
INLINE void NoteOnTriggerModule::process(Module *module, double *out, int voiceIndex)
{
  *out = voiceAllocator.getNoteOnTriggerFlag(voiceIndex);
}
CREATE_AND_ASSIGN_PROCESSING_FUNCTIONS_0(NoteOnTriggerModule);

//-------------------------------------------------------------------------------------------------

void NoteOffTriggerModule::initialize()
{
  initOutputPins({ "" });
  hasHeaderFlag = false;
}
INLINE void NoteOffTriggerModule::process(Module *module, double *out, int voiceIndex)
{
  *out = voiceAllocator.getNoteOffTriggerFlag(voiceIndex);
}
CREATE_AND_ASSIGN_PROCESSING_FUNCTIONS_0(NoteOffTriggerModule);

//-------------------------------------------------------------------------------------------------

void VoiceKillerModule::initialize()
{
  initInputPins({ "In" });
  addParameter(rosic::rsString("Threshold"), "-100.0");
  addParameter(rosic::rsString("TimeOut"), "0.01");
  parameterChanged(0);   // to init internal variables threshold, timeOut
}
INLINE void VoiceKillerModule::process(Module *module, double *in, double *out, int voiceIndex)
{
  if(!voiceAllocator.isVoicePlaying(voiceIndex))
    return;  // because it is always called for voice 0 in monophonic mode

  VoiceKillerModule *voiceKiller = static_cast<VoiceKillerModule*> (module);

  if(fabs(*in) >= voiceKiller->threshold)
    voiceKiller->sampleCounters[voiceIndex] = 0;
  else
  {
    voiceKiller->sampleCounters[voiceIndex] += 1;
    if((double)voiceKiller->sampleCounters[voiceIndex] > 
      voiceKiller->timeOut * processingStatus.getSystemSampleRate())
    {
      voiceAllocator.killVoice(voiceIndex);
      //voiceKiller->sampleCounters[voiceIndex] = 0;
      module->getTopLevelModule()->resetVoiceState(voiceIndex);
    }
  }
}
void VoiceKillerModule::resetVoiceState(int voiceIndex)
{
  AtomicModule::resetVoiceState(voiceIndex);
  sampleCounters[voiceIndex] = 0;
}
void VoiceKillerModule::parameterChanged(int index)
{
  threshold = RAPT::rsDbToAmp(parameters[0].value.asDouble());
  timeOut   = parameters[1].value.asDouble();
}
void VoiceKillerModule::allocateMemory()
{
  AtomicModule::allocateMemory();
  sampleCounters = new unsigned int[getNumVoices()];
}
void VoiceKillerModule::freeMemory()
{
  AtomicModule::freeMemory();
  delete[] sampleCounters;
  sampleCounters = NULL;
}
CREATE_AND_ASSIGN_PROCESSING_FUNCTIONS_1(VoiceKillerModule);

//-------------------------------------------------------------------------------------------------

void VoiceCombinerModule::initialize()
{
  initInputPins({ "" });
  initOutputPins({ "" });
  hasHeaderFlag = false;
}
INLINE void VoiceCombinerModule::process(Module *module, double *ins, double *outs, int voiceIndex)
{
  outs[0] = ins[0];  // function actually not used
}
void VoiceCombinerModule::processMonoFrame(Module *module, int voiceIndex)
{
  module->audioOutputs[0] = *(module->inputPins[0].outputPointer);
}
void VoiceCombinerModule::processPolyFrame(Module *module, int voiceIndex)
{
  module->audioOutputs[0] = 0.0;
  for(int playIndex = 0; playIndex < voiceAllocator.getNumPlayingVoices(); playIndex++)
  {
    int voiceIndex  = voiceAllocator.getPlayingVoiceIndices()[playIndex];
    module->audioOutputs[0] += *(module->inputPins[0].outputPointer 
      + voiceIndex * module->inputPins[0].outputVoiceStride);
    // can be streamlined (like in the macros)
  }
}
void VoiceCombinerModule::processMonoBlock(Module *module, int voiceIndex, int blockSize)
{
  for(int frameIndex = 0; frameIndex < blockSize; frameIndex++)
    module->audioOutputs[frameIndex] = *(module->inputPins[0].outputPointer 
      + module->inputPins[0].outputFrameSize);
}
void VoiceCombinerModule::processPolyBlock(Module *module, int voiceIndex, int blockSize)
{
  for(int frameIndex = 0; frameIndex < blockSize; frameIndex++)
  {
    module->audioOutputs[frameIndex] = 0.0;
    double *sourcePointer = module->inputPins[0].outputPointer 
      + frameIndex * module->inputPins[0].outputFrameSize;
    for(int playIndex = 0; playIndex < voiceAllocator.getNumPlayingVoices(); playIndex++)
    {
      int voiceIndex  = voiceAllocator.getPlayingVoiceIndices()[playIndex];
      module->audioOutputs[frameIndex] += 
        *(sourcePointer + voiceIndex * module->inputPins[0].outputVoiceStride);
    }
  }
}
void VoiceCombinerModule::allocateMemory()
{
  AtomicModule::allocateMemory();
   // \todo we should always allocate only memory for a single voice in the outputs
}
void VoiceCombinerModule::freeMemory()
{
  AtomicModule::freeMemory();
}
void VoiceCombinerModule::setPolyphonic(bool shouldBePolyphonic)
{
  // do nothing - voice-combiners are always monophonic (on their output side)
}
void VoiceCombinerModule::clearVoiceBuffer(int voiceIndex)
{
  // do nothing - clearing the buffer would lead to gaps in the output signal (with the length of 
  // bufferSize) whenever voice 0 gets killed
}
void VoiceCombinerModule::connectInputPinTo(int inputPinIndex, Module *sourceModule, 
  int sourceOutputPinIndex)
{
  AtomicModule::connectInputPinTo(inputPinIndex, sourceModule, sourceOutputPinIndex);
  assignProcessingFunctions();
}
void VoiceCombinerModule::updateInputPointersAndInFrameStrides()
{
  // this function is called whenever some potential input-source module has changed its output 
  // polyphony - we use this callback to set up our processing functions:
  AtomicModule::updateInputPointersAndInFrameStrides();
  assignProcessingFunctions();
}
void VoiceCombinerModule::assignProcessingFunctions()
{
  // when the source is monophonic, we use our monophonic processing function (that just passes 
  // through the 0-th voice signal), when the source is polyphonic, we use our polyphonic procesing
  // function (which sums over the voices):
  if(inputPins[0].sourceModule != NULL)
  {
    if(inputPins[0].sourceModule->isPolyphonic())
    {
      processFrame = &VoiceCombinerModule::processPolyFrame;
      processBlock = &VoiceCombinerModule::processPolyBlock;
    }
    else
    {
      processFrame = &VoiceCombinerModule::processMonoFrame;
      processBlock = &VoiceCombinerModule::processMonoBlock;
    }
  }
  else
  {
    processFrame = &VoiceCombinerModule::processMonoFrame;
    processBlock = &VoiceCombinerModule::processMonoBlock;
  }
}

//-------------------------------------------------------------------------------------------------

void NoteFrequencyModule::initialize()
{
  initOutputPins({ "" });
  hasHeaderFlag = false;
}
INLINE void NoteFrequencyModule::process(Module *module, double *out, int voiceIndex)
{
  *out = (double)voiceAllocator.getFrequencyOfVoice(voiceIndex);
}
CREATE_AND_ASSIGN_PROCESSING_FUNCTIONS_0(NoteFrequencyModule);

//-------------------------------------------------------------------------------------------------

void NoteVelocityModule::initialize()
{
  initOutputPins({ "" });
  hasHeaderFlag = false;
}
INLINE void NoteVelocityModule::process(Module *module, double *out, int voiceIndex)
{
  *out = voiceAllocator.getNormalizedVelocityOfVoice(voiceIndex);
}
CREATE_AND_ASSIGN_PROCESSING_FUNCTIONS_0(NoteVelocityModule);

//-------------------------------------------------------------------------------------------------

void ParameterModule::initialize()
{
  ConstantModule::initialize();
  minValue           = 0.0;
  maxValue           = 1.0;
  value              = 0.5;
  defaultValue       = 0.5;
  quantization       = 0.0;
  mappingFunction    = LINEAR_MAPPING;
  assignedController = -1; // none

                                          // index in parameters array
  addParameter("Value", "0.5");           // 0
  addParameter("MinValue", "0.0");        // 1
  addParameter("MaxValue", "1.0");        // 2
  addParameter("DefaultValue", "0.5");    // 3
  addParameter("Scaling", "Linear");      // 4 other possible values: exponential, power, etc.
  addParameter("Quantization", "0.0");    // 5
  addParameter("Controller", "0");        // 6    
  addParameter("ControlRangeMin", "0.0"); // 7
  addParameter("ControlRangeMax", "0.0"); // 8
  addParameter("Smoothing", "0.0");       // 9
  addParameter("HelpText", "Parameter");  // 10

  parameterChanged(0); // updates internal variables that depend on the 0th ("Value") parameter
  parameterChanged(1); // causes all parameter-dependent variables except the 0th ("Value") to be 
                       // updated at once
}
void ParameterModule::resetVoiceState(int voiceIndex)
{
  AtomicModule::resetVoiceState(voiceIndex);
  // later: set current-value to target value here (smoothing)
}
void ParameterModule::parameterChanged(int index)
{
  switch(index)
  {
  case 0:
  {
    value = RAPT::rsClip(parameters[0].value.asDouble(), minValue, maxValue);
    setParameter(0, rosic::rsString(value), false); // make the string reflect the possibly modified value
  }
  break;
  case 1:
  {
    minValue     = RAPT::rsMin(parameters[1].value.asDouble(), maxValue);
    value        = RAPT::rsClip(value, minValue, maxValue);
    defaultValue = RAPT::rsClip(defaultValue, minValue, maxValue);
    setParameter(1, rosic::rsString(minValue), false);
    setParameter(0, rosic::rsString(value), false);
    setParameter(3, rosic::rsString(defaultValue), false);
  }
  break;
  case 2:
  {
    maxValue     = RAPT::rsMax(parameters[2].value.asDouble(), minValue);
    value        = RAPT::rsClip(value, minValue, maxValue);
    defaultValue = RAPT::rsClip(defaultValue, minValue, maxValue);
    setParameter(2, rosic::rsString(maxValue), false);
    setParameter(0, rosic::rsString(value), false);
    setParameter(3, rosic::rsString(defaultValue), false);
  }
  break;
  case 3:
  {
    defaultValue = RAPT::rsClip(parameters[3].value.asDouble(), minValue, maxValue);
    setParameter(3, rosic::rsString(defaultValue), false);
  }
  break;
  case 4:
  {
    if(parameters[4].value == "Exponential")
      setMinMaxAndMapping(minValue, maxValue, EXPONENTIAL_MAPPING);
    else
      setMinMaxAndMapping(minValue, maxValue, LINEAR_MAPPING);
  }
  break;

  default:
  {
    //DEBUG_BREAK; // unknown parameter index
  }

  }

  int dummy = 0;
}
void ParameterModule::setModuleName(const std::string& newName)
{
  // overriden to avoid setting the value according to the name (which the ConstantModule 
  // baseclass would do)
  AtomicModule::setModuleName(newName);
}

/*
void ParameterModule::setValue(double newValue)
{
  value = clip(newValue, minValue, maxValue);
  setParameter(0, rosic::rsString(value), false);  // maybe make this optional - bool parameter to the function
}
void ParameterModule::setMinValue(double newMinValue)
{
  minValue = newMinValue;
  enforceConsistencyOfValues();
}
void ParameterModule::setMaxValue(double newMaxValue)
{
  maxValue = newMaxValue;
  enforceConsistencyOfValues();
}
void ParameterModule::setDefaultValue(double newDefaultValue)
{
  defaultValue = newDefaultValue;
  enforceConsistencyOfValues();
}
void ParameterModule::setQuantization(double newQuantization)
{
  DEBUG_BREAK; // optional value quantization not yet implemented
}
*/

void ParameterModule::setMinMaxAndMapping(double newMin, double newMax, int newMapping)
{
  if(newMin <= 0.0 && newMapping == EXPONENTIAL_MAPPING)
    newMapping = LINEAR_MAPPING; // min must be > 0.0 for exponential mapping
  mappingFunction = newMapping;

  if(newMin <= newMax)
  {
    minValue = newMin;
    maxValue = newMax;
  }
  else
    minValue = maxValue = newMin;

  value        = RAPT::rsClip(value, minValue, maxValue);
  defaultValue = RAPT::rsClip(defaultValue, minValue, maxValue);

  setParameter(0, rosic::rsString(value), false);
  setParameter(1, rosic::rsString(minValue), false);
  setParameter(2, rosic::rsString(maxValue), false);
  setParameter(3, rosic::rsString(defaultValue), false);
  if(mappingFunction == EXPONENTIAL_MAPPING)
    setParameter(4, "Exponential", false);
  else
    setParameter(4, "Linear", false);
}

void ParameterModule::setValueFromController(double controllerValue)
{
  DEBUG_BREAK; // controller stuff not yet implemented
}
void ParameterModule::setValueFromSnapshots(int topLeftIndex, int topRightIndex, 
  int bottomLeftIndex, int bottomRightIndex, double x, double y)
{
  DEBUG_BREAK; // snapshot morphing not yet implemented
}

void ParameterModule::enforceConsistencyOfValues()
{
  maxValue     = RAPT::rsMax(minValue, maxValue);
  minValue     = RAPT::rsMin(minValue, maxValue);
  value        = RAPT::rsClip(value, minValue, maxValue);
  defaultValue = RAPT::rsClip(defaultValue, minValue, maxValue);

  /*
  // ... we also need to quantize the value when interval != 0.0
  // update the parameter value strings, too:
  setParameter(0, rosic::rsString(value),        false);
  setParameter(1, rosic::rsString(minValue),     false);
  setParameter(2, rosic::rsString(maxValue),     false);
  setParameter(3, rosic::rsString(defaultValue), false);
  */
}

double ParameterModule::mapNormalizedValue(double normalizedValue)
{
  switch(mappingFunction)
  {
  case EXPONENTIAL_MAPPING: return RAPT::rsLinToExp(normalizedValue, 0.0, 1.0, minValue, maxValue);
  default:                  return RAPT::rsLinToLin(normalizedValue, 0.0, 1.0, minValue, maxValue);
  }
}

double ParameterModule::unmapToNormalizedValue(double mappedValue)
{
  switch(mappingFunction)
  {
  case EXPONENTIAL_MAPPING: return RAPT::rsExpToLin(mappedValue, minValue, maxValue, 0.0, 1.0);
  default:                  return RAPT::rsLinToLin(mappedValue, minValue, maxValue, 0.0, 1.0);
  }
}
