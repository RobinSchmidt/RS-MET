namespace rosic { namespace Sampler {

//===============================================================================================
// rsSamplerEngine::SignalProcessorChain

void SignalProcessorChain::processFrame(rsFloat64x2& inOut)
{
  for(size_t i = 0; i < processors.size(); i++)
    processors[i]->processFrame(inOut);
}

size_t SignalProcessorChain::getNumProcessors(DspType type) const
{
  size_t count = 0;
  for(size_t i = 0; i < processors.size(); i++) {
    if(processors[i]->getType() == type)
      count++; }
  return count;
}

SignalProcessor* SignalProcessorChain::getProcessor(
  DspType type, int index)
{
  int count = 0;  // counts, how many DSPs of given type we have iterated over - why not size_t?
  for(int i = 0; i < (int) processors.size(); i++) {
    SignalProcessor* dsp = getProcessor(i);
    if(dsp->getType() == type) {
      if(count == index)
        return dsp;
      else
        count++;   }}
  return nullptr;
}


/*
void SignalProcessorChain::resetState()
{
for(size_t i = 0; i < processors.size(); i++)
processors[i]->resetState();
}
*/

//===============================================================================================



}}