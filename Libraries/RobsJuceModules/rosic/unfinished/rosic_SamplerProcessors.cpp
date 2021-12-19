namespace rosic {
namespace Sampler {

//=================================================================================================

void rsSamplerFilter::setup(rsSamplerFilter::Type type, float cutoff, float resonance)
{

}

void rsSamplerFilter::initCoeffs()
{
  /*
  switch(type) // use switch type.topology
  {

  }
  */
}

void rsSamplerFilter::processFrame(float& L, float& R)
{

}

void rsSamplerFilter::resetState()
{

}



//=================================================================================================



//=================================================================================================

SignalProcessorPool::SignalProcessorPool()
{
  allocateProcessors();
  // Maybe don't do this on construction. Maybe client code should explicitly request this
}

SignalProcessorPool::~SignalProcessorPool()
{

}

void SignalProcessorPool::allocateProcessors()
{
  filters.resize(128);
  waveShapers.resize(64);
  // These numbers are preliminary. We need to do something more sensible here later. Perhaps, this 
  // function should be called when a new sfz is loaded and it should have arguments for how many
  // objects of each type are needed. The engine should analyze, how many filters, waveshapers, 
  // etc. could be needed in the worst case, fill a suitable data structure with that information
  // and pass it in.
}

SignalProcessor* SignalProcessorPool::grabProcessor(SignalProcessorType type)
{
  using SPT = SignalProcessorType;
  SignalProcessor* p = nullptr;

  switch(type)
  {
  case SPT::Filter:
  {
    if(!filters.empty())
    {
      p = &filters[filters.size()-1];
      filters.resize(filters.size()-1);
    }
  }
  break;
  case SPT::WaveShaper:
  {
    if(!waveShapers.empty())
    {
      p = &waveShapers[waveShapers.size()-1];
      waveShapers.resize(waveShapers.size()-1);
    }
  }
  break;

  };
  // ToDo: Find a way to reduce the boilerplate! This is unbearable! For each type, we should have
  // one line of code! ..at most. Maybe with a macro? or maybe a template?

  return p;
}
// ToDo: 
// -maybe we should somehow report whether it was successful - if we have no processor of given
//  type available anymore, we may either return a nullptr or a pointer to some sort of dummy
//  processor...a dummy should probably mute the output - bypassing could lead to undesirable 
//  results, like getting a sound through unfiltered or unattenuated which is supposed to be 
//  attenuated via the processor. Using such a muting dummy would gracefully handle overload. 
//  The layer would just be muted, if not enough DSP objects are available
// -maybe somewhere we should let the user control, how many of each processor type should be 
//  pre-allocated - filters are needed a lot, more exotic processors much less so
// -hmm - i think, it's better to return a nullptr, if no processor of desired type is available
//  anymore - the calling code can then forego the whole RegionPlayer


void SignalProcessorPool::returnProcessor(SignalProcessor* p)
{
  using SPT = SignalProcessorType;
  switch(p->getType())
  {
  case SPT::Filter:     
  { 
    filters.resize(filters.size()+1); 
  } break;
  case SPT::WaveShaper: 
  { 
    waveShapers.resize(waveShapers.size()+1); 
  } break;
  }


  int dummy = 0;
  // maybe we should return something?
}








}} // namespaces


//=================================================================================================
/*

Notes:

rsSamplerFilter:
 -Maybe make also a struct for very basic 1-pole/1-zero filter. It can be realized by the biquad 
 (and by the other structures, too), but maybe it's more efficient to do it like that. I expect 
 that a simple 1st order lowpass is a quite common thing to use.
 -Maybe
 -Use the filter also for the equalizer opcode. No need to define a different class for that. Maybe
 extend sfz to support 4 instead of 3 bands when we later can realize 2 bands per filter...



*/
