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

  // Maybe instead of manually resizing everything, do it programmatically, maybe using a
  // SignalProcessorFactory to create the objects. But then: how do we keep track of the different
  // kinds of modules?
}



template<class T> // Grow v by given amount
inline void rsGrow(std::vector<T>& v, size_t amount = 1)
{
  v.resize(v.size() + amount);
}
template<class T> // Shrink v by given amount
inline void rsShrink(std::vector<T>& v, size_t amount = 1)
{
  RAPT::rsAssert(v.size() >= amount);
  v.resize(v.size() - amount);
}
template<class T> // Pointer to last element in v
inline T* rsLastPointer(std::vector<T>& v)
{
  RAPT::rsAssert(!v.empty());
  return &v[v.size()-1];
}
template<class T> // Pointer to last element, shrink by 1
inline T* rsGetLastPtrAndShrink(std::vector<T>& v)
{
  if(v.empty())
    return nullptr;
  T* p = rsLastPointer(v);
  rsShrink(v);
  return p;
}

// maybe move to rapt

SignalProcessor* SignalProcessorPool::grabProcessor(SignalProcessorType type)
{
  using SPT = SignalProcessorType;
  SignalProcessor* p = nullptr;
  switch(type)
  {
  case SPT::Filter:     p = rsGetLastPtrAndShrink(filters);     break;
  case SPT::WaveShaper: p = rsGetLastPtrAndShrink(waveShapers); break;
  };
  return p;
}

void SignalProcessorPool::repositProcessor(SignalProcessor* p)
{
  using SPT = SignalProcessorType;
  switch(p->getType())
  {
  case SPT::Filter:     rsGrow(filters);     break;
  case SPT::WaveShaper: rsGrow(waveShapers); break;
  }
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
