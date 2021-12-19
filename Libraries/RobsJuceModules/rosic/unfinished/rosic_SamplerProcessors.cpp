namespace rosic {
namespace Sampler {


SignalProcessorPool::SignalProcessorPool()
{
  allocateProcessors();
  // Maybe don't do this on construction. Maybe client code should explicitly request this
}

SignalProcessorPool::~SignalProcessorPool()
{
  deAllocateProcessors();
}

void SignalProcessorPool::deAllocateProcessors()
{
  for(size_t i = 0; i < pool.size(); i++)
    delete pool[i];
  pool.clear();
}

SignalProcessor* SignalProcessorPool::grabProcessor(SignalProcessorType type)
{

  return nullptr;
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


  int dummy = 0;
}

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

/** Class where all the boilerplate for making DSP processors available in the sampler goes. */

class rsSamplerProcessors
{

public:

  class Filter : public SignalProcessor
  {

  public:

    void processFrame(rsFloat64x2& inOut) override {}
    void processBlock(rsFloat64x2* inOut, int N) override {}
    void resetState() override { core.resetState(); }
    void resetSettings() override { core.initCoeffs(); }

  protected:

    rsSamplerFilter core;

  };


  class WaveShaper : public SignalProcessor
  {

  public:

    void processFrame(rsFloat64x2& inOut) override {}
    void processBlock(rsFloat64x2* inOut, int N) override {}
    void resetState() override {}
    void resetSettings() override {}

  protected:

    rsSamplerWaveShaper core;

  };

  // before writing too much boilerplate, fix the API:
  // -processFrame should work either with 2 floats or rsFloat32x4...or maybe make a class
  //  rsFloat32x2...could be just a synonym for rsFloat32x4 (because that can be simdified), but 
  //  the name makes clear that we use only 2 of the 4 available slots
  // -can we avoid the need for the boilerplate? ...or at least reduce the amount? maybe with 
  //  similar strategies as in romos, using macros?



};



void SignalProcessorPool::allocateProcessors()
{
  deAllocateProcessors();    // ...just in case

  // These numbers are preliminary - do something more sensible here:
  int numFilters     = 128;
  int numWaveShapers =  64;
  int numProcessors  = numFilters + numWaveShapers;


  // The actual allocation of the DSP objects:
  using SP = rsSamplerProcessors;
  pool.reserve(numProcessors);
  int i;
  for(i = 0; i < numFilters;     ++i) 
    pool.push_back(new SP::Filter);
  for(i = 0; i < numWaveShapers; ++i) 
    pool.push_back(new SP::WaveShaper);

  int dummy = 0;

  // ToDo: 
  // -Try to find a better way to allocate the memory. We don't want to allocate each object 
  //  seperately. Maybe have a bunch of arrays of direct objects of each specific type of 
  //  processor and the pool only stores pointers to these. we'll see...
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
