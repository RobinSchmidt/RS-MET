namespace rosic {
namespace Sampler {

void SignalProcessor::addParameter(Opcode opcode)
{
  params.push_back(Parameter(opcode));
}

void SignalProcessor::setParameter(Opcode opcode, float value)
{
  size_t i = 0;
  for(i = 0; i < params.size(); i++) {
    if(params[i].getOpcode() == opcode) {
      params[i].setValue(value);
      return; }}
  RAPT::rsError("Parameter not found in SignalProcessor::setParameter");
}

//=================================================================================================

void rsSamplerFilter::setup(rsSamplerFilter::Type type, float w, float reso)
{
  this->type = type;
  using FO = RAPT::rsOnePoleFilter<float, float>;
  using BQ = RAPT::rsBiquadDesigner;  // maybe it should have a template parameter?

  // Preliminary to cater for the API of rsBiquadDesigner - which should really be changed...
  static const float s = float(1/(2*PI));
  float Q = 1.f / sqrt(2.f);    // preliminary - todo: find a conversion formula reso2q

  FilterImpl& i = impl;  // as abbreviation
  switch(type)
  {
  case Type::FO_Lowpass:  FO::coeffsLowpassIIT( w, &i.fo.b0, &i.fo.b1, &i.fo.a1); return;
  case Type::FO_Highpass: FO::coeffsHighpassMZT(w, &i.fo.b0, &i.fo.b1, &i.fo.a1); return;

    // This API sucks - fix it! The functions should take w for the frequency (not freq in Hz and 
    // sample-rate), their names should be much shorter (e.g. coeffsLowpassRBJ or just lowpassRBJ),
    // output params should come last and be passed as pointers.
  case Type::BQ_Lowpass:  BQ::calculateCookbookLowpassCoeffs(
    i.bqd.b0, i.bqd.b1, i.bqd.b2, i.bqd.a1, i.bqd.a2, 1.f, s*w, Q); return;
  case Type::BQ_Highpass: BQ::calculateCookbookHighpassCoeffs(
    i.bqd.b0, i.bqd.b1, i.bqd.b2, i.bqd.a1, i.bqd.a2, 1.f, s*w, Q); return;


  }
  RAPT::rsError("Unknown filter type in rsSamplerFilter::setup");

  // ToDo:
  // -Organize the Type enum in such a way that we can retrieve the filter topology from it via 
  //  bitmasking such that we do not need a branch for every type but only one for every topology.
  //  This will reduce the boilerplate a lot.
  // -Make a consistent choice for all of RAPT whether a recursion coeffs of filters should have a 
  //  minus sign or not and update all code accordingly. Make benchmarks what is faster, ask at KVR
  //  what others do.
}

/*
void rsSamplerFilter::initCoeffs()
{

}

void rsSamplerFilter::updateCoeffs()
{

}
*/

void rsSamplerFilter::processFrame(float& L, float& R)
{
  TSig io(L, R);
  FilterImpl& i = impl;
  switch(type)
  {
  case Type::FO_Lowpass:  io = i.fo.getSample(io); break;
  case Type::FO_Highpass: io = i.fo.getSample(io); break;

  case Type::BQ_Lowpass:  io = i.bqd.getSample(io); break;
  case Type::BQ_Highpass: io = i.bqd.getSample(io); break;
  };
  L = io.x; // Preliminary - as long as we are abusing rsVector2D for the signal
  R = io.y;

  // ...later, we want to use a simd type and retrieve the elements like so:
  //L = io[0]; R = io[1];
}

void rsSamplerFilter::resetState()
{
  FilterImpl& i = impl;
  switch(type)
  {
  case Type::FO_Lowpass:  i.fo.resetState();  return;
  case Type::FO_Highpass: i.fo.resetState();  return;

  case Type::BQ_Lowpass:  i.bqd.resetState(); return;
  case Type::BQ_Highpass: i.bqd.resetState(); return;
  }
  RAPT::rsError("Unknown filter type in rsSamplerFilter::resetState");
}

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
                        // Debug  Release  ...these values make sense for development
  filters.init(4);      //   4      64
  waveShapers.init(4);  //   8     128
  // These numbers are preliminary. We need to do something more sensible here later. Perhaps, this 
  // function should be called when a new sfz is loaded and it should have arguments for how many
  // objects of each type are needed. The engine should analyze, how many filters, waveshapers, 
  // etc. could be needed in the worst case, fill a suitable data structure with that information
  // and pass it in.

  // Maybe instead of manually resizing everything, do it programmatically, maybe using a
  // SignalProcessorFactory to create the objects. But then: how do we keep track of the different
  // kinds of modules?
}


/*
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
*/
// maybe move to rapt ...hmm...not needed anymore...perhaps delete...not sure, if they are useful
// enough in general to include in the library

SignalProcessor* SignalProcessorPool::grabProcessor(DspType type)
{
  using SPT = DspType;
  SignalProcessor* p = nullptr;
  switch(type)
  {
  case SPT::Filter:     p = filters.grabItem();     break;
  case SPT::WaveShaper: p = waveShapers.grabItem(); break;
  };
  return p;
}

void SignalProcessorPool::repositProcessor(SignalProcessor* p)
{
  using SPT = DspType;
  int i = -1;
  switch(p->getType())
  {
  case SPT::Filter:     i = filters.repositItem(p);     break;
  case SPT::WaveShaper: i = waveShapers.repositItem(p); break;
  }
  RAPT::rsAssert(i != -1, "Reposited processor was not in pool");
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
