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
  RAPT::rsError("Parameter not found in DSP");
}

//=================================================================================================

void rsSamplerFilter::setup(rsSamplerFilter::Type type, float w, float reso)
{
  this->type = type;

  using P1Z1 = RAPT::rsOnePoleFilter<float, float>;

  switch(type)
  {
  case Type::Lowpass_6:
  {
    P1Z1::coeffsLowpassIIT(w, &vars.p1z1.b0, &vars.p1z1.b1, &vars.p1z1.a1);
  } break;

  }


  int dummy = 0;
}

void rsSamplerFilter::initCoeffs()
{

}

void rsSamplerFilter::updateCoeffs()
{

}

void rsSamplerFilter::processFrame(float& L, float& R)
{
  TSig io(L, R);
  switch(type)
  {
  case Type::Lowpass_6: io = vars.p1z1.getSample(io); break;
  };

  // Preliminary - as long as we are abusing rsVector2D for the signal:
  L = io.x;
  R = io.y;

  // ...later, we want to use a simd type and retrieve the elements like so:
  //L = io[0]; R = io[1];


  int dummy = 0;
}

void rsSamplerFilter::resetState()
{
  switch(type)
  {
  case Type::Lowpass_6: { vars.p1z1.resetState(); } break;
  }
  int dummy = 0;
}

//-------------------------------------------------------------------------------------------------





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
  filters.init(4);
  waveShapers.init(4);

                          // Debug  Release  ...these values make sense for development
  //waveShapers.resize(4);  //   4      64
  //filters.resize(8);      //   8     128
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
// maybe move to rapt

SignalProcessor* SignalProcessorPool::grabProcessor(SignalProcessorType type)
{
  using SPT = SignalProcessorType;
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
  using SPT = SignalProcessorType;
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
