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

void SignalProcessor::resetSettings(int index)
{
  SfzCodeBook* cb = SfzCodeBook::getInstance();
  for(size_t i = 0; i < params.size(); i++) {
    Opcode op = params[i].getOpcode();
    float defVal = cb->opcodeDefaultValue(op, index);
    params[i].setValue(defVal); }
}

//=================================================================================================

void rsSamplerFilter::setupCutRes(rsSamplerFilter::Type type, float w, float resoGainDb)
{
  if(type == Type::Unknown || w == 0.f)
    type = Type::Bypass;
    // When the cutoff is not defined, it defaults to zero. In this case, we switch into bypass
    // mode. We also default to bypass if mode is not set...hmm...maybe we should default to
    // lpf_12 in this case? ...but only if cutoff is nonzero...we'll see...

  using namespace RAPT;

  this->type = type;
  using FO = rsOnePoleFilter<float, float>;
  using BQ = rsBiquadDesigner;  // maybe it should have a template parameter?

  static const float s = float(1/(2*PI));
  // Preliminary to cater for the API of rsBiquadDesigner - ToDo: change API (maybe write a new 
  // class fo that and deprecate the old)

  // Compute the desired filter quality factor Q from the resonance gain in dB:
  float A = rsDbToAmp(resoGainDb);  // Raw resonance amplitude
  float Q = A;                      // This is correct for 2nd order bandpass filters...
  if(type == Type::BQ_Lowpass || type == Type::BQ_Highpass) // ..low- and highpass filters need..
    Q = rsBandwidthConverter::lowpassResoGainToQ(A);        // ..a more complicated formula

  FilterImpl& i = impl;  // as abbreviation
  switch(type)
  {
  case Type::Bypass: FO::coeffsBypass(&i.fo.b0, &i.fo.b1, &i.fo.a1); return;
    // maybe use this as default branch

  case Type::FO_Lowpass:  FO::coeffsLowpassIIT( w, &i.fo.b0, &i.fo.b1, &i.fo.a1); return;
  case Type::FO_Highpass: FO::coeffsHighpassMZT(w, &i.fo.b0, &i.fo.b1, &i.fo.a1); return;

  // This API sucks - fix it! The functions should take w for the frequency (not freq in Hz and 
  // sample-rate), their names should be much shorter (e.g. coeffsLowpassRBJ or just lowpassRBJ),
  // output params should come last and be passed as pointers.
  case Type::BQ_Lowpass:  BQ::calculateCookbookLowpassCoeffs(
    i.bqd.b0, i.bqd.b1, i.bqd.b2, i.bqd.a1, i.bqd.a2, 1.f, s*w, Q); return;
  case Type::BQ_Highpass: BQ::calculateCookbookHighpassCoeffs(
    i.bqd.b0, i.bqd.b1, i.bqd.b2, i.bqd.a1, i.bqd.a2, 1.f, s*w, Q); return;
  case Type::BQ_Bandpass_Skirt: BQ::calculateCookbookBandpassConstSkirtCoeffsViaQ(
    i.bqd.b0, i.bqd.b1, i.bqd.b2, i.bqd.a1, i.bqd.a2, 1.f, s*w, Q); return;
  case Type::BQ_Bandstop: BQ::calculateCookbookBandrejectCoeffsViaQ(
    i.bqd.b0, i.bqd.b1, i.bqd.b2, i.bqd.a1, i.bqd.a2, 1.f, s*w, Q); return;

  }
  RAPT::rsError("Unknown filter type in rsSamplerFilter::setupCutRes");

  // ToDo:
  // -Make a consistent choice for all of RAPT whether recursion coeffs of filters should have a 
  //  minus sign or not and update all code accordingly. Be careful - this change ripples through 
  //  all products - maybe introduce new names for the functions and deprecate the old ones instead
  //  of just changing their code. Make benchmarks what is faster, ask at KVR what others do.
  // -Optimize: Compute resonance related stuff only when applicable. ...but maybe we should have 
  //  a separate class for first order filters anyway to save memory
  // -Figure out, if sfz+ and other implementations also use the exact formula for Q for lowpass 
  //  and highpass filters. It's quite expensive to compute and the simple identity function that 
  //  works  perfectly for bandpass gives a reasonable approximation for low- and highpass, too, 
  //  especially as Q gets larger. At low Q, it would give a little bit of extra resonance.
  // -Figure out, if the Q = A formula should also be used for bandreject filters. At the moment,
  //  we just do it, but i'm not sure, if that's the right thing to do. It seems plausible, though.
}

void rsSamplerFilter::setupGainFreqBw(Type type, float gainDb, float w, float bw)
{
  using namespace RAPT;
  rsAssert(type == Type::BQ_Bell);
  rsAssert(bw > 0.f, "Bandwidth must be positive in setupGainFreqBw");
  rsAssert( w > 0.f, "Omega must be positive in setupGainFreqBw");
  // ToDo: allow w=0 later

  this->type = type;
  static const float s = float(1/(2*PI));
  float rawGain = rsDbToAmp(gainDb);
  FilterImpl& i = impl;  // as abbreviation
  rsBiquadDesigner::calculatePrescribedNyquistGainEqCoeffs(
    i.bqd.b0, i.bqd.b1, i.bqd.b2, i.bqd.a1, i.bqd.a2, 1.f, s*w, bw, rawGain, 1.f);


  /*
  float k = powf(2.f, 0.5f*bw);
  float Q = k / (k*k - 1.f);       // Q = 2^(bo/2) / (2^bo - 1)
  // ToDo: verify formula - if correct, move to RAPT::rsBandwidthConverter and call it like:
  //float Q = rsBandwidthConverter::octavesToQ(bw);


  using BQ = rsBiquadDesigner;

  rsBiquadDesigner::calculateCookbookPeakFilterCoeffsViaQ(
    i.bqd.b0, i.bqd.b1, i.bqd.b2, i.bqd.a1, i.bqd.a2, 1.f, s*w, Q, rawGain);
    */

  //rsError("Unknown filter type in rsSamplerFilter::setupGainFreqBw");
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
  case Type::Bypass: break;

  case Type::FO_Lowpass:  io = i.fo.getSample(io); break;
  case Type::FO_Highpass: io = i.fo.getSample(io); break;

  case Type::BQ_Lowpass:        io = i.bqd.getSample(io); break;
  case Type::BQ_Highpass:       io = i.bqd.getSample(io); break;
  case Type::BQ_Bandpass_Skirt: io = i.bqd.getSample(io); break;
  case Type::BQ_Bandstop:       io = i.bqd.getSample(io); break;
  case Type::BQ_Bell:           io = i.bqd.getSample(io); break;

  };
  L = io.x; // Preliminary - as long as we are abusing rsVector2D for the signal
  R = io.y;

  // ...later, we want to use a simd type and retrieve the elements like so:
  //L = io[0]; R = io[1];

  // ToDo:
  // -Organize the Type enum in such a way that we can retrieve the filter topology from it via 
  //  bitmasking such that we do not need a branch for every type but only one for every topology.
  //  This will reduce the boilerplate a lot.
}

void rsSamplerFilter::resetState()
{
  FilterImpl& i = impl;
  switch(type)
  {
  case Type::Bypass: return;

  case Type::FO_Lowpass:  i.fo.resetState();  return;
  case Type::FO_Highpass: i.fo.resetState();  return;

  case Type::BQ_Lowpass:        i.bqd.resetState(); return;
  case Type::BQ_Highpass:       i.bqd.resetState(); return;
  case Type::BQ_Bandpass_Skirt: i.bqd.resetState(); return;
  case Type::BQ_Bandstop:       i.bqd.resetState(); return;
  case Type::BQ_Bell:           i.bqd.resetState(); return;

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
  filters.init(64);
  equalizers.init(64); 
  waveShapers.init(8);
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
// Maybe move to rapt ...hmm...not needed anymore...perhaps delete...not sure, if they are useful
// enough in general to include in the library. Maybe create a file where we can deposit code that
// may become useful later

SignalProcessor* SignalProcessorPool::grabProcessor(DspType type)
{
  using SPT = DspType;
  SignalProcessor* p = nullptr;
  switch(type)
  {
  case SPT::Filter:     p = filters.grabItem();     break;
  case SPT::Equalizer:  p = equalizers.grabItem();  break;
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
  case SPT::Equalizer:  i = equalizers.repositItem(p);  break;
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

 -to figure out the correct design formulas, check source code of sfizz and linuxsampler
  https://github.com/sfztools/sfizz/blob/develop/src/sfizz/dsp/filters/filters_modulable.dsp
  https://github.com/sfztools/sfizz/blob/develop/src/sfizz/dsp/filters/rbj_filters.dsp

  https://github.com/sfztools/sfizz/blob/develop/src/sfizz/dsp/filters/sfz_filters.dsp
  sfzGetQFromSlope
  cutoff = hslider("[01] Cutoff [unit:Hz] [scale:log]", 440.0, 50.0, 10000.0, 1.0) : max(1.0) : min(20000.0);
  Q = vslider("[02] Resonance [unit:dB]", 0.0, 0.0, 40.0, 0.1) : max(-60.0) : min(60.0) : ba.db2linear;
  ...

  https://github.com/sfztools/sfizz/blob/develop/src/sfizz/gen/filters/sfzPeq.hxx
  https://github.com/sfztools/sfizz/blob/develop/src/sfizz/gen/filters/sfzLpf2p.hxx
  https://github.com/sfztools/sfizz/blob/develop/src/sfizz/gen/filters/sfzHpf2p.hxx

  https://github.com/sfztools/sfizz/blob/develop/src/sfizz/effects/Filter.cpp
  _filter.prepare(_desc.cutoff, _desc.resonance, _desc.gain);
  ..look for that function ...  _filter is of class sfz::Filter

  https://github.com/sfztools/sfizz/blob/acd866fd3d247d2fc659593cac96e88e801c29e2/src/sfizz/SfzFilter.h
  https://github.com/sfztools/sfizz/blob/acd866fd3d247d2fc659593cac96e88e801c29e2/src/sfizz/SfzFilter.cpp


  https://www.electronics-tutorials.ws/filter/second-order-filters.html


  https://github.com/linuxsampler/linuxsampler/blob/master/src/engines/common/BiquadFilter.h


rsSamplerWaveShaper:
-Try to find a similar parametrization for each of the available shapes.
-One parameter should be interpretable as "hardness". it should go from 0 to 100% where 100% maps 
 to a hard-clipper and 0% should perhaps be the default, suing the prototype sigmoid without 
 modification like using e.g. tanh(drive * x) directly without applying a "hardening" 
 transformation.
 -for tanh, maybe try functions like tanh(a*x^b), tanh(a*x+b*x^3), tanh(sinh(a*x)/x), 
  tanh((a*x+b*x^3)/(c*x+d*x^3))


ToDo:
-Maybe at some later stage, generalize the dspChain in the sampler engine to a more flexible 
 dspGraph where the routing is not necessarily fixed to x -> 1 -> 2 -> 3 -> y, say (where x,y stand
 for input and output and the numbers stand for the dsps). Maybe the routing could be controlled in 
 the sfz file via a new opcode dsp_graph_... where the ... could be a string like x>1_1>2_2>3_3>y 
 to build the simple chain. The string x>1_x>2_1>3_2>3_3>y could mean: input goes into 1,2 in 
 paralell, 1,2 get summed and go into 3, 3 goes to out. Maybe using > could be problematic when 
 embedding such strings in xml files? maybe use x-1_x-2_1-3_2-3_3-y instead? Or some other symbol?
 Maybe we need a constraint that the left number (source node) must always be less than the right 
 number (target node). This will avoid loops and presumably make parsing easier. Maybe the 
 restriction can be lifted later...
-We actually already need a similar syntax to specify the wiring of the modulations. Maybe 
 something like lfo2_to_eq5_freq, eg2_to_cutoff3, etc. ...so, maybe for the DSPs, we should have 
 one opcode per wire instead of specifying the complete wiring in a single opcode? Maybe something
 like in_to_fil1=1, fil1_to_fil2=1, fil2_to_out=1 specifies the chain in -> fil1 -> fil2 -> out and 
 that wiring is implicitly assumed if nothing else is said?
-Test the computational load incurred by the on-demand assembling of the DSP chain. Maybe 
 pre-assemble the chains for all notes, if it turns out to be a performance problem. That would 
 increase the RAM usage beause we would have to keep around all those pre-built DSP chains.


Ideas:
-Maybe use std::variant for the DSP chain, see here: https://www.youtube.com/watch?v=h-zy1hBqT74
 towards the end. The same thing should actually also be possible with unions and then maybe with
 a struct consisting of such a union and a function-pointer (or a bunch of them) for the process
 etc. functions.
 -This may let us do away with the pre-allocation of specific DSP devices. We could perhaps use 
  an array or vector of direct objects for the DSPs instead of pointers (the objects would be the
  variants or unions). That would also increase the locality of the data used in the DSP chain. 
  -> good!
 -But we would perhaps have to commit to a certain maximum length of the DSP chain and the memory 
  usage of a variant DSP would be dictated by the maximum size of all of them which may be 
  wasteful. It would also be wasteful because the actual length may always be equal to the maximum
  length.
  -> bad!
-We could also have a general pool of pre-allocated memory and let the DSP chain use a std::vector
 with a custom allocator tapping into that pool.
 -That would require to write a custom memory allocator. Such an allocator could pre-allocate a 
  certain amount of memory (in bytes), then maybe split that into buckets of different sizes

Notes:
-Different effects have different requirements on the resources. At least two relevant resources
 to consider are CPU usage and RAM usage mapping to the classical time- and space-requirements in
 computer science. Here is a table of what I would expect from theoretical considerations:

   DSP Process                        CPU   RAM

   simple filter (biquad, etc.)        L     L     L: low, H: high, M: medium
   gain, pan, width, tremolo, wah      L     L
   naive distortion (waveshaping)      L     L
   flanger, vibrato                    L     M
   echo, delay                         L     H

   anti-aliased distortion             M     L
   ringmod, freq-shift, dynamics       M     L
   chorus, pitch-shift                 M     M
   dynamics with look-ahead            M     M
   adaptive filters (LMS, etc)         M     M
   simple reverb                       M     H

   complex filter (VA, ZDF, etc.)      H     L
   spectral stuff (FFT)                H     M
   good reverb, convolution            H     H

-Another thing to consider is the possible introduction of latency and/or the spikeyness of 
 CPU-load distribution due to block-based processing. Anything with look-ahead will introduce 
 latency. Pitch-shifters will also. FFT-based stuff also will, except for convolution, if a 
 zero-delay algo is used.

*/
