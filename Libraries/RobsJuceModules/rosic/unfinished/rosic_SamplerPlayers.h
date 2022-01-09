#ifndef rosic_SamplerPlayers_h
#define rosic_SamplerPlayers_h

namespace rosic { namespace Sampler {


//===============================================================================================

class SignalProcessorChain
{
public:

  void processFrame(rsFloat64x2& inOut);
  void processBlock(rsFloat64x2* inOut, int N);
  void prepareToPlay(double fs) { for(auto & p : processors) p->prepareToPlay(fs); }
  //void resetState()    { for(auto & p : processors) p->resetState();    }
  //void resetSettings() { for(auto & p : processors) p->resetSettings(); }
  //void reset() { resetState(); resetSettings(); }

  void reserve(size_t num) { processors.reserve(num); }
  void addProcessor(SignalProcessor* p) { processors.push_back(p); }
  void clear() { processors.clear(); }

  bool isEmpty() const { return processors.empty(); }

  /** Returns the total number of processors in the chain. */
  size_t getNumProcessors() const { return processors.size(); }

  /** Returns the number of processors of given type in the chain. */
  size_t getNumProcessors(DspType type) const;


  SignalProcessor* getProcessor(int i) { return processors[i]; }

  /** Returns the index-th processor of the given type within the chain or nullptr, if there are
  not enough (i.e. less than i+1, counting starts at zero) processors of the given type in the 
  chain. To get the 3rd filter, you would pass type = SignalProcessorType::Filter, index = 2.  */
  SignalProcessor* getProcessor(DspType type, int index);

protected:

  std::vector<SignalProcessor*> processors;

};

//===============================================================================================




}}

#endif