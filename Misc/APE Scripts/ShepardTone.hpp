#include "../../../../RS-MET/Libraries/RobsJuceModules/rapt/rapt_for_ape.cpp"
// relative path from "Audio Programming Environment/includes"

#include <effect.h> // or maybe we should use generator.h?


//=================================================================================================
// The core DSP object

template<class T>
class rsShepardToneGenerator
{
  
public:
  

  /** Produces one output sample at a time. */
  inline T getSample()
  {
    return T(0);
  }
  
  
protected:




};


//=================================================================================================
// The APE effect


GlobalData(ShepardTone, "Endless glissando generator");

class ShepardToneAPE : public ape::Effect
{
public:


  ShepardTone() {}
  
  // Shorthands for convenience to reduce boilerplate:
  using Par = ape::Param<float>;
  using Rng = ape::Range;
  using Map = Rng::Mapping;

  Par parGain{   "Gain",     Rng(-48,  12)              }; // in dB
  Par parFreqLo{ "FreqLo",   Rng(20,   20000, Map::Exp) }; // in Hz
  Par parFreqHi{ "FreqHi",   Rng(20,   20000, Map::Exp) };  

private:   

  rsShepardToneGenerator<float> core;

  /** Resets the internal state. */
  void reset()
  {

  }

  //-----------------------------------------------------------------------------------------------
  // \name Overriden callbacks (why are they private?) 

  void start(const ape::IOConfig& cfg) override
  {
    reset();
  }

  void process(ape::umatrix<const float> ins, ape::umatrix<float> outs, size_t numFrames) override
  {
    // Pull the values of the parameters for this block and update the DSP object:
    const float gain = RAPT::rsDbToAmp((float)parGain);
    float s  = 2*PI / config().sampleRate;
    float wl = s * (float)parFreqLo;
    float wh = s * (float)parFreqHi;
    core.setOmegas(wl, wh);

    // Loop over the sample frames:
    const auto numChannels = sharedChannels();     // shared? with whom?
    for(size_t n = 0; n < numFrames; ++n)
    {
      float y = gain * core.getSample();      // compute output
      for(size_t c = 0; c < numChannels; ++c) // write generated sample into all channels
        outs[c][n] = y;                 
    }

    clear(outs, numChannels); // what does this do? clear unused channels?
  }
};

