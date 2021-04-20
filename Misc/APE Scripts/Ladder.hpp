#include "../../../../RS-MET/Misc/APE Scripts/rapt_for_ape.cpp"
// relative path from "Audio Programming Environment/includes"

#include <effect.h>


GlobalData(LadderFilter, "Moog-style Ladder Filter");

class LadderFilter : public ape::Effect
{
  
public:

  LadderFilter() {}

  // Shorthands for convenience to reduce boilerplate:
  using Par = ape::Param<float>;
  using Rng = ape::Range;
  using Map = Rng::Mapping;

  Par parFreq{ "Freq", Rng(20, 20000, Map::Exp) };  // in Hz
  Par parReso{ "Reso", Rng(0,  100)             };  // in %
  // ToDo: Mode, Drive (pre-gain), Gain (post-gain), DC, SatShape
  
  // How can we set the parameters to reasonable default values?


private:

  RAPT::rsLadderFilter<float, float> flt;
  float sampleRate;

  
  /** Resets the internal state. */
  void reset()
  {
    flt.reset();
  }

  //-----------------------------------------------------------------------------------------------
  // \name Overriden callbacks (why are they private?) 

  void start(const ape::IOConfig& cfg) override
  {
    sampleRate = (float)cfg.sampleRate;
    flt.setSampleRate(sampleRate);
    reset();
  }

  void process(ape::umatrix<const float> ins, ape::umatrix<float> outs, size_t numFrames) override
  {
    // Pull the values of the parameters for this block and update the DSP objects:
    const float freq = (float)parFreq;
    const float reso = (float)parReso;
    flt.setCutoff(freq);
    flt.setResonance(reso);
    // todo: set them per sample
    
    // Loop over the sample frames:
    const auto numChannels = sharedChannels();
    for(size_t n = 0; n < numFrames; ++n)
    {
      float x = ins[0][n];                    // read input - we just use the 1st channel
      float y = flt.getSample(x);             // compute output
      for(size_t c = 0; c < numChannels; ++c) // write generated sample into all channels
        outs[c][n] = y;                 
    }

    clear(outs, numChannels); // what does this do? clear unused channels?
  }
};

// ToDo: factor out a class for the pure DSP process