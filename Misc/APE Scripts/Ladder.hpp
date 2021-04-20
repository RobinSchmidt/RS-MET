#include "../../../../RS-MET/Misc/APE Scripts/rapt_for_ape.cpp"
// relative path from "Audio Programming Environment/includes"

#include <effect.h>

using namespace ape;

GlobalData(LadderFilter, "Moog-style Ladder Filter");

class LadderFilter : public ape::Effect
{
  
public:

  using LDR = RAPT::rsLadderFilter<float, float>;

  LadderFilter() {}
  

  static constexpr typename Param<LDR::modes>::Names ModeNames = 
  { 
	"Bypass",
    "Lowpass_6"
    "Lowpass_12"
    "Lowpass_18"
    "Lowpass_24"
    "Highpass_6",
    "Highpass_12",
    "Highpass_18",
    "Highpass_24"
  };
  

  // Shorthands for convenience to reduce boilerplate:
  using Par = ape::Param<float>;
  using Rng = ape::Range;
  using Map = Rng::Mapping;

  Par parFreq { "Freq",  Rng(20,  20000, Map::Exp) }; // in Hz
  Par parReso { "Reso",  Rng(0,   100)             }; // in %
  Par parDrive{ "Drive", Rng(-12, 48)              }; // in dB
  Par parGain { "Gain",  Rng(-48, 12)              }; // in dB
  
  // ToDo: Mode, Drive (pre-gain), Gain (post-gain), DC, SatShape
  
  // How can we set the parameters to reasonable default values?


private:

  LDR flt;
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
    const float freq  = (float)parFreq;
    const float reso  = (float)parReso * 0.01f;             // convert from percent to raw
    const float drive = RAPT::rsDbToAmp((float)parDrive);
    const float gain  = RAPT::rsDbToAmp((float)parGain);
    flt.setCutoff(freq);
    flt.setResonance(reso);
    // todo: set them per sample
    
    // Loop over the sample frames:
    const auto numChannels = sharedChannels();
    for(size_t n = 0; n < numFrames; ++n)
    {
      float x = drive * ins[0][n];            // read input - we just use the 1st channel
      float y = gain  * flt.getSample(x);     // compute output
      for(size_t c = 0; c < numChannels; ++c) // write generated sample into all channels
        outs[c][n] = y;                 
    }

    clear(outs, numChannels); // what does this do? clear unused channels?
  }
};

// Notes:
// -with high cutoff (> 12 kHz) and 100% resonance, it starts behaving chatocially
//  -> oversample by factor 2 (factor out a DSP class that wraps that)
// -with even more resonance (> 100%), the artifacts begin to show already at lower frequencies
// -When a filter should support resonances > 100, we need a nonlinear parameter mapping, 
//  maybe 100% should occur at 2/3 or 3/4 of the full range. actually, that would be good anyway
//  because from 0..50%, not much seems to be happening and before 100, a lot is happening. be 
//  sure to test this with settings where the drive does not diminish the resonance





